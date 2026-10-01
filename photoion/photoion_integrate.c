/* Adaptive numerical integration, GSL-backed.
 *
 * Replaces the Numerical Recipes qromb (with trapzd and polint, which were
 * reachable only from it) by gsl_integration_qag. All 19 live call sites moved
 * in one commit, because NR's trapzd carried
 *
 *     static double s; static int it;
 *
 * across calls: those two statics are the running trapezoid sum and the current
 * panel count, so any call whose integrand itself integrated would corrupt
 * both. Nothing in the package nests today, but a half-converted tree would
 * have had two integrators sharing that hazard asymmetrically, which is harder
 * to reason about than either end state.
 *
 * WHY qag AND NOT gsl_integration_romberg
 *
 * GSL 2.7 does ship a Romberg routine, and being the algorithmic look-alike it
 * reads as the conservative swap. Measurement says otherwise: on the same
 * integrands gsl_integration_romberg and qag agree with EACH OTHER to ~1e-11
 * while both differ from NR's qromb by ~1e-6. NR is the outlier, not GSL.
 *
 * Its convergence test is the reason. qromb extrapolates the trapezoid sequence
 * with polint and stops when the extrapolation's own error estimate dss
 * satisfies |dss| < EPS*|ss|. That estimates the error of the extrapolation,
 * not of the result, and it stops optimistically early. Where a closed form is
 * available the direction is unambiguous: for E*L(E) with L a power law over
 * 0.05-20 keV the integral is exactly ln(400) = 5.99146454710798, which qag
 * returns to 1.5e-16 while qromb gives 5.99162461529 -- off by 2.7e-5, having
 * been asked for 1e-5.
 *
 * So the output moving by ~1e-5 at these call sites is NR's error being
 * removed, not error being introduced.
 *
 * qag is also far cheaper, because it refines where the integrand needs it
 * instead of doubling a uniform grid. Evaluations to convergence:
 *
 *     integrand     NR qromb   GSL romberg   GSL qag
 *     smooth              33           129        21
 *     peaked           8 193         4 097       483
 *     power law        4 097       262 145       357
 *
 * The extreme case is in photoion itself, where qromb(pion_integrand, EMIN,
 * EMAX) exhausts JMAX = 22 without converging: 2 097 152 integrand
 * evaluations, each one a spline evaluation, for a result that reaches nothing
 * but the "Radiation Pressure" printf.
 *
 * WHY GSL_INTEG_GAUSS61
 *
 * The 21-point rule is the usual recommendation and was the one this migration
 * was planned around. It is not safe here. An adaptive rule can only see
 * structure its nodes land on, and on a test integrand of 20 Gaussian spikes of
 * width 0.002 over [0.3, 10] -- deliberately shaped like a line spectrum --
 * GAUSS21 returns 0.0097 against a true value of ~0.110 and reports
 * GSL_SUCCESS. It stepped over every spike and converged confidently on the
 * continuum alone. GAUSS15 and GAUSS31 also fail (0.025 and 0.055); GAUSS61
 * recovers 0.10997 for 10 553 evaluations.
 *
 * That shape is not hypothetical: pion_integrand interpolates the emitted
 * spectrum, emission lines included, across the whole band. A silently wrong
 * answer is the one failure mode worth paying for, and the price here is a few
 * thousand evaluations against the 2 097 152 being removed.
 *
 * ERROR HANDLING
 *
 * GSL's default handler calls abort(), which would take XSPEC down with it, and
 * the handler is process-global while this package is dlopen'ed alongside other
 * libraries that may use GSL. So it is suspended across the call and restored,
 * never switched off for good.
 *
 * On non-convergence qag still writes its best estimate, and that is returned
 * with a report -- the same contract as qromb, which fell through to its error
 * stub (whose body was commented out) and returned its last unconverged value.
 *
 * NOT re-entrant: one workspace is shared by every call, in keeping with the
 * rest of the package. None of the four live integrands integrates.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>

#include "photoion_alloc.h"   /* pion_error */
#include "photoion_integrate.h"

/* Subinterval limit. The workspace is sized to match, so qag can bisect 1000
 * times before giving up -- up to 61 000 evaluations at GAUSS61, which is still
 * two orders of magnitude below what the uniform refinement it replaces spent
 * on photoion's worst integrand. */
#define LIMIT   (1000)

/* Relative target. NR asked for EPS = 1e-5 and, as above, did not always
 * deliver it; 1e-8 is comfortably better while staying far from the roundoff
 * floor where qag starts reporting GSL_EROUND. epsabs is 0: every integrand
 * here is non-negative, so there is no cancellation for a relative-only
 * criterion to mishandle, and a genuinely zero integral returns 0 cleanly. */
#define EPSREL  (1.0e-8)
#define EPSABS  (0.0)

#define KEY     GSL_INTEG_GAUSS61

/* Allocated on first use and then DELIBERATELY NEVER FREED. There is no
 * pion_integrate_free(), and that is not an oversight.
 *
 * It is one bounded block for the lifetime of the process: at LIMIT = 1000 the
 * workspace is four double arrays plus two size_t arrays, 48 088 bytes measured
 * on this target, allocated once no matter how many times XSPEC evaluates a
 * model. Freeing it per call would defeat the only reason it is static, since
 * xiabs makes hundreds of pion_integrate calls per evaluation.
 *
 * There is also nowhere to free it from. The glue initpackage generates has no
 * unload, finalize or atexit hook, and XSPEC keeps a dlopen'ed model package
 * loaded for the whole session, so process exit is the only teardown point --
 * and there the kernel reclaims the address space regardless.
 *
 * Two things that are easy to assume here and are both wrong. The pointer below
 * has static storage duration, so it never "goes out of scope" during the run;
 * and C frees nothing on its own -- the block is not collected, it is simply
 * never returned. `static` buys reuse, not reclamation: a malloc leaked from an
 * ordinary local would behave the same way at exit, only unreachable in the
 * meantime.
 *
 * The one case where this would become a real leak is dlclose() on the package:
 * the pointer's storage would be unmapped while the block stayed allocated and
 * permanently unreachable. XSPEC does not unload model packages within a
 * session. If that ever changes, this needs an explicit teardown -- and so does
 * the spline cache in photoion_spline.c, which holds the same policy. */
static gsl_integration_workspace *ws = NULL;

/* gsl_function passes a void* through to the integrand, so the bare
 * double(*)(double) the call sites use is carried in a struct rather than cast
 * to void* and back -- that cast is not something C guarantees for function
 * pointers. The struct is a local of pion_integrate, so this adapter adds no
 * state of its own. */
struct wrapped {
   double (*f)(double);
};

static double adapter(double x, void *params)
{
   return ((struct wrapped *) params)->f(x);
}

double pion_integrate(double (*func)(double), double a, double b)
{
   struct wrapped w;
   gsl_function F;
   gsl_error_handler_t *old;
   double result = 0.0, abserr;
   int status;

   if (a == b) return 0.0;      /* qag agrees, but say so without the setup */

   w.f = func;
   F.function = adapter;
   F.params   = &w;

   old = gsl_set_error_handler_off();

   if (!ws) ws = gsl_integration_workspace_alloc((size_t) LIMIT);
   if (!ws) {
      gsl_set_error_handler(old);
      pion_error("integrate: gsl_integration_workspace_alloc failed");
      return 0.0;
   }

   status = gsl_integration_qag(&F, a, b, EPSABS, EPSREL, (size_t) LIMIT,
                                KEY, ws, &result, &abserr);

   gsl_set_error_handler(old);

   /* Best estimate plus a report, as qromb did on JMAX exhaustion. */
   if (status != GSL_SUCCESS)
      pion_error("integrate: gsl_integration_qag did not converge");

   return result;
}
