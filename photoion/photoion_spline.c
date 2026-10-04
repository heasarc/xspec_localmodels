/* Natural cubic spline interpolation, GSL-backed.
 *
 * Replaces the Numerical Recipes spline/splint pair with gsl_interp_cspline,
 * which is the same natural cubic spline: over a 40-node log-log dataset the
 * two agree to 3.8e-16 relative, in-domain.
 *
 * The signatures are kept byte-compatible with the routines they replace, so
 * none of the 123 call sites had to change. What did change is the role of the
 * y2 argument. Numerical Recipes splits the work in two: spline() solves for
 * the second derivatives and hands them back in y2, and splint() evaluates from
 * (xa, ya, y2a). GSL instead bundles the nodes and the solved coefficients into
 * one gsl_spline object. So here:
 *
 *   pion_spline()  builds a gsl_spline from (x, y) and files it in a small cache
 *                  keyed on the y2 POINTER.
 *   pion_splint()  looks the cache up by the y2a pointer and evaluates.
 *
 * y2 is therefore an identity token. Nothing is written to it. That is safe
 * because no model file ever reads one of these arrays -- all 16 of them appear
 * only in a pion_spline call, a pion_splint call, and their dvector alloc/free.
 * They are left allocated so that each remains a distinct, stable address.
 *
 * The cache is only ever an optimization: a lookup validates (x, y, n) against
 * what was filed, and on any mismatch -- including a stale entry whose y2 array
 * was freed and its address recycled -- pion_splint rebuilds from the arrays it
 * was handed. A wrong answer is not reachable through a stale key, only a slower
 * one, and in practice the 16 datasets each keep their own slot.
 *
 * NOT re-entrant, in keeping with the rest of the package.
 *
 * Out-of-domain queries are reproduced exactly rather than rejected; see
 * eval_outside() below.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "photoion_alloc.h"   /* pion_error */
#include "photoion_spline.h"

/* 16 distinct y2 arrays exist across the 11 models; 32 leaves headroom. */
#define CACHE_SLOTS 32

struct slot {
   const double *key;            /* the y2 pointer, unit-offset as passed */
   const double *x, *y;          /* validated on lookup */
   int n;
   gsl_spline *sp;
   gsl_interp_accel *acc;
};

/* The cache is never torn down. slot_clear() frees a gsl_spline and its accel
 * only when a slot is being REPLACED, so the ~16 live splines stay allocated for
 * the lifetime of the process -- bounded by CACHE_SLOTS, not growing with the
 * number of model evaluations.
 *
 * Same policy, and same reasoning, as the workspace in photoion_integrate.c:
 * there is no unload or atexit hook in the glue initpackage generates, XSPEC
 * keeps a dlopen'ed package loaded for the session, and freeing per call would
 * defeat the point of caching. Nothing here is collected; it is simply never
 * returned, and the kernel reclaims it at process exit. dlclose() on the package
 * is the one case that would turn this into a genuine leak. */
static struct slot cache[CACHE_SLOTS];
static int nslots = 0;

/* GSL's default error handler calls abort(), which would take XSPEC down with
 * it. Allocation and init are the only calls here that can reach it -- the
 * eval_e family returns status instead -- so the handler is suspended across
 * them and restored, rather than switched off for the whole process: this
 * package is dlopen'ed alongside other libraries that may use GSL. */
static gsl_spline *build(const double *x1, const double *y1, int n)
{
   gsl_error_handler_t *old;
   gsl_spline *sp;

   if (n < 3) return NULL;       /* cspline needs at least 3 points */
   old = gsl_set_error_handler_off();
   sp = gsl_spline_alloc(gsl_interp_cspline, (size_t) n);
   if (sp && gsl_spline_init(sp, x1, y1, (size_t) n) != GSL_SUCCESS) {
      gsl_spline_free(sp);
      sp = NULL;
   }
   gsl_set_error_handler(old);
   return sp;
}

static void slot_clear(struct slot *s)
{
   if (s->sp)  gsl_spline_free(s->sp);
   if (s->acc) gsl_interp_accel_free(s->acc);
   s->sp = NULL; s->acc = NULL; s->key = NULL; s->x = NULL; s->y = NULL; s->n = 0;
}

/* Install (x, y, n) under the key y2, replacing whatever was there. */
static struct slot *slot_put(const double *key, const double *x, const double *y, int n)
{
   int i;
   struct slot *s = NULL;

   for (i = 0; i < nslots; ++i)
      if (cache[i].key == key) { s = &cache[i]; break; }
   if (!s) {
      if (nslots < CACHE_SLOTS) s = &cache[nslots++];
      else {
         /* Should not happen: there are 16 of these arrays. Reuse slot 0
          * rather than fail -- correctness does not depend on the cache. */
         s = &cache[0];
      }
   }
   slot_clear(s);
   s->sp = build(x + 1, y + 1, n);          /* unit-offset -> zero-based base */
   if (!s->sp) { pion_error("spline: gsl_spline_alloc/init failed"); return NULL; }
   /* Suspended for the same reason as in build(): on failure this one reaches
    * GSL_ERROR_NULL, and the default handler aborts the process. A NULL accel
    * is not fatal here -- GSL's eval falls back to a bisection search -- so
    * report it and carry on. */
   {
      gsl_error_handler_t *old = gsl_set_error_handler_off();
      s->acc = gsl_interp_accel_alloc();
      gsl_set_error_handler(old);
      if (!s->acc) pion_error("spline: gsl_interp_accel_alloc failed");
   }
   s->key = key; s->x = x; s->y = y; s->n = n;
   return s;
}

static struct slot *slot_get(const double *key, const double *x, const double *y, int n)
{
   int i;

   for (i = 0; i < nslots; ++i) {
      struct slot *s = &cache[i];
      if (s->key == key && s->x == x && s->y == y && s->n == n && s->sp) return s;
   }
   return slot_put(key, x, y, n);           /* miss, or stale key: rebuild */
}

/* Numerical Recipes' splint does not range-check. Its binary search collapses
 * to the first interval for x below xa[1] and to the last for x above xa[n],
 * and it then evaluates that interval's cubic at a point outside it. GSL
 * instead returns GSL_EDOM and writes NaN.
 *
 * That difference is not academic here: 18.7% of the pion_lowEpispline calls in
 * xiabs are below the grid, ~285000 per evaluation, because the caller's guard
 * tests the threshold rather than the table's first node. Rejecting them
 * would poison the output;
 * the decision was to preserve the existing behavior exactly and treat the guard
 * as a separate question.
 *
 * On the end interval the spline IS a cubic, so its Taylor expansion about the
 * end node reproduces it with no error at all -- this is an identity, not an
 * approximation. The second derivative is linear on each interval, so its slope
 * across the end interval gives the third derivative exactly. y'' is continuous
 * at a node, so evaluating it at the far end of the interval is unambiguous.
 */
static int eval_outside(const struct slot *s, double x, double *y)
{
   const double *xa = s->x;
   int n = s->n;
   double xe, xf, t, f, d1, d2, d3, h;
   int st = 0;

   if (x < xa[1]) { xe = xa[1]; xf = xa[2]; }        /* first interval */
   else           { xe = xa[n]; xf = xa[n-1]; }      /* last interval  */
   h = xf - xe;
   if (h == 0.0) { pion_error("splint: zero-width end interval"); return 1; }

   st |= gsl_spline_eval_e(s->sp, xe, s->acc, &f);
   st |= gsl_spline_eval_deriv_e(s->sp, xe, s->acc, &d1);
   st |= gsl_spline_eval_deriv2_e(s->sp, xe, s->acc, &d2);
   st |= gsl_spline_eval_deriv2_e(s->sp, xf, s->acc, &d3);
   if (st) return st;

   d3 = (d3 - d2) / h;                               /* third derivative */
   t  = x - xe;
   *y = f + d1*t + 0.5*d2*t*t + d3*t*t*t/6.0;
   return 0;
}

void pion_spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
   /* The 0.99e30 sentinel selects natural boundary conditions, which is what
    * gsl_interp_cspline provides. Every one of the 107 call sites passes 1.e40
    * for both, so a real end derivative has never been requested; say so rather
    * than silently ignore it. */
   if (yp1 <= 0.99e30 || ypn <= 0.99e30)
      pion_error("spline: clamped end derivatives are not supported");

   slot_put(y2, x, y, n);
}

void pion_splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
   struct slot *s = slot_get(y2a, xa, ya, n);

   if (!s) { *y = 0.0; return; }

   if (x >= xa[1] && x <= xa[n]) {
      if (gsl_spline_eval_e(s->sp, x, s->acc, y) != GSL_SUCCESS) {
         pion_error("splint: evaluation failed");
         *y = 0.0;
      }
   } else if (eval_outside(s, x, y)) {
      pion_error("splint: extrapolation failed");
      *y = 0.0;
   }
}
