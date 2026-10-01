/* Unit-offset array helpers shared by all photoion models.
 *
 * Written fresh for this package. These replace the Numerical Recipes nrutil
 * allocators that had been copied into every model file; with them gone,
 * nothing in photoion derives from Numerical Recipes any more.
 *
 * What is kept is the unit-offset CONVENTION -- v[nl] is the first element --
 * because ~12 000 lines of model code loop `for (k=1;k<=SPECBINS;++k)`. A
 * convention is an idea; the licence covers expression.
 *
 * Deliberately not gsl_vector/gsl_matrix. Those are zero-based, so adopting
 * them would force rebasing every index in the package: the single largest
 * off-by-one risk available here, for no benefit, since these are allocation
 * helpers and not algorithms.
 *
 * Three differences from the routines they replace, all intentional:
 *
 * ZEROING, which is behavior the models depend on. photoion.c's variants zeroed
 * and the other ten files' did not, and several models accumulate with += onto
 * a freshly allocated array; sharing the zeroing version fixed that class of
 * bug library-wide, and this preserves it. calloc() does the work rather than a
 * store loop, which also lets the allocator hand back pre-zeroed pages for the
 * large SPECBINS arrays instead of dirtying them. That relies on all-bits-zero
 * being the zero value: guaranteed by C for the integer types, and by IEEE-754
 * for double, hence true on every platform HEASOFT builds on.
 *
 * NO POINTER UNDERFLOW. Numerical Recipes offsets by returning `v - nl`, which
 * for nl > 0 forms a pointer before the start of the allocation -- undefined
 * behavior, however reliably it works in practice, and the construct
 * -Wfree-nonheap-object exists to flag when the matching free() unwinds it.
 * These allocate indices 0..nh instead and return the base, so every pointer
 * formed stays inside its allocation and free() receives exactly the pointer
 * malloc returned. The cost is nl unused elements per array, and all 366
 * allocation call sites in the package use lower bound 1 -- 306 dvector, 10
 * ivector, 50 dmatrix with nrl = ncl = 1 -- so that is one element. It does
 * mean 0 <= nl is
 * required, and that a large nl would over-allocate rather than offset -- so
 * the helpers reject out-of-range bounds instead of quietly obliging.
 *
 * A REAL ERROR PATH. pion_error, which these used to call, had a commented-out
 * body: an allocation failure was silent and the caller went on to dereference
 * NULL. pion_error now reports. The callers still do not check -- that would
 * mean touching all 366 sites -- so a failure is still fatal, but it now says
 * what failed first instead of only producing a segfault.
 *
 * pion_dmatrix additionally unwinds a partial failure: the old version left
 * already-allocated rows leaked and NULL rows in the array for the caller to
 * walk into.
 */

#include <stdlib.h>
#include <stdio.h>

#include "photoion_alloc.h"

/* Reported on stderr rather than through XSPEC's own channel: the models all
 * use plain printf for their output, so there is no established alternative
 * here, and errors belong on stderr regardless. */
void pion_error(const char *msg)
{
   fprintf(stderr, "photoion: %s\n", msg ? msg : "(unspecified error)");
   fflush(stderr);
}

/* Shared bound check. Returns the element count to allocate (nh + 1, covering
 * the unused 0..nl-1 slots), or 0 if the bounds are unusable. */
static size_t span(int nl, int nh, const char *what)
{
   if (nl < 0 || nh < nl) {
      pion_error(what);
      return 0;
   }
   return (size_t) nh + 1;
}

int *pion_ivector(int nl, int nh)
{
   size_t n = span(nl, nh, "ivector: bad bounds");
   int *v;

   if (n == 0) return NULL;
   v = calloc(n, sizeof *v);
   if (!v) pion_error("ivector: out of memory");
   return v;
}

double *pion_dvector(int nl, int nh)
{
   size_t n = span(nl, nh, "dvector: bad bounds");
   double *v;

   if (n == 0) return NULL;
   v = calloc(n, sizeof *v);
   if (!v) pion_error("dvector: out of memory");
   return v;
}

double **pion_dmatrix(int nrl, int nrh, int ncl, int nch)
{
   size_t nrow = span(nrl, nrh, "dmatrix: bad row bounds");
   size_t ncol = span(ncl, nch, "dmatrix: bad column bounds");
   double **m;
   int i;

   if (nrow == 0 || ncol == 0) return NULL;

   m = calloc(nrow, sizeof *m);
   if (!m) { pion_error("dmatrix: out of memory"); return NULL; }

   for (i = nrl; i <= nrh; ++i) {
      m[i] = calloc(ncol, sizeof **m);
      if (!m[i]) {                 /* unwind rather than hand back holes */
         while (--i >= nrl) free(m[i]);
         free(m);
         pion_error("dmatrix: out of memory");
         return NULL;
      }
   }
   return m;
}

/* nl and nh are unused -- the returned pointer is the one malloc gave us. They
 * stay in the signature so that none of the ~300 call sites had to change. */
void pion_free_ivector(int *v, int nl, int nh)
{
   (void) nl; (void) nh;
   free(v);
}

void pion_free_dvector(double *v, int nl, int nh)
{
   (void) nl; (void) nh;
   free(v);
}

void pion_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   (void) ncl; (void) nch;
   if (!m) return;
   for (i = nrh; i >= nrl; --i) free(m[i]);
   free(m);
}
