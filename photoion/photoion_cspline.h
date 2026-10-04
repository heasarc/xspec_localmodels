/* A compact natural cubic spline for the per-edge photoionization tables.
 *
 * Written from the defining equations (see photoion_cspline.c), not adapted
 * from any library. It exists because those tables are cached (one per FAC
 * record, up to a few thousand), so each must be small. It stores only the
 * values y and the second derivatives M; the abscissae are shared and
 * referenced, not copied. All other splines in the package stay on GSL
 * through photoion_spline.c.
 *
 * Indexing is 0-based, n >= 2. Not re-entrant, like the rest of the package:
 * the interval hint is updated on every evaluation.
 */
#ifndef PHOTOION_CSPLINE_H
#define PHOTOION_CSPLINE_H

#include <stddef.h>

struct pion_cspline {
  const double *x;   /* n abscissae, strictly increasing; not owned */
  double *y;         /* n values, owned */
  double *M;         /* n second derivatives, owned; M[0] = M[n-1] = 0 */
  int n;
  int hint;          /* last interval used, to start the search there */
};

/* Fit a natural spline through (x[i], y[i]), i = 0..n-1. Copies y; keeps the
 * pointer x, which must stay valid and unchanged while the spline is used.
 * Returns 0 on success, nonzero on bad input or allocation failure (and then
 * leaves the spline empty). */
int pion_cspline_init(struct pion_cspline *s, const double *x, const double *y, int n);

/* Release y and M and empty the spline. Safe on an empty spline. */
void pion_cspline_free(struct pion_cspline *s);

/* The spline at x. Outside [x[0], x[n-1]] the end interval's cubic is
 * continued, which is what the package has always done. Callers do evaluate
 * below x[0]: the edge tables between threshold and their first node. */
double pion_cspline_eval(struct pion_cspline *s, double x);

/* Bytes owned by a spline of n nodes, for cache accounting. */
size_t pion_cspline_bytes(int n);

#endif
