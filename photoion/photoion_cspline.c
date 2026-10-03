/* A compact natural cubic spline. See the header for why it exists.
 *
 * The mathematics. On each interval [x_i, x_{i+1}], of width h_i, the spline
 * is the cubic with values y_i, y_{i+1} and second derivatives M_i, M_{i+1}
 * at its ends. Continuity of the first derivative at the interior nodes
 * gives, for i = 1..n-2,
 *
 *   h_{i-1} M_{i-1} + 2 (h_{i-1} + h_i) M_i + h_i M_{i+1}
 *       = 6 [ (y_{i+1} - y_i) / h_i  -  (y_i - y_{i-1}) / h_{i-1} ],
 *
 * and a natural spline closes the system with M_0 = M_{n-1} = 0. The matrix
 * is symmetric, tridiagonal and strictly diagonally dominant, so forward
 * elimination followed by back substitution, without pivoting, is stable.
 *
 * Written relative to x_i, with d = x - x_i, the interval's cubic is
 *
 *   S(x) = y_i + d (b_i + d (M_i/2 + d (M_{i+1} - M_i) / (6 h_i))),
 *   b_i  = (y_{i+1} - y_i) / h_i  -  h_i (2 M_i + M_{i+1}) / 6,
 *
 * which also gives the continuation beyond either end when evaluated with
 * the first or last interval's coefficients.
 */
#include <stdlib.h>
#include <string.h>

#include "photoion_cspline.h"
#include "photoion_alloc.h"   /* pion_error */

/* Elimination workspace: two arrays of n, grown as needed and kept. Like the
 * package's other caches it is allocated once and deliberately not freed. */
static double *work = NULL;
static int work_n = 0;

static int ensure_work(int n)
{
  double *w;
  if (n <= work_n) return 0;
  w = realloc(work, 2 * (size_t) n * sizeof *w);
  if (w == NULL) return 1;
  work = w;
  work_n = n;
  return 0;
}

size_t pion_cspline_bytes(int n)
{
  return 2 * (size_t) n * sizeof(double);
}

void pion_cspline_free(struct pion_cspline *s)
{
  free(s->y);
  free(s->M);
  s->y = s->M = NULL;
  s->x = NULL;
  s->n = 0;
  s->hint = 0;
}

int pion_cspline_init(struct pion_cspline *s, const double *x, const double *y, int n)
{
  int i;
  double *cp, *rp, a, b, c, r, den;

  s->x = NULL; s->y = s->M = NULL; s->n = 0; s->hint = 0;
  if (n < 2) { pion_error("pion_cspline_init: fewer than 2 nodes"); return 1; }
  for (i = 1; i < n; ++i)
    if (!(x[i] > x[i-1])) { pion_error("pion_cspline_init: abscissae not increasing"); return 1; }

  s->y = malloc((size_t) n * sizeof *s->y);
  s->M = malloc((size_t) n * sizeof *s->M);
  if (s->y == NULL || s->M == NULL || ensure_work(n)) {
    free(s->y); free(s->M); s->y = s->M = NULL;
    pion_error("pion_cspline_init: out of memory");
    return 1;
  }
  memcpy(s->y, y, (size_t) n * sizeof *y);
  s->x = x;
  s->n = n;
  s->M[0] = s->M[n-1] = 0.;
  if (n == 2) return 0;

  /* Forward elimination over the unknowns M_1..M_{n-2}. cp holds the
   * normalized superdiagonal, rp the normalized right-hand side. */
  cp = work;
  rp = work + n;
  for (i = 1; i <= n-2; ++i) {
    a = x[i] - x[i-1];                 /* h_{i-1}: subdiagonal */
    c = x[i+1] - x[i];                 /* h_i:     superdiagonal */
    b = 2. * (a + c);                  /* diagonal */
    r = 6. * ((y[i+1] - y[i]) / c - (y[i] - y[i-1]) / a);
    if (i == 1) {
      den = b;
      rp[i] = r / den;
    } else {
      den = b - a * cp[i-1];
      rp[i] = (r - a * rp[i-1]) / den;
    }
    cp[i] = c / den;
  }
  /* Back substitution; M_{n-1} = 0 closes it. */
  s->M[n-2] = rp[n-2];
  for (i = n-3; i >= 1; --i) s->M[i] = rp[i] - cp[i] * s->M[i+1];
  return 0;
}

/* The interval i, 0 <= i <= n-2, with x[i] <= x < x[i+1]; the first interval
 * below the table and the last at or above its end. */
static int interval(struct pion_cspline *s, double x)
{
  const double *xa = s->x;
  int lo, hi, mid, i = s->hint;

  if (x < xa[1]) return 0;
  if (x >= xa[s->n-2]) return s->n-2;
  if (xa[i] <= x && x < xa[i+1]) return i;
  if (i+2 < s->n && xa[i+1] <= x && x < xa[i+2]) return i+1;
  lo = 0; hi = s->n-1;               /* invariant: xa[lo] <= x < xa[hi] */
  while (hi - lo > 1) {
    mid = (lo + hi) / 2;
    if (xa[mid] <= x) lo = mid; else hi = mid;
  }
  return lo;
}

double pion_cspline_eval(struct pion_cspline *s, double x)
{
  int i = interval(s, x);
  const double *xa = s->x, *y = s->y, *M = s->M;
  double h = xa[i+1] - xa[i];
  double d = x - xa[i];
  double b = (y[i+1] - y[i]) / h - h * (2. * M[i] + M[i+1]) / 6.;

  s->hint = i;
  return y[i] + d * (b + d * (0.5 * M[i] + d * (M[i+1] - M[i]) / (6. * h)));
}
