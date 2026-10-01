/* Unit-offset array helpers shared by all photoion models.
 *
 * Written fresh for this package; see photoion_alloc.c. A vector allocated
 * [nl..nh] is indexed v[nl]..v[nh], and a matrix [nrl..nrh][ncl..nch] is
 * indexed m[nrl][ncl].. m[nrh][nch]. The memory comes back zeroed.
 *
 * Requires 0 <= nl <= nh. Every call site in the package passes nl = 1.
 */

#ifndef PHOTOION_ALLOC_H
#define PHOTOION_ALLOC_H

void pion_error(const char *msg);

int    *pion_ivector(int nl, int nh);
double *pion_dvector(int nl, int nh);
double **pion_dmatrix(int nrl, int nrh, int ncl, int nch);

void pion_free_ivector(int *v, int nl, int nh);
void pion_free_dvector(double *v, int nl, int nh);
void pion_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

#endif /* PHOTOION_ALLOC_H */
