/* Numerical Recipes style array allocators shared by all photoion models.
 * See photoion_nr_alloc.c. All routines use the NR unit-offset convention
 * and zero the memory they return.
 *
 * Only the allocators the models actually call are here; note that vector and
 * dvector are identical, both having been widened from NR's float.
 */

#ifndef PHOTOION_NR_ALLOC_H
#define PHOTOION_NR_ALLOC_H

void pion_nrerror(char error_text[]);
int *pion_ivector(int nl, int nh);
double *pion_dvector(int nl, int nh);
double **pion_dmatrix(int nrl, int nrh, int ncl, int nch);
void pion_free_ivector(int *v, int nl, int nh);
void pion_free_dvector(double *v, int nl, int nh);
void pion_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);

#endif /* PHOTOION_NR_ALLOC_H */
