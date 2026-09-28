/* Numerical Recipes style array allocators shared by all photoion models.
 * See photoion_nr_alloc.c. All routines use the NR unit-offset convention
 * and zero the memory they return.
 */

#ifndef PHOTOION_NR_ALLOC_H
#define PHOTOION_NR_ALLOC_H

void pion_nrerror(char error_text[]);
double *pion_vector(int nl, int nh);
int *pion_ivector(int nl, int nh);
double *pion_dvector(int nl, int nh);
double **pion_matrix(int nrl, int nrh, int ncl, int nch);
double **pion_dmatrix(int nrl, int nrh, int ncl, int nch);
int **pion_imatrix(int nrl, int nrh, int ncl, int nch);
double **pion_submatrix(double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);
void pion_free_vector(double *v, int nl, int nh);
void pion_free_ivector(int *v, int nl, int nh);
void pion_free_dvector(double *v, int nl, int nh);
void pion_free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
void pion_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
void pion_free_imatrix(int **m, int nrl, int nrh, int ncl, int nch);
void pion_free_submatrix(double **b, int nrl, int nrh, int ncl, int nch);
double **pion_convert_matrix(double *a, int nrl, int nrh, int ncl, int nch);
void pion_free_convert_matrix(double **b, int nrl, int nrh, int ncl, int nch);

#endif /* PHOTOION_NR_ALLOC_H */
