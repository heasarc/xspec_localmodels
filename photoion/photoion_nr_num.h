/* Numerical Recipes numerical routines shared by all photoion models.
 * See photoion_nr_num.c. Arrays use the NR unit-offset convention.
 */

#ifndef PHOTOION_NR_NUM_H
#define PHOTOION_NR_NUM_H

double pion_qromb(double (*func)(double), double a, double b);
double pion_trapzd(double (*func)(double), double a, double b, int n);
void pion_polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void pion_spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void pion_splint(double xa[], double ya[], double y2a[], int n, double x, double *y);
void pion_gaussj(double **a, int n, double **b, int m);

#endif /* PHOTOION_NR_NUM_H */
