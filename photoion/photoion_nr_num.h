/* Numerical Recipes numerical routines shared by all photoion models.
 * See photoion_nr_num.c. Arrays use the NR unit-offset convention.
 */

#ifndef PHOTOION_NR_NUM_H
#define PHOTOION_NR_NUM_H

double pion_qromb(double (*func)(double), double a, double b);
double pion_trapzd(double (*func)(double), double a, double b, int n);
void pion_polint(double xa[], double ya[], int n, double x, double *y, double *dy);

#endif /* PHOTOION_NR_NUM_H */
