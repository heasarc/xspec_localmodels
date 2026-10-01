/* Natural cubic spline interpolation, GSL-backed.
 *
 * Replaces the Numerical Recipes spline/splint pair. The signatures are
 * unchanged so that no call site had to move, but the meaning of the y2
 * argument has: it is now an opaque identity token, not an array of second
 * derivatives. See photoion_spline.c.
 */

#ifndef PHOTOION_SPLINE_H
#define PHOTOION_SPLINE_H

void pion_spline(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void pion_splint(double xa[], double ya[], double y2a[], int n, double x, double *y);

#endif /* PHOTOION_SPLINE_H */
