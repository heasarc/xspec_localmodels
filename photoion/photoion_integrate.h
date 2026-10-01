/* Adaptive numerical integration, GSL-backed.
 *
 * Replaces the Numerical Recipes qromb/trapzd/polint trio. The signature
 * matches the qromb it replaces, so no call site changed shape; only the name
 * did, since the method is no longer Romberg. See photoion_integrate.c.
 */

#ifndef PHOTOION_INTEGRATE_H
#define PHOTOION_INTEGRATE_H

double pion_integrate(double (*func)(double), double a, double b);

#endif /* PHOTOION_INTEGRATE_H */
