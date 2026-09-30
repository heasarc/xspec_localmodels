/* Numerical Recipes numerical routines shared by all photoion models.
 *
 * Extracted verbatim from the per-model copies, which were byte-identical
 * apart from a name suffix. Romberg integration (qromb/trapzd) and polynomial
 * interpolation (polint).
 *
 * gaussj (Gauss-Jordan elimination) was carried over too, but had no callers
 * in any of the 11 models, and is gone. trapzd and polint are reached only
 * from qromb. spline/splint have been replaced by gsl_interp_cspline; see
 * photoion_spline.c.
 *
 * Arrays follow the Numerical Recipes unit-offset convention; see
 * photoion_nr_alloc.h.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "photoion_nr_alloc.h"
#include "photoion_nr_num.h"

#define EPS   (1.0e-5)
#define JMAX  (22)
#define JMAXP (JMAX+1)
#define K     (5)

#define pion_FUNC(x) ((*func)(x))

double pion_qromb(double (*func)(double), double a, double b)
{
   double ss, dss;
   double s[JMAXP+1],h[JMAXP+1];
   int j;

   h[1]=1.0;
   for (j=1;j<=JMAX;j++) {
      s[j]=pion_trapzd(func,a,b,j);
      if (j >= K) {
         pion_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
         if (fabs(dss) < EPS*fabs(ss)) return ss;
      }
      /*printf("%2d     %e\n",j,ss);*/
      s[j+1]=s[j];
      h[j+1]=0.25*h[j];
   }
   pion_nrerror("Too many steps in routine QROMB");
   /*     printf("not accurate\n");*/
   return ss;
}

double pion_trapzd(double (*func)(double), double a, double b, int n)
{
   double x,tnm,sum,del;
   static double s;
   static int it;
   int j;

   if (n == 1) {
      it=1;
      return (s=0.5*(b-a)*(pion_FUNC(a)+pion_FUNC(b)));
   } else {
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += pion_FUNC(x);
      it *= 2;
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}

void pion_polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   int i,m,ns=1;
   double den,dif,dift,ho,hp,w;
   double *c, *d;

   dif=fabs(x-xa[1]);
   c=pion_vector(1,n);
   d=pion_vector(1,n);
   for (i=1;i<=n;i++) {
      if ( (dift=fabs(x-xa[i])) < dif) {
         ns=i;
         dif=dift;
      }
      c[i]=ya[i];
      d[i]=ya[i];
   }
   *y=ya[ns--];
   for (m=1;m<n;m++) {
      for (i=1;i<=n-m;i++) {
         ho=xa[i]-x;
         hp=xa[i+m]-x;
         w=c[i+1]-d[i];
         if ( (den=ho-hp) == 0.0) pion_nrerror("Error in routine POLINT");
         den=w/den;
         d[i]=hp*den;
         c[i]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   pion_free_vector(d,1,n);
   pion_free_vector(c,1,n);
}
