/* Numerical Recipes numerical routines shared by all photoion models.
 *
 * Extracted verbatim from the per-model copies, which were byte-identical
 * apart from a name suffix. Romberg integration (qromb/trapzd), polynomial
 * interpolation (polint) and cubic splines (spline/splint).
 *
 * gaussj (Gauss-Jordan elimination) was carried over too, but had no callers
 * in any of the 11 models, and is gone. trapzd and polint are reached only
 * from qromb.
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

void pion_spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
   int i,k;
   double p, qn, sig, un, *u;

   u=pion_vector(1,n-1);
   if (yp1 > 0.99e30)
      y2[1]=u[1]=0.0;
   else {
      y2[1] = -0.5;
      u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
   }
   for (i=2;i<=n-1;i++) {
      sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
      p=sig*y2[i-1]+2.0;
      y2[i]=(sig-1.0)/p;
      u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
      u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
   }
   if (ypn > 0.99e30)
      qn=un=0.0;
   else {
      qn=0.5;
      un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
   }
   y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
   for (k=n-1;k>=1;k--)
      y2[k]=y2[k]*y2[k+1]+u[k];
   pion_free_vector(u,1,n-1);
}

void pion_splint(double xa[], double ya[], double y2a[], int n, double x, double *y)
{
   int klo,khi,k;
   double h,b,a;

   klo=1;
   khi=n;
   while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (xa[k] > x) khi=k;
      else klo=k;
   }
   h=xa[khi]-xa[klo];
   if (h == 0.0) pion_nrerror("Bad XA input to routine SPLINT");
   a=(xa[khi]-x)/h;
   b=(x-xa[klo])/h;
   *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}
