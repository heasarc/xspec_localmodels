#include "cfortran.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int neutral
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(neutral,NEUTRAL,neutral,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))
#define SMALL (1.e-6)
#define FINE_STRUCTURE (1./137.0359895)
#define eVtoHartree (1./(2.*13.6056981))
#define Rydberg (13.6056981)        /* 1 Rydberg [eV] */
#define a0 (5.29177249e-9)          /* Bohr radius [cm] */
#define ccc (2.99792458e10)         /* speed of light [cm/s] */
#define eee (4.8032068e-10)         /* electron charge [esu] */
#define hhh (6.6260755e-27)         /* Planck's constant [cgs] */
#define H_0 (2.301e-18) /* Hubble constant [1/s]: (WMAP h = 71) 71*1e5/3.085678e18/1e6*/
#define me (9.1093897e-28)          /* electron mass [g] */
#define re (2.81794092e-13)         /* classical electron radius: e^2/m c^2 */
#define ge (2.)                     /* gyromagnetic ratio for the electron */
#define meeV (5.1099906e5)          /* electron mass [eV] */
#define sigmaT (6.6525e-25)    /* Thomson cross-section [cm^2]: 8*Pi/3*re^2 */
#define eVtoergs (1.60217733e-12)   /* convert eV to ergs */
#define ergstoeV (1./eVtoergs)      /* convert ergs to eV */
#define AMconv (4.13413733134e16)   /* convert Einstein A_ij A.U.'s -> s^-1 */
#define lambdatokeV (12.3984244)    /* lambda=12.3984244/E_keV */
#define FACTOR (66.784)     /* For calculating line center optical depth */
#define parsectocm (3.085678e18)  /* parsecs to cm */
#define TEMPERATURES (30)
#define SMALL (1.e-6)
#define EPS 1.0e-5
#define JMAX 22
#define JMAXP JMAX+1
#define K 5
#define FUNC_ne(x) ((*func)(x))
#define SWAP_ne(a,b) {double temp=(a);(a)=(b);(b)=temp;}
#define PI (3.141592653589793)

/* START */

/* Forward declarations (converted from K&R style). */
void nrerror_ne(char error_text[]);
double *vector_ne(int nl, int nh);
int *ivector_ne(int nl, int nh);
double *dvector_ne(int nl, int nh);
double **matrix_ne(int nrl, int nrh, int ncl, int nch);
double **dmatrix_ne(int nrl, int nrh, int ncl, int nch);
int **imatrix_ne(int nrl, int nrh, int ncl, int nch);
double **submatrix_ne(double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl);
void free_vector_ne(double *v, int nl, int nh);
void free_ivector_ne(int *v, int nl, int nh);
void free_dvector_ne(double *v, int nl, int nh);
void free_matrix_ne(double **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix_ne(double **m, int nrl, int nrh, int ncl, int nch);
void free_imatrix_ne(int **m, int nrl, int nrh, int ncl, int nch);
void free_submatrix_ne(double **b, int nrl, int nrh, int ncl, int nch);
double **convert_matrix_ne(double *a, int nrl, int nrh, int ncl, int nch);
void free_convert_matrix_ne(double **b, int nrl, int nrh, int ncl, int nch);
double qromb_ne(double (*func)(double), double a, double b);
double trapzd_ne(double (*func)(double), double a, double b, int n);
void polint_ne(double xa[], double ya[], int n, double x, double *y, double *dy);
void spline_ne(double x[], double y[], int n, double yp1, double ypn, double y2[]);
void splint_ne(double xa[], double ya[], double y2a[], int n, double x, double *y);
void gaussj_ne(double **a, int n, double **b, int m);

void nrerror_ne(char error_text[])
{

   /*fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);*/
}

double *vector_ne(int nl, int nh)
{
   double *v;

   v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
   if (!v) nrerror_ne("allocation failure in vector()");
   return v-nl;
}

int *ivector_ne(int nl, int nh)
{
   int *v;

   v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
   if (!v) nrerror_ne("allocation failure in ivector()");
   return v-nl;
}

double *dvector_ne(int nl, int nh)
{
   double *v;

   v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
   if (!v) nrerror_ne("allocation failure in dvector()");
   return v-nl;
}

double **matrix_ne(int nrl, int nrh, int ncl, int nch)
{
   int i;
   double **m;

   m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
   if (!m) nrerror_ne("allocation failure 1 in matrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
      if (!m[i]) nrerror_ne("allocation failure 2 in matrix()");
      m[i] -= ncl;
   }
   return m;
}

double **dmatrix_ne(int nrl, int nrh, int ncl, int nch)
{
   int i;
   double **m;

   m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
   if (!m) nrerror_ne("allocation failure 1 in dmatrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
      if (!m[i]) nrerror_ne("allocation failure 2 in dmatrix()");
      m[i] -= ncl;
   }
   return m;
}

int **imatrix_ne(int nrl, int nrh, int ncl, int nch)
{
   int i,**m;

   m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
   if (!m) nrerror_ne("allocation failure 1 in imatrix()");
   m -= nrl;

   for(i=nrl;i<=nrh;i++) {
      m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
      if (!m[i]) nrerror_ne("allocation failure 2 in imatrix()");
      m[i] -= ncl;
   }
   return m;
}

double **submatrix_ne(double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl)
{
   int i,j;
   double **m;

   m=(double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(double*));
   if (!m) nrerror_ne("allocation failure in submatrix()");
   m -= newrl;

   for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

   return m;
}

void free_vector_ne(double *v, int nl, int nh)
{
   free((char*) (v+nl));
}

void free_ivector_ne(int *v, int nl, int nh)
{
   free((char*) (v+nl));
}

void free_dvector_ne(double *v, int nl, int nh)
{
   free((char*) (v+nl));
}

void free_matrix_ne(double **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void free_dmatrix_ne(double **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void free_imatrix_ne(int **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void free_submatrix_ne(double **b, int nrl, int nrh, int ncl, int nch)
{
   free((char*) (b+nrl));
}

double **convert_matrix_ne(double *a, int nrl, int nrh, int ncl, int nch)
{
   int i,j,nrow,ncol;
   double **m;

   nrow=nrh-nrl+1;
   ncol=nch-ncl+1;
   m = (double **) malloc((unsigned) (nrow)*sizeof(double*));
   if (!m) nrerror_ne("allocation failure in convert_matrix()");
   m -= nrl;
   for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
   return m;
}

void free_convert_matrix_ne(double **b, int nrl, int nrh, int ncl, int nch)
{
   free((char*) (b+nrl));
}

double qromb_ne(double (*func)(double), double a, double b)
{
   double ss, dss;
   double s[JMAXP+1],h[JMAXP+1];
   int j;

   h[1]=1.0;
   for (j=1;j<=JMAX;j++) {
      s[j]=trapzd_ne(func,a,b,j);
      if (j >= K) {
         polint_ne(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
         if (fabs(dss) < EPS*fabs(ss)) return ss;
      }
      /*printf("%2d     %e\n",j,ss);*/
      s[j+1]=s[j];
      h[j+1]=0.25*h[j];
   }
   nrerror_ne("Too many steps in routine QROMB");
   /*     printf("not accurate\n");*/
   return ss;
}

double trapzd_ne(double (*func)(double), double a, double b, int n)
{
   double x,tnm,sum,del;
   static double s;
   static int it;
   int j;

   if (n == 1) {
      it=1;
      return (s=0.5*(b-a)*(FUNC_ne(a)+FUNC_ne(b)));
   } else {
      tnm=it;
      del=(b-a)/tnm;
      x=a+0.5*del;
      for (sum=0.0,j=1;j<=it;j++,x+=del) sum += FUNC_ne(x);
      it *= 2;
      s=0.5*(s+(b-a)*sum/tnm);
      return s;
   }
}

void polint_ne(double xa[], double ya[], int n, double x, double *y, double *dy)
{
   int i,m,ns=1;
   double den,dif,dift,ho,hp,w;
   double *c, *d;

   dif=fabs(x-xa[1]);
   c=vector_ne(1,n);
   d=vector_ne(1,n);
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
         if ( (den=ho-hp) == 0.0) nrerror_ne("Error in routine POLINT");
         den=w/den;
         d[i]=hp*den;
         c[i]=ho*den;
      }
      *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
   }
   free_vector_ne(d,1,n);
   free_vector_ne(c,1,n);
}

void spline_ne(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
   int i,k;
   double p, qn, sig, un, *u;

   u=vector_ne(1,n-1);
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
   free_vector_ne(u,1,n-1);
}

void splint_ne(double xa[], double ya[], double y2a[], int n, double x, double *y)
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
   if (h == 0.0) nrerror_ne("Bad XA input to routine SPLINT");
   a=(xa[khi]-x)/h;
   b=(x-xa[klo])/h;
   *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

void gaussj_ne(double **a, int n, double **b, int m)
{
   int *indxc,*indxr,*ipiv;
   int i, icol=1, irow=1, j, k, l, ll;
   double big,dum,pivinv;

   indxc=ivector_ne(1,n);
   indxr=ivector_ne(1,n);
   ipiv=ivector_ne(1,n);
   for (j=1;j<=n;j++) ipiv[j]=0;
   for (i=1;i<=n;i++) {
      big=0.0;
      for (j=1;j<=n;j++)
         if (ipiv[j] != 1)
            for (k=1;k<=n;k++) {
               if (ipiv[k] == 0) {
                  if (fabs(a[j][k]) >= big) {
                     big=fabs(a[j][k]);
                     irow=j;
                     icol=k;
                  }
               } else if (ipiv[k] > 1) nrerror_ne("GAUSSJ: Singular Matrix-1");
            }
      ++(ipiv[icol]);
      if (irow != icol) {
         for (l=1;l<=n;l++) SWAP_ne(a[irow][l],a[icol][l])
         for (l=1;l<=m;l++) SWAP_ne(b[irow][l],b[icol][l])
      }
      indxr[i]=irow;
      indxc[i]=icol;
      if (a[icol][icol] == 0.0) nrerror_ne("GAUSSJ: Singular Matrix-2");
      pivinv=1.0/a[icol][icol];
      a[icol][icol]=1.0;
      for (l=1;l<=n;l++) a[icol][l] *= pivinv;
      for (l=1;l<=m;l++) b[icol][l] *= pivinv;
      for (ll=1;ll<=n;ll++)
         if (ll != icol) {
            dum=a[ll][icol];
            a[ll][icol]=0.0;
            for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
            for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
         }
   }
   for (l=n;l>=1;l--) {
      if (indxr[l] != indxc[l])
         for (k=1;k<=n;k++)
            SWAP_ne(a[k][indxr[l]],a[k][indxc[l]]);
   }
   free_ivector_ne(ipiv,1,n);
   free_ivector_ne(indxr,1,n);
   free_ivector_ne(indxc,1,n);
}

/* FINISH */


struct VERNER_STRUCT {
  double Eth;
  double Emax;
  double Ezero;
  double s0;
  double ya;
  double P;
  double yw;
  double y0;
  double y1;
};

struct VERNER_PARTIAL_STRUCT {
  int electron;
  int principal;
  int angular;
  double Eth;
  double Ezero;
  double s0;
  double ya;
  double P;
  double yw;
};

static double *E_array;

static double *EGRID,*PIGRID,*PIGRID_2,*RRGRID,*RRGRID_2,THRESHOLD,ANGULAR;
static int GRIDNUM=20000;

static double *LOWE_EGRID,*LOWE_PIGRID,*LOWE_PIGRID_2,*LOWE_RRGRID,*LOWE_RRGRID_2;
static int LOWE_GRIDNUM=6;

static double EBIN,EMIN,EMAX; /* range of spectrum in [eV] */

static double voigt_lim; /* Voigt function parameter */
static double tau_lim; /* Voigt function opacity limit */

static double doppler_rad;

static double *ionizsigmatemp,*tau_exc,*tau_edge,*tau;

static int SPECBINS;

/* XSPEC subroutines */
char* FGMSTR(char* name);

/* Subroutines from PHOTOION */
double dfdE_ne(double g_i, double p0, double p1, double p2, double p3, double E);
double pisigma_ne(double g_i, double p0, double p1, double p2, double p3, double E);
double pispline_ne(double E);
double lowEpispline_ne(double E);
double vernerph_ne(struct VERNER_STRUCT verner, double E); /* verner photoionization sigma */
double vernerpartialph_ne(struct VERNER_PARTIAL_STRUCT verner, double E); /* verner partial photoionization sigma */
double excitsigma_ne(double E0, double OSCILLATOR,double DELTANUD,double ALPHA, double E);
double gauss_ne(double s,double x);   /* norm. gaussian: 1/sqrt(2*PI)/s e^(-x^2/2*s^2) */
double doppler_ne(double v /*[km/s]*/);
double voigt_ne(double ALPHA,double v);

/* NEUTRAL Subroutines */
double fac_ionizsigma_ne(double THRESHOLD, double E);
void line_limits_ne(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double *Elo,double *Ehi);
void line_opacity_ne(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double tau_exc_p[]);
void fac_edge_opacity_ne(double Nion_column_density, double THRESHOLD, double tau_edge_p[]) ;
void verner_full_edge_opacity_ne(double Nion_column_density,double THRESHOLD,struct VERNER_STRUCT verner,double tau_edge_p[]);
void verner_partial_edge_opacity_ne(double Nion_column_density,double THRESHOLD,struct VERNER_PARTIAL_STRUCT verner,double tau_edge_p[]);

int neutral
(float *ear,int ne,float *param,int ifl,float *photar,float *photer)
{
  /*ear[0] -> ear[ne]   in keV
    photar[0] -> photar[ne-1]  
    photar[i] = [photons/cm^2/s] for bin ear[i] -> ear[i+1]
    param[0] -> param[TOTAL-1]   gives XSPEC parameters from 1 to TOTAL*/

  struct VERNER_STRUCT vernerionizsigma[31][31];

  struct VERNER_PARTIAL_STRUCT partialsigma[29][125];

  struct HYDROGEN_STRUCT {
    double lambda[8];
    double f[8];
    double A[8];
    double b[8]; /* branching ratios? */
  };
  struct HYDROGEN_STRUCT hydrogen[27];
  
  struct HELIUM_STRUCT {
    double lambda[11];
    double f[11];
    double A[11];
    double b[11]; /* branching ratios? */
  };
  struct HELIUM_STRUCT helium[27];
  
  struct HIGHER_ORDER_STRUCT {
    double lambda[101];
    double f[101];
  };
  struct HIGHER_ORDER_STRUCT highn[27][3];

  double redshift,v_rad,sigmav_rad,N_e,**Nion;
  int verbose;

  double Elo,Ehi; /* Voigt function param.'s */
  double earBIN,Ewidth;

  double *E_bin,bin;
  
  /* junk values for strings, ints, and floats */
  char *sjunk,*line,*sjunk1,*sjunk2,*sjunk3,*sjunk4;
  double djunk;

  double g_i,g_j;

  double ELO=1.e-3,EHI=1.e6;  /* MAKE SURE THIS RANGE IS OK!!! CHECK HERE!!! */

  double en,Atemp,ftemp,ga;

  char *root,*element_name,*ext,*temp;
  int i,j,k,n;
  /* *.en file */
  double energy;

  double grounddeg[3],leveldeg[3][11],freedeg[3];

  double earlo,earhi;
  
  int HIGHN=50;

  int ELEMENTS=12;

  double E0=0.,DELTANUD,ALPHA,OSCILLATOR;

  int element,electron,element_prev=0,principal,angular;
  int LINE;
  double WAVE,AMtemp,fMtemp;
  double Emax,Ezero,Eth,s0,ya,P,yw,y0,y1;
  double Etemp;

  int *list;

  char name[]="PHOTOION_DIR";
  char *DATADIR;

  double **H_excite, **He_excite;

  double *ABUND,*oshe,oscillatornorm=0.;

  int klo,khi,jlow,jhigh;

  int DIST,lines;

  FILE *vernerphoto,*vernerpartial,*linedat,*highnfile;
  char *vernerphoto_name,*vernerpartial_name,*linedat_name,*highnfile_name;
  FILE *input;

  /* initialize photar array */
  for (i=0;i<ne;++i) photar[i]=0.;

  /* FILE NAMES */
  vernerphoto_name=malloc(200);
  vernerpartial_name=malloc(200);
  linedat_name=malloc(200);
  highnfile_name=malloc(200);

  root=malloc(200);
  element_name=malloc(3);  /* 2-char symbols ("Ne") need 3 bytes with the NUL */
  ext=malloc(30);
  temp=malloc(130);
  sjunk=malloc(50);
  sjunk1=malloc(50);
  sjunk2=malloc(50);
  sjunk3=malloc(50);
  sjunk4=malloc(50);
  line=malloc(400);

  Nion=dmatrix_ne(1,28,1,28);
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) Nion[i][j]=0.;

  DATADIR=FGMSTR(name);

  ABUND=dvector_ne(1,30);
  for (i=1;i<=30;++i) ABUND[i]=0.;
  sprintf(temp,"%s/photoion_dat/abundance.dat",DATADIR);
  input=fopen(temp,"r");
  if ( input == NULL ) {
    printf("NEUTRAL: Failed to open %s\n", temp);
    return 1;
  }
  while (fscanf(input,"%d%lf",&element,&djunk)!=EOF) {
    ABUND[element]=djunk;
  }
  fclose(input);

  oshe=dvector_ne(1,30);
  sprintf(temp,"%s/photoion_dat/oscillator_he.dat",DATADIR);
  input=fopen(temp,"r");
  while (fscanf(input,"%d%lf",&element,&djunk)!=EOF) {
    oshe[element]=cube(10.)*djunk;
  }
  fclose(input);

  Nion[1][1]=param[0];
  /*djunk=param[1]; ABUND[2]*=djunk;*/  Nion[2][2]=ABUND[2]*Nion[1][1];
  /*djunk=param[2]; ABUND[6]*=djunk;*/  Nion[6][6]=ABUND[6]*Nion[1][1];
  /*djunk=param[3]; ABUND[7]*=djunk;*/  Nion[7][7]=ABUND[7]*Nion[1][1];
  /*djunk=param[4]; ABUND[8]*=djunk;*/  Nion[8][8]=ABUND[8]*Nion[1][1];
  /*djunk=param[5]; ABUND[10]*=djunk;*/  Nion[10][10]=ABUND[10]*Nion[1][1];
  /*djunk=param[6]; ABUND[12]*=djunk;*/  Nion[12][12]=ABUND[12]*Nion[1][1];
  /*djunk=param[7]; ABUND[13]*=djunk;*/  Nion[13][13]=ABUND[13]*Nion[1][1];
  /*djunk=param[8]; ABUND[14]*=djunk;*/  Nion[14][14]=ABUND[14]*Nion[1][1];
  /*djunk=param[9]; ABUND[16]*=djunk;*/  Nion[16][16]=ABUND[16]*Nion[1][1];
  /*djunk=param[10]; ABUND[18]*=djunk;*/  Nion[18][18]=ABUND[18]*Nion[1][1];
  /*djunk=param[11]; ABUND[20]*=djunk;*/  Nion[20][20]=ABUND[20]*Nion[1][1];
  /*djunk=param[12]; ABUND[26]*=djunk;*/  Nion[26][26]=ABUND[26]*Nion[1][1];
  /*djunk=param[13]; ABUND[28]*=djunk;*/  Nion[28][28]=ABUND[28]*Nion[1][1];
  sigmav_rad = 1.e5*param[1];   /* gaussian (1-sigma) velocity width: convert [km/s] to [cm/s] */
  if (sigmav_rad) lines=1;
  else lines=0;
  redshift=0.;     /* cosmological redshift */
  v_rad = 0.;  /* velocity shift: convert [km/s] to [cm/s] */
  EMIN=1000.*1.e-3; /* convert keV to eV */
  EMAX=1000.*15.; /* convert keV to eV */
  SPECBINS=100000;

  verbose=0;

  for (element=1;element<=28;++element) {
    for (electron=1;electron<=28;++electron) {
      if (Nion[element][electron]) {
        Nion[element][electron]=fabs(Nion[element][electron]);
	if (verbose) printf("Nion[%2d][%2d]=%e\n",element,electron,Nion[element][electron]);
      }
    }
  }

  /* To put EMIN and EMAX in source rest-frame energy units [eV]
     We convert back at the very end */  
  EMIN=EMIN*(1.+redshift);
  EMAX=EMAX*(1.+redshift);
 
  /* Decrease these to improve accuracy */
  voigt_lim=1.e-4;
  tau_lim=1.e-5;

  doppler_rad=doppler_ne(v_rad);

  EGRID=dvector_ne(1,GRIDNUM);  
  PIGRID=dvector_ne(1,GRIDNUM); 
  PIGRID_2=dvector_ne(1,GRIDNUM);
  RRGRID=dvector_ne(1,GRIDNUM);  
  RRGRID_2=dvector_ne(1,GRIDNUM);

  LOWE_EGRID=dvector_ne(1,LOWE_GRIDNUM);  
  LOWE_PIGRID=dvector_ne(1,LOWE_GRIDNUM); 
  LOWE_PIGRID_2=dvector_ne(1,LOWE_GRIDNUM);
  LOWE_RRGRID=dvector_ne(1,LOWE_GRIDNUM); 
  LOWE_RRGRID_2=dvector_ne(1,LOWE_GRIDNUM);

  E_array=dvector_ne(1,SPECBINS);       /* energy axis */
  tau=dvector_ne(1,SPECBINS);           /* total opacity in all ions */
  tau_exc=dvector_ne(1,SPECBINS);        
  tau_edge=dvector_ne(1,SPECBINS);       
  ionizsigmatemp=dvector_ne(1,SPECBINS);
  H_excite=dmatrix_ne(1,28,1,6);
  He_excite=dmatrix_ne(1,28,1,9);
  list=ivector_ne(1,ELEMENTS);

  for (i=1;i<=SPECBINS;++i) {
    tau[i]=0.;
    tau_exc[i]=0.;
    tau_edge[i]=0.;
  }

  /* Atomic number for C,N,O,Ne,Mg,Al,Si,S,Ar,Ca,Fe */
  list[1]=6;
  list[2]=7;
  list[3]=8;
  list[4]=10;
  list[5]=12;
  list[6]=13;
  list[7]=14;
  list[8]=16;
  list[9]=18;
  list[10]=20;
  list[11]=26;
  list[12]=28;

  /* uniformly-spaced energy bin size */
  EBIN=(EMAX-EMIN)/((double) SPECBINS);
  for (k=1;k<=SPECBINS;++k) {
    E_array[k]=((double) k-0.5)*EBIN+EMIN;
  }

  /* for velocity convolution out to >= 4 sigma */
  DIST=(int) (4.*sigmav_rad/ccc*EMAX/EBIN); /*  +/-4 sigma */
  if (verbose) printf("DIST = %d\n",DIST);

  /* logarithmic energy grid for Fe-L RR and PI cross sections*/
  for (k=1;k<=GRIDNUM;++k) {
    /*EGRID=log10(energy)*/
    EGRID[k]=((double) k-1.)/((double) GRIDNUM)*log10(EHI/ELO)+log10(ELO);
  }

  /*  Reading in Total Photoionization Cross Sections (Verner)  */
  sprintf(vernerphoto_name,"%s/photoion_dat/verner_photo.dat",DATADIR);
  vernerphoto=fopen(vernerphoto_name,"r");
  while (fscanf(vernerphoto,"%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf",&element,&electron,&Eth,&Emax,&Ezero,&s0,&ya,&P,&yw,&y0,&y1)!=EOF) {
    vernerionizsigma[element][electron].Eth=Eth;
    vernerionizsigma[element][electron].Emax=Emax;
    vernerionizsigma[element][electron].Ezero=Ezero;
    vernerionizsigma[element][electron].s0=s0;
    vernerionizsigma[element][electron].ya=ya;
    vernerionizsigma[element][electron].P=P;
    vernerionizsigma[element][electron].yw=yw;
    vernerionizsigma[element][electron].y0=y0;
    vernerionizsigma[element][electron].y1=y1;
  }
  fclose(vernerphoto);
  free(vernerphoto_name);

  /*  Reading in Partial Photoionization Cross Sections (Verner)  */
  sprintf(vernerpartial_name,"%s/photoion_dat/verner_partial_PIsigmas.dat",DATADIR);
  vernerpartial=fopen(vernerpartial_name,"r");
  k=0;
  while (fscanf(vernerpartial,"%d%d%d%d%lf%lf%lf%lf%lf%lf",&element,&electron,&principal,&angular,&Eth,&Ezero,&s0,&ya,&P,&yw)!=EOF) {
    if (element!=element_prev) k=0;
    partialsigma[element][k].electron=electron;
    partialsigma[element][k].principal=principal;
    partialsigma[element][k].angular=angular;
    partialsigma[element][k].Eth=Eth;
    partialsigma[element][k].Ezero=Ezero;
    partialsigma[element][k].s0=s0;
    partialsigma[element][k].ya=ya;
    partialsigma[element][k].P=P;
    partialsigma[element][k].yw=yw;
    element_prev=element;
    k=k+1;
  }
  fclose(vernerpartial);
  free(vernerpartial_name);

  /* Total electron Thomson depth */
  for (k=1;k<=SPECBINS;++k) {
    tau[k]+=N_e*sigmaT;
  }

  if (lines) {
    if (verbose) printf("Determining LOW-n Photoexcitation Cross Sections & Opacity for C to Fe...\n");
    grounddeg[1]=2.;
    freedeg[1]=1.;
    for (i=1;i<=6;++i) leveldeg[1][i]=6.;
    grounddeg[2]=1.;
    freedeg[2]=2.;
    for (i=3;i<=8;++i) leveldeg[2][i]=3.;
    sprintf(linedat_name,"%s/photoion_dat/line.dat",DATADIR);
    linedat=fopen(linedat_name,"r");
    while(fscanf(linedat,"%s%d%d%d%s%lf%s%lf%lf",sjunk,&element,&electron,&LINE,sjunk,&WAVE,sjunk,&AMtemp,&fMtemp) != EOF) {
      if (electron == 1) {
	hydrogen[element].lambda[LINE]=WAVE;
	hydrogen[element].A[LINE]=AMtemp;
	hydrogen[element].f[LINE]=fMtemp;
	if (Nion[element][electron] && fMtemp && LINE <= 4) {
	  E0=1000.*lambdatokeV/hydrogen[element].lambda[LINE];
	  E0=E0*doppler_rad;
	  OSCILLATOR=hydrogen[element].f[LINE];
	  g_j=leveldeg[electron][LINE];
	  g_i=grounddeg[electron];
	  Atemp=hydrogen[element].A[LINE];
	  ga=Atemp/*g_j*Atemp/(3.*g_i*OSCILLATOR)*/;
	  DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	  ALPHA=ga/(4.*PI*DELTANUD);
	  line_opacity_ne(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc);
	}
      } else if (electron == 2) {
	helium[element].lambda[LINE]=WAVE;
	helium[element].A[LINE]=AMtemp;
	helium[element].f[LINE]=fMtemp;
	if (Nion[element][electron] && fMtemp  && LINE <= 6) {
	  E0=1000.*lambdatokeV/helium[element].lambda[LINE];
	  E0=E0*doppler_rad;
	  OSCILLATOR=helium[element].f[LINE];
	  g_j=leveldeg[electron][LINE];
	  g_i=grounddeg[electron];
	  Atemp=helium[element].A[LINE];
	  ga=Atemp/*g_j*Atemp/(3.*g_i*OSCILLATOR)*/;
	  DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	  ALPHA=ga/(4.*PI*DELTANUD);
	  line_opacity_ne(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc);
	}
      }
    }
    fclose(linedat);
    free(linedat_name);
    
    if (verbose) printf("Determining HIGH-n Photoexcitation Cross Sections & Opacity for C to Fe\n");
    sprintf(highnfile_name,"%s/photoion_dat/highn.dat",DATADIR);
    highnfile=fopen(highnfile_name,"r");
    while (fscanf(highnfile,"%d%d%d%lf%lf",&element,&electron,&n,&Etemp,&ftemp)!=EOF) {
      highn[element][electron].lambda[n]=lambdatokeV/Etemp*1000.;
      highn[element][electron].f[n]=ftemp;
    }
    fclose(highnfile);
    free(highnfile_name);
    for (i=1;i<=ELEMENTS;++i) {
      element=list[i];
      for (electron=1;electron<=2;++electron) {
	if (Nion[element][electron]) {
	  if (electron==1) oscillatornorm=1.6; /* Bethe-Salpeter p. 265 */
	  if (electron==2) oscillatornorm=oshe[element]; /* defined above */
	  for (n=6;n<=HIGHN;++n) {
	    OSCILLATOR=oscillatornorm/cube((double) n);
	    if (element!=28) {
	      E0=lambdatokeV/highn[element][electron].lambda[n]*1000.;
	    } else if (electron==1) {/* Use Fe numbers for Ni */
	      E0=lambdatokeV/(highn[26][electron].lambda[n]/1.1614)*1000.;
	    } else if (electron==2) {/* Use Fe numbers for Ni */
	      E0=lambdatokeV/(highn[26][electron].lambda[n]/1.165)*1000.;
	    }
	    E0=E0*doppler_rad;
	    DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	    /* taking classical value for "ga=gamma" from p. 112-114 B+D */
	    ga=2.47*pow(10.,-22.)*sqr(E0*eVtoergs/hhh);
	    ALPHA=ga/(4.*PI*DELTANUD);
	    line_opacity_ne(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc);
	  }
	}
      }
    }
    
    
    /* File root for L-shell input files for Ne through Ca */
    if (verbose) printf("L-shell ions:  Ne through Ni\n");
    for (element=10;element<=28;++element) {
      for (electron=1;electron<=10;++electron) {
	/* Need H- and He-like photoionization cross-sections */
	if (Nion[element][electron] && (electron>=3 || element==28)) {
	  if (element == 10) sprintf(element_name,"Ne");
	  if (element == 12) sprintf(element_name,"Mg");
	  if (element == 13) sprintf(element_name,"Al");
	  if (element == 14) sprintf(element_name,"Si");
	  if (element == 16) sprintf(element_name,"S");
	  if (element == 18) sprintf(element_name,"Ar");
	  if (element == 20) sprintf(element_name,"Ca");
	  if (element == 26) sprintf(element_name,"Fe");
	  if (element == 28) sprintf(element_name,"Ni");
	  sprintf(root,"%s/photoion_dat/L_shell/%s",DATADIR,element_name);
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  if (electron>=3) {
	    ext="tr_short"; 
	    if (electron<10) sprintf(temp,"%s0%da.%s",root,electron,ext);
	    else sprintf(temp,"%s%da.%s",root,electron,ext);
	    input=fopen(temp,"r");
	    while (fscanf(input,"%d%lf%d%lf%lf%lf%lf",&j,&g_j,&i,&g_i,&en,&ftemp,&Atemp) != EOF) {
	      j=j+1;
	      i=i+1;
	      g_i=g_i+1.;
	      g_j=g_j+1.;
	      ftemp=ftemp/g_i;     /* CHECK THIS - VERY IMPORTANT!!! */
	      E0=en;
	      E0=E0*doppler_rad;
	      DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	      if ((E0>=EMIN && E0<=EMAX) && ftemp>tau_lim*FACTOR*DELTANUD/Nion[element][electron]) {
		OSCILLATOR=ftemp;
		ga=Atemp/*g_j*Atemp/(3.*g_i*OSCILLATOR)*/;
		ALPHA=ga/(4.*PI*DELTANUD);
		line_opacity_ne(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc);
	      }
	    }
	    fclose(input);
	  }
	}
      }
    }
    
    /* File root for L-shell input files for Ne through Ca */
    if (verbose) printf("L-shell ions:  C through O\n");
    for (element=6;element<=8;++element) {
      for (electron=3;electron<=8;++electron) {
	if (Nion[element][electron]) {
	  if (element == 6) sprintf(element_name,"C");
	  if (element == 7) sprintf(element_name,"N");
	  if (element == 8) sprintf(element_name,"O");
	  sprintf(root,"%s/photoion_dat/L_shell/%s",DATADIR,element_name);
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  ext="tr_short"; 
	  if (electron<10) sprintf(temp,"%s0%da.%s",root,electron,ext);
	  else sprintf(temp,"%s%da.%s",root,electron,ext);
	  input=fopen(temp,"r");
	  while (fscanf(input,"%d%lf%d%lf%lf%lf%lf",&j,&g_j,&i,&g_i,&en,&ftemp,&Atemp) != EOF) {
	    j=j+1;
	    i=i+1;
	    g_i=g_i+1.;
	    g_j=g_j+1.;
	    ftemp=ftemp/g_i;     /* CHECK THIS - VERY IMPORTANT!!! */
	    E0=en;
	    E0=E0*doppler_rad;
	    DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	    if ((E0>=EMIN && E0<=EMAX) && ftemp>tau_lim*FACTOR*DELTANUD/Nion[element][electron]) {
	      OSCILLATOR=ftemp;
	      ga=Atemp/*g_j*Atemp/(3.*g_i*OSCILLATOR)*/;
	      ALPHA=ga/(4.*PI*DELTANUD);
	      line_opacity_ne(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc);
	    }
	  }
	  fclose(input);
    	}
      }
    }
    
    
    /* M-shell ions */
    if (verbose) printf("M-shell ions:  Mg through Ni\n");
    for (element=12;element<=28;++element) {
      for (electron=11;electron<=28;++electron) {
	if (Nion[element][electron]) {
	  if (element == 12) sprintf(element_name,"Mg");
	  if (element == 13) sprintf(element_name,"Al");
	  if (element == 14) sprintf(element_name,"Si");
	  if (element == 16) sprintf(element_name,"S");
	  if (element == 18) sprintf(element_name,"Ar");
	  if (element == 20) sprintf(element_name,"Ca");
	  if (element == 26) sprintf(element_name,"Fe");
	  if (element == 28) sprintf(element_name,"Ni");
	  sprintf(root,"%s/photoion_dat/M_shell/%s",DATADIR,element_name);
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  
	  /* Transitions */
	  ext="tr_short"; 
	  sprintf(temp,"%s%2da.%s",root,electron,ext);
	  input=fopen(temp,"r");
	  while (fscanf(input,"%d%lf%d%lf%lf%lf%lf",&j,&g_j,&i,&g_i,&en,&ftemp,&Atemp) != EOF) {
	    g_j=g_j+1.;
	    g_i=g_i+1.;
	    ftemp=ftemp/g_i;     /* CHECK THIS - VERY IMPORTANT!!! */
	    E0=en;
	    E0=E0*doppler_rad;
	    DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	    if ((E0>=EMIN && E0<=EMAX) && ftemp>tau_lim*FACTOR*DELTANUD/Nion[element][electron]) {
	      OSCILLATOR=ftemp;
	      ga=Atemp/*g_j*Atemp/(3.*g_i*OSCILLATOR)*/;
	      ALPHA=ga/(4.*PI*DELTANUD);
	      line_opacity_ne(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc);
	    }
	  }
	  fclose(input);
	}      
      }
    }
  }
  
  /* Read in edge opacity */
  sprintf(temp,"%s/photoion_dat/neutral.tau",DATADIR);
  input=fopen(temp,"r");
  j=1;

  E_bin=dvector_ne(1,SPECBINS);
  /*  for (k=1;k<=3;++k) fgets(line,LENGTH,input);    */
  while (fgets(line,SPECBINS,input) != NULL) {
    sscanf(line,"%lf%lf%lf",&energy /* keV */,&bin /* keV */,&djunk /* opacity */);
    
    bin=2.*bin;
    /*    E_array[j]=1000.*energy;*/
    E_bin[j]=1000.*bin;
    tau_edge[j]=Nion[1][1]*djunk;
    ++j;
  }
  fclose(input);

  if (0 && sigmav_rad) {
    if (verbose) printf("Convolving spectrum with appropriate velocity distribution...");
    klo=1;
    khi=SPECBINS;
    for (k=klo;k<=khi;++k) {
      if (k<=DIST) {jlow = 1;      jhigh = k+DIST;} 
      else         {jlow = k-DIST; jhigh = k+DIST;}
      if (jlow<1) jlow=1;
      if (jhigh>SPECBINS) jhigh=SPECBINS;
      for (n=jlow;n<=jhigh;++n) {
	djunk=EBIN/(E_array[k])*gauss_ne(sigmav_rad/ccc,(E_array[n]-E_array[k])/E_array[k])*tau_edge[n];
	tau[k]+=djunk;
      }
    }
    if (verbose) printf("done.\n");
  } else {
    for (k=1;k<=SPECBINS;++k) tau[k]+=tau_edge[k];
  }


  for (k=1;k<=SPECBINS;++k) tau[k]+=tau_exc[k];

  if (verbose) printf("Calculating final spectrum...");
  /*    E_redshift=(1.+redshift)*1000.*(ear[i]+ear[i+1])/2.; */
  j=1;
  i=0;
  earlo=1000.*(1.+redshift)*ear[i];
  earhi=1000.*(1.+redshift)*ear[i+1];
  earBIN=earhi-earlo;
  Elo=E_array[j]-EBIN/2.;
  Ehi=E_array[j]+EBIN/2.;
  while (earhi<=Elo && i<ne-1) {
    ++i;
    earlo=1000.*(1.+redshift)*ear[i];
    earhi=1000.*(1.+redshift)*ear[i+1];
    earBIN=earhi-earlo;
  }
  while (Ehi<earlo && j<SPECBINS) {
    ++j;
    Elo=E_array[j]-EBIN/2.;
    Ehi=E_array[j]+EBIN/2.;
  }
  while (j<SPECBINS && i<ne-1) {
    if (Elo<earlo) {
      if (Ehi-earlo<earBIN) Ewidth=Ehi-earlo;
      else Ewidth=earBIN;
    } else {
      if (earhi-Elo<EBIN) Ewidth=earhi-Elo;
      else Ewidth=EBIN;
    }
    photar[i]+=Ewidth/earBIN*exp(-tau[j]);
    if (earhi<Ehi) {
      ++i;
      earlo=1000.*(1.+redshift)*ear[i];
      earhi=1000.*(1.+redshift)*ear[i+1];
      earBIN=earhi-earlo;
    } else {
      ++j;
      Elo=E_array[j]-EBIN/2.;
      Ehi=E_array[j]+EBIN/2.;
    }
  }
  if (verbose) printf("done.\n");


  if (verbose) printf("Freeing memory...");
  /* Free all the memory */
  free(root);
  free(temp);
  free(sjunk);
  free(sjunk1);
  free(sjunk2);
  free(sjunk3);
  free(sjunk4);
  free(line);
  /*  free(element_name);
      free(ext);*/
  free_dvector_ne(LOWE_EGRID,1,LOWE_GRIDNUM);  
  free_dvector_ne(LOWE_PIGRID,1,LOWE_GRIDNUM); 
  free_dvector_ne(LOWE_PIGRID_2,1,LOWE_GRIDNUM);
  free_dvector_ne(LOWE_RRGRID,1,LOWE_GRIDNUM); 
  free_dvector_ne(LOWE_RRGRID_2,1,LOWE_GRIDNUM);
  free_dvector_ne(ABUND,1,30);
  free_dvector_ne(oshe,1,30);
  free_dvector_ne(E_array,1,SPECBINS);       
  free_dvector_ne(tau,1,SPECBINS);           
  free_dvector_ne(tau_exc,1,SPECBINS);           
  free_dvector_ne(tau_edge,1,SPECBINS);           
  free_dvector_ne(ionizsigmatemp,1,SPECBINS);
  free_dvector_ne(E_bin,1,SPECBINS);       
  free_dmatrix_ne(Nion,1,28,1,28);
  free_dmatrix_ne(H_excite,1,28,1,6);
  free_dmatrix_ne(He_excite,1,28,1,9);
  free_ivector_ne(list,1,ELEMENTS);
  free_dvector_ne(EGRID,1,GRIDNUM);  
  free_dvector_ne(PIGRID,1,GRIDNUM); 
  free_dvector_ne(PIGRID_2,1,GRIDNUM);
  free_dvector_ne(RRGRID,1,GRIDNUM);  
  free_dvector_ne(RRGRID_2,1,GRIDNUM);
  if (verbose) printf("...done!\n");
  
  return 0;
}

double doppler_ne(double v /*[cm/s]*/) /*[unitless] red: v>0.,d(v)>1, blue: v<0.,d(v)<1*/
{
  double answer;

  answer=sqrt((1.-v/ccc)/(1.+v/ccc));
  return answer;
}

/* For convolution: Normalized gaussian profile centered at 0. with s=sigma */
double gauss_ne(double s,double x) 
{
  double answer=0.;

  answer=1./sqrt(2.*PI)/s*exp(-sqr((x)/s)/2.);  
  return answer;
}

/* Photoexcitation cross section [cm^2] from ground.
   global: E0, DELTANUD, OSCILLATOR, ALPHA */
double excitsigma_ne(double E0, double OSCILLATOR,double DELTANUD,double ALPHA, double E)
{
  double answer;
  double V;

  V=(E-E0)/(hhh*DELTANUD*ergstoeV);  
  answer=PI*re*ccc*OSCILLATOR/sqrt(PI)/DELTANUD*voigt_ne(ALPHA,V)*doppler_rad;
  return answer;
}

/* Differential oscillator strength in atomic units: 1/Hartree = 1/(2.*13.6 eV) */
double dfdE_ne(double g_i,double p0, double p1, double p2, double p3, double Te)
{
  double x,y,answer,E;

  E=Te+THRESHOLD;
  x=(Te+p3)/p3;
  y=(1.+p2)/(sqrt(x)+p2);
  answer=E/(Te+p3)*p0*pow(x,-3.5-ANGULAR+0.5*p1)*pow(y,p1)/g_i;
  return answer;
}

/* Photoionization cross-section in atomic units */
double pisigma_ne(double g_i, double p0, double p1, double p2, double p3, double E)
{
  double answer;

  answer=sqr(a0)*2.*PI*FINE_STRUCTURE*dfdE_ne(g_i,p0,p1,p2,p3,E-THRESHOLD);
  return answer;
}

double pispline_ne(double E)
{
  double answer;
  
  E=log10(E);
  splint_ne(EGRID,PIGRID,PIGRID_2,GRIDNUM,E,&answer);
  return pow(10.,answer);
}

double lowEpispline_ne(double E)
{
  double answer;
  
  E=log10(E);
  splint_ne(LOWE_EGRID,LOWE_PIGRID,LOWE_PIGRID_2,LOWE_GRIDNUM,E,&answer);
  return pow(10.,answer);
}

void line_limits_ne(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double *Elo,double *Ehi)
{
  double Vlim;
  
  /*tau_lim=Nion_column_density*excitsigma=Nion_column_density*PI*re*ccc*OSCILLATOR/sqrt(PI)/DELTANUD*voigt_ne(ALPHA,V);*/
  
  Vlim=sqrt(ALPHA*Nion_column_density*re*ccc*OSCILLATOR/DELTANUD/tau_lim);
  if (voigt_ne(ALPHA,Vlim)>voigt_lim) {
    Vlim=sqrt(ALPHA/sqrt(PI)/(voigt_lim*voigt_ne(ALPHA,0.)));
  }
  if (E0-Vlim*(hhh*DELTANUD*ergstoeV)>EMIN) {
    *Elo=E0-Vlim*(hhh*DELTANUD*ergstoeV);
  } else *Elo=(1.+SMALL)*EMIN;
  if (E0+Vlim*(hhh*DELTANUD*ergstoeV)<EMAX) {
    *Ehi=E0+Vlim*(hhh*DELTANUD*ergstoeV);
  } else *Ehi=(1.-SMALL)*EMAX;

}

void line_opacity_ne(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double tau_exc_p[])
{
  int SUMlo,SUMhi;
  double Ehi,Elo,djunk;
  int k,klo,khi;

  line_limits_ne(Nion_column_density,E0,OSCILLATOR,ALPHA,DELTANUD,&Elo,&Ehi);
  SUMlo=(int) ((Elo-EMIN+0.5*EBIN)/EBIN);
  if ((Elo-EMIN+0.5*EBIN)/EBIN-(double) SUMlo >= 0.5) ++SUMlo;
  SUMhi=(int) ((Ehi-EMIN+0.5*EBIN)/EBIN);
  if ((Ehi-EMIN+0.5*EBIN)/EBIN-(double) SUMlo >= 0.5) ++SUMhi;
  if (SUMlo < 1) SUMlo=1;
  if (SUMhi > SPECBINS) SUMhi=SPECBINS;
  klo=SUMlo;
  khi=SUMhi;

  for (k=klo;k<=khi;++k) {
    djunk=Nion_column_density*excitsigma_ne(E0,OSCILLATOR,DELTANUD,ALPHA,E_array[k]);
    tau_exc_p[k]+=djunk;
  }
}

/* ground-state photoionization cross-section from Verner */
double vernerph_ne(struct VERNER_STRUCT verner, double E)
{
  double answer=0.;
  double Eth,Emax,Ezero,s0,ya,P,yw,y0,y1;
  double F,x,y;

  Eth=verner.Eth;
  Emax=verner.Emax;
  if (E>=Eth && E<=Emax) {
    Ezero=verner.Ezero;
    s0=verner.s0;
    ya=verner.ya;
    P=verner.P;
    yw=verner.yw;
    y0=verner.y0;
    y1=verner.y1;

    x=E/Ezero-y0;
    y=sqrt(sqr(x)+sqr(y1));
    F=(sqr(x-1.)+sqr(yw))*pow(y,0.5*P-5.5)*pow(1.+sqrt(y/ya),-P);
    answer=1.e-18*s0*F;
  }
  return answer;
}

/* ground-state photoionization cross-section from Verner */
double vernerpartialph_ne(struct VERNER_PARTIAL_STRUCT verner, double E)
{
  double answer=0.;
  double Eth,Ezero,s0,ya,P,yw,l;
  int angular;
  double F,y;

  Eth=verner.Eth;
  if (E>=Eth) {
    Ezero=verner.Ezero;
    s0=verner.s0;
    ya=verner.ya;
    P=verner.P;
    yw=verner.yw;
    angular=verner.angular;
    l=(double) angular;

    y=E/Ezero;
    F=(sqr(y-1.)+sqr(yw))*pow(y,-(5.5+l-0.5*P))*pow(1.+sqrt(y/ya),-P);
    answer=1.e-18*s0*F;
  }
  return answer;
}

double fac_ionizsigma_ne(double THRESHOLD, double E)
{
  double answer=0.;
  
  if (E>=THRESHOLD) {
    if (E<=pow(10.,LOWE_EGRID[1])) answer=lowEpispline_ne(E);
    else answer=pispline_ne(E);
  }
  return answer;
}

void fac_edge_opacity_ne(double Nion_column_density, double THRESHOLD, double tau_p[]) 
{
  int k,klo,kth;
  double taujunk;
  
  klo=(int) ((0.95*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.05*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  taujunk=Nion_column_density*fac_ionizsigma_ne(THRESHOLD,E_array[k]/doppler_rad);
  while ((taujunk > tau_lim || k <= kth) && k < SPECBINS) {
    k=k+1;
    taujunk=Nion_column_density*fac_ionizsigma_ne(THRESHOLD,E_array[k]/doppler_rad);
    tau_p[k]+=taujunk;
  }
}

void verner_full_edge_opacity_ne(double Nion_column_density, double THRESHOLD, struct VERNER_STRUCT verner, double tau_p[]) 
{
  int k,klo,kth;
  double taujunk;
  
  klo=(int) ((0.95*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.05*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  taujunk=Nion_column_density*vernerph_ne(verner, E_array[k]/doppler_rad);
  while ((taujunk > tau_lim || k <= kth) && k < SPECBINS) {
    k=k+1;
    taujunk=Nion_column_density*vernerph_ne(verner, E_array[k]/doppler_rad);
    tau_p[k]+=taujunk;
  }
}

void verner_partial_edge_opacity_ne(double Nion_column_density, double THRESHOLD, struct VERNER_PARTIAL_STRUCT verner, double tau_p[]) 
{
  int k,klo,kth;
  double taujunk;
  
  klo=(int) ((0.95*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.05*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  taujunk=Nion_column_density*vernerpartialph_ne(verner, E_array[k]/doppler_rad);
  while ((taujunk > tau_lim || k <= kth) && k < SPECBINS) {
    k=k+1;
    taujunk=Nion_column_density*vernerpartialph_ne(verner, E_array[k]/doppler_rad);
    tau_p[k]+=taujunk;
  }
}

double voigt_ne(double alpha,double v)
{
  int i;
  double *a,*b,*c;
  double v2,v3,fac1,fac2;
  double p1,p2,p3,p4,p5,p6,p7;
  double o1,o2,o3,o4,o5,o6,o7;
  double q1,q2;
  double r1,r2;
  double H;
 
  a=dvector_ne(1,7);
  b=dvector_ne(1,7);
  c=dvector_ne(1,7);

  a[1]=122.607931777104326;
  a[2]=214.382388694706425;
  a[3]=181.928533092181549;
  a[4]=93.155580458134410;
  a[5]=30.180142196210589;   
  a[6]=5.912626209773153;
  a[7]=0.564189583562615;

  b[1]=122.607931773875350;
  b[2]=352.730625110963558;
  b[3]=457.334478783897737;
  b[4]=348.703917719495792;
  b[5]=170.354001821091472;
  b[6]=53.992906912940207;
  b[7]=10.479857114260399;
  
  c[1]=0.5641641;
  c[2]=0.8718681;
  c[3]=1.474395;
  c[4]=-19.57862;
  c[5]=802.4513;
  c[6]=-4850.316;
  c[7]=8031.468;
  
  if (alpha <= .001 && v >= 2.5) {
    v2   = v * v;
    v3   = 1.0;
    fac1 = c[1];
    fac2 = c[1] * (v2 - 1.0);
      
    for (i=1;i<=7;++i) {
      v3     = v3 * v2;
      fac1 = fac1 + c[i] / v3;
      fac2 = fac2 + c[i] / v3 * (v2 - (double) i);
    }
    
    H = exp(-v2) * (1. + sqr(alpha) * (1. - 2.*v2)) + fac1 * (alpha/v2);
    
  } else { 
    p1 = alpha;
    o1 = -v;
    p2 = (p1 * alpha + o1 * v);
    o2 = (o1 * alpha - p1 * v);
    p3 = (p2 * alpha + o2 * v);
    o3 = (o2 * alpha - p2 * v);
    p4 = (p3 * alpha + o3 * v);
    o4 = (o3 * alpha - p3 * v);
    p5 = (p4 * alpha + o4 * v);
    o5 = (o4 * alpha - p4 * v);
    p6 = (p5 * alpha + o5 * v);
    o6 = (o5 * alpha - p5 * v);
    p7 = (p6 * alpha + o6 * v);
    o7 = (o6 * alpha - p6 * v);
    
    q1 = a[1] + p1 * a[2] + p2 * a[3] + p3 * a[4] +
      p4 * a[5] + p5 * a[6] + p6 * a[7];
    r1 =        o1 * a[2] + o2 * a[3] + o3 * a[4] +
      o4 * a[5] + o5 * a[6] + o6 * a[7];
    q2 = b[1] + p1 * b[2] + p2 * b[3] + p3 * b[4] +
      p4 * b[5] + p5 * b[6] + p6 * b[7] + p7;
    r2 =        o1 * b[2] + o2 * b[3] + o3 * b[4] +
      o4 * b[5] + o5 * b[6] + o6 * b[7] + o7;
    
    H = (q1 * q2 + r1 * r2) / (q2 * q2 + r2 * r2);
  }

  free_dvector_ne(a,1,7);
  free_dvector_ne(b,1,7);
  free_dvector_ne(c,1,7);
  
  return H;
}
