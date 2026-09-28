#include "cfortran.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void nrerror_ex(char error_text[]);

int mulext
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(mulext,MULEXT,mulext,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define ccc (2.99792458e10)         /* speed of light [cm/s] */
#define lambdatokeV (12.39841856)    /* lambda=12.39841856/E_keV */

/* 
Reads in file either in energy (E_or_l = 0) or wavelength (E_or_l = 1) units.  
This file is either generated within XSPEC ("wdata" command) or with the 
"photoion" model, which automatically generates "E_spectrum_#.qdp" or 
"l_spectrum_#.qdp".

The input file must be in the same directory in which XSPEC is running and must
have the name "mulext.qdp":

   If E_or_l == 0
     Format of file:
                     (SKIP first three lines)
     Energy [keV]     Half-Bin [keV]     Spectrum [photons/cm^2/s/keV]

     If E_or_l == 1
     Format of file:
                     (SKIP first three lines)
     lambda [A]       Half-Bin [A]       Spectrum [photons/cm^2/s/A]
*/  

/* Forward declarations (converted from K&R style). */
void nrerror_me(char error_text[]);
double *dvector_me(int nl, int nh);
void free_dvector_me(double *v, int nl, int nh);

void nrerror_me(char error_text[])
{

   /*fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);*/
}

double *dvector_me(int nl, int nh)
{
   double *v;

   v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
   if (!v) nrerror_ex("allocation failure in dvector()");
   return v-nl;
}

void free_dvector_me(double *v, int nl, int nh)
{
   free((char*) (v+nl));
}

int mulext
(float *ear,int ne,float *param,int ifl,float *photar,float *photer)
{
  FILE *specfile;
  double energy_or_lambda,energy=0.,bin,spectrum,*E_array,*E_bin,*E_spectrum;
  double Elo,Ehi,earlo,earhi,earBIN,Ewidth,redshift,v,doppler_v;
  int i,j,k,LENGTH,E_or_l,SPECBINS;
  char *line;

  LENGTH=1000;
  line=malloc(LENGTH);

  for (i=0;i<ne;++i) photar[i]=0.;

  E_or_l = param[0];
  redshift = param[1];
  v = 1.e5*param[2]; /* put into cm/s */

  doppler_v = sqrt((1.+v/ccc)/(1.-v/ccc));

  SPECBINS=0;
  specfile=fopen("mulext.qdp","r");
  for (k=1;k<=3;++k) fgets(line,LENGTH,specfile);    
  while (fgets(line,LENGTH,specfile) != NULL) ++SPECBINS;
  fclose(specfile);

  E_array=dvector_me(1,SPECBINS);
  E_bin=dvector_me(1,SPECBINS);
  E_spectrum=dvector_me(1,SPECBINS);

  specfile=fopen("mulext.qdp","r");
  j=1;
  for (k=1;k<=3;++k) fgets(line,LENGTH,specfile);    
  while (fgets(line,LENGTH,specfile) != NULL) {
    sscanf(line,"%lf%lf%lf",&energy_or_lambda /* keV or A*/,&bin /* keV or A*/,&spectrum /* photons/cm^2/s/keV  or  photons/cm^2/s/A */);
    
    if (E_or_l == 0) {
      energy=energy_or_lambda;
      bin=2.*bin;
    }
    else if (E_or_l == 1) {
      bin=lambdatokeV*(1./(energy_or_lambda-bin)-1./(energy_or_lambda+bin)); /* convert from A to keV */
      energy=lambdatokeV/energy_or_lambda;
    }

    E_array[j]=1000.*energy;
    E_bin[j]=1000.*bin;
    E_spectrum[j]=spectrum;
    ++j;
  }
  fclose(specfile);

  /*    E_redshift=(1.+redshift)*1000.*(ear[i]+ear[i+1])/2.; */
  j=1;
  i=0;
  earlo=1000.*(1.+redshift)*doppler_v*ear[i];
  earhi=1000.*(1.+redshift)*doppler_v*ear[i+1];
  earBIN=earhi-earlo;
  Elo=E_array[j]-E_bin[j]/2.;
  Ehi=E_array[j]+E_bin[j]/2.;
  while (earhi<=Elo && i<ne-1) {
    ++i;
    earlo=1000.*(1.+redshift)*doppler_v*ear[i];
    earhi=1000.*(1.+redshift)*doppler_v*ear[i+1];
    earBIN=earhi-earlo;
  }
  while (Ehi<earlo && j<SPECBINS) {
    ++j;
    Elo=E_array[j]-E_bin[j]/2.;
    Ehi=E_array[j]+E_bin[j]/2.;
  }
  while (j<SPECBINS && i<ne-1) {
    if (Elo<earlo) {
      if (Ehi-earlo<earBIN) Ewidth=Ehi-earlo;
      else Ewidth=earBIN;
    } else {
      if (earhi-Elo<E_bin[j]) Ewidth=earhi-Elo;
      else Ewidth=E_bin[j];
    }
    photar[i]+=Ewidth/earBIN*E_spectrum[j];
    if (earhi<Ehi) {
      ++i;
      earlo=1000.*(1.+redshift)*doppler_v*ear[i];
      earhi=1000.*(1.+redshift)*doppler_v*ear[i+1];
      earBIN=earhi-earlo;
    } else {
      ++j;
	Elo=E_array[j]-E_bin[j]/2.;
	Ehi=E_array[j]+E_bin[j]/2.;
    }
  }

  free_dvector_me(E_array,1,SPECBINS);
  free_dvector_me(E_bin,1,SPECBINS);
  free_dvector_me(E_spectrum,1,SPECBINS);

  return 0;
}

