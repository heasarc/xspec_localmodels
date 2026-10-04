#include "cfortran.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "photoion_phys.h"
#include "photoion_atomdata.h"
#include "photoion_const.h"

#include "photoion_spline.h"
#include "photoion_state.h"

#include "photoion_alloc.h"

int vneutral
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(vneutral,VNEUTRAL,vneutral,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))

/* START */

/* Forward declarations (converted from K&R style). */


/* FINISH */


char* FGMSTR(char* name);

/* Subroutines from PHOTOION */

/* VNEUTRAL Subroutines */
void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_edge_p[]) ;


int vneutral
(float *ear,int ne,float *param,int ifl,float *photar,float *photer)
{
  /*ear[0] -> ear[ne]   in keV
    photar[0] -> photar[ne-1]  
    photar[i] = [photons/cm^2/s] for bin ear[i] -> ear[i+1]
    param[0] -> param[TOTAL-1]   gives XSPEC parameters from 1 to TOTAL*/

  const struct VERNER_STRUCT (*vernerionizsigma)[31];  /* cached, see photoion_atomdata.h */

  const struct VERNER_PARTIAL_STRUCT (*partialsigma)[125];
  const int *npartial;
  const struct PION_LINE_ROW *line_rows;   /* line.dat, cached */
  int nline_rows;

  
  

  double redshift,v_rad,sigmav_rad,N_e,**Nion;
  int verbose;

  double Elo,Ehi; /* Voigt function param.'s */
  double earBIN,Ewidth;

  /* junk values for strings, ints, and floats */
  double djunk;


  double ELO=1.e-3,EHI=1.e6;  /* MAKE SURE THIS RANGE IS OK!!! CHECK HERE!!! */


  int i,j,k,n;
  /* *.rr file */


  double earlo,earhi;
  
  int HIGHN=50;

  int ELEMENTS=12;


  int element,electron;

  int *list;

  char name[]="PHOTOION_DIR";
  char *DATADIR;

  double **H_excite, **He_excite;

  double *ABUND,*oshe;

  int klo,khi,jlow,jhigh;

  int DIST,lines;


  /* initialize photar array */
  for (i=0;i<ne;++i) photar[i]=0.;

  DATADIR=FGMSTR(name);
  pion_ad_begin(DATADIR);

  /* Probe for the atomic data before anything is allocated. The open below
   * is checked, but it returns with ~15 allocations already live, and this
   * path is re-entered on every evaluation of a misconfigured model. Message
   * and return value are unchanged; the check below stays as a guard. */
  {
    char probepath[1024];
    FILE *probe;
    snprintf(probepath,sizeof probepath,"%s/photoion_dat/abundance.dat",DATADIR);
    probe=fopen(probepath,"r");
    if (probe == NULL) {
      printf("VNEUTRAL: Failed to open %s\n", probepath);
      return 1;
    }
    fclose(probe);
  }

  /* FILE NAMES */


  Nion=pion_dmatrix(1,28,1,28);
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) Nion[i][j]=0.;


  ABUND=pion_dvector(1,30);
  pion_ad_abundance(ABUND);

  oshe=pion_dvector(1,30);
  pion_ad_oscillator_he(oshe);

  Nion[1][1]=param[0];
  /* No electron scattering is added through N_e: it was never assigned in
   * the original code, so it held whatever the stack did (in practice
   * about 0). It is now 0 by definition, so vneutral has no scattering
   * term, as before. */
  N_e=0.;
  sigmav_rad = 1.e5*param[1];   /* gaussian (1-sigma) velocity width: convert [km/s] to [cm/s] */
  djunk=param[2];  ABUND[2]*=djunk;  Nion[2][2]=ABUND[2]*Nion[1][1];
  djunk=param[3];  ABUND[6]*=djunk;  Nion[6][6]=ABUND[6]*Nion[1][1];
  djunk=param[4];  ABUND[7]*=djunk;  Nion[7][7]=ABUND[7]*Nion[1][1];
  djunk=param[5];  ABUND[8]*=djunk;  Nion[8][8]=ABUND[8]*Nion[1][1];
  djunk=param[6];  ABUND[10]*=djunk;  Nion[10][10]=ABUND[10]*Nion[1][1];
  djunk=param[7];  ABUND[12]*=djunk;  Nion[12][12]=ABUND[12]*Nion[1][1];
  djunk=param[8];  ABUND[13]*=djunk;  Nion[13][13]=ABUND[13]*Nion[1][1];
  djunk=param[9];  ABUND[14]*=djunk;  Nion[14][14]=ABUND[14]*Nion[1][1];
  djunk=param[10];  ABUND[16]*=djunk;  Nion[16][16]=ABUND[16]*Nion[1][1];
  djunk=param[11];  ABUND[18]*=djunk;  Nion[18][18]=ABUND[18]*Nion[1][1];
  djunk=param[12];  ABUND[20]*=djunk;  Nion[20][20]=ABUND[20]*Nion[1][1];
  djunk=param[13];  ABUND[26]*=djunk;  Nion[26][26]=ABUND[26]*Nion[1][1];
  djunk=param[14];  ABUND[28]*=djunk;  Nion[28][28]=ABUND[28]*Nion[1][1];
  redshift = param[15];     /* cosmological redshift */
  v_rad = 1.e5*param[16];   /* velocity shift: convert [km/s] to [cm/s] */
  EMIN=1000.*param[17]; /* convert keV to eV */
  EMAX=1000.*param[18]; /* convert keV to eV */
  SPECBINS=(int) param[19];

  verbose=0;

  if (sigmav_rad) lines=1;
  else lines=0;

  for (element=1;element<=28;++element) {
    for (electron=1;electron<=28;++electron) {
      if (Nion[element][electron]) {
        Nion[element][electron]=fabs(Nion[element][electron]);
	if (verbose) printf("Nion[%2d][%2d]=%e\n",element,electron,Nion[element][electron]);
      }
    }
  }

  /* SAME AS MIABS FROM HERE ON OUT? (except for moving ABUND up )*/

  /* To put EMIN and EMAX in source rest-frame energy units [eV]
     We convert back at the very end */  
  EMIN=EMIN*(1.+redshift);
  EMAX=EMAX*(1.+redshift);
 
  /* Decrease these to improve accuracy */
  voigt_lim=1.e-4;
  tau_lim=1.e-5;

  doppler_rad=pion_doppler(v_rad);

  EGRID=pion_dvector(1,GRIDNUM);  
  PIGRID=pion_dvector(1,GRIDNUM); 

  LOWE_EGRID=pion_dvector(1,LOWE_GRIDNUM);  
  LOWE_PIGRID=pion_dvector(1,LOWE_GRIDNUM); 

  E_array=pion_dvector(1,SPECBINS);       /* energy axis */
  tau=pion_dvector(1,SPECBINS);           /* total opacity in all ions */
  tau_exc=pion_dvector(1,SPECBINS);        
  tau_edge=pion_dvector(1,SPECBINS);       
  ionizsigmatemp=pion_dvector(1,SPECBINS);
  H_excite=pion_dmatrix(1,28,1,6);
  He_excite=pion_dmatrix(1,28,1,9);
  list=pion_ivector(1,ELEMENTS);

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
  if (verbose) if (EBIN>0.5) {
    printf("****************************************\n");
    printf("* Warning: EBIN = %6.2e eV is big!  *\n",EBIN);
    printf("* You might want to increase SPECBINS. *\n");
    printf("****************************************\n");
  }
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

  /* Verner total and partial photoionization fits, parsed once per process */
  vernerionizsigma=pion_ad_verner_full();
  partialsigma=pion_ad_verner_partial();
  npartial=pion_ad_verner_partial_count();

  /* Total electron Thomson depth */
  for (k=1;k<=SPECBINS;++k) {
    tau[k]+=N_e*sigmaT;
  }

  /* H- and He-like cross sections for C, N, and O */
  if (verbose) printf("H- and He-like edge cross sections for H,He,C to Ni...\n");
  pion_verner_edges(Nion,vernerionizsigma,partialsigma,npartial,26,tau_edge);

  if (lines) {
    if (verbose) printf("Determining LOW-n Photoexcitation Cross Sections & Opacity for C to Fe...\n");
    line_rows=pion_ad_line_rows(&nline_rows);
    pion_lown_line_opacity(Nion,sigmav_rad,line_rows,nline_rows,tau_exc,1.0,1.0,0);

    if (verbose) printf("Determining HIGH-n Photoexcitation Cross Sections & Opacity for C to Fe\n");
    pion_highn_line_opacity(Nion,sigmav_rad,pion_ad_highn(),list,ELEMENTS,HIGHN,oshe,tau_exc,1.0,1.0,0);
  }
    
  if (verbose) printf("L-shell ions:  Ne through Ni\n");
  pion_fac_shell_opacity(PION_FAC_L_NE_NI,Nion,sigmav_rad,1,lines,verbose,tau_edge,tau_exc,1.0,1.0,0);


  if (verbose) printf("L-shell ions:  C through O\n");
  pion_fac_shell_opacity(PION_FAC_L_C_O,Nion,sigmav_rad,1,lines,verbose,tau_edge,tau_exc,1.0,1.0,0);
  

  /* M-shell ions */
  if (verbose) printf("M-shell ions:  Mg through Ni\n");
  pion_fac_shell_opacity(PION_FAC_M,Nion,sigmav_rad,1,lines,verbose,tau_edge,tau_exc,1.0,1.0,0);
  

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
	djunk=EBIN/(E_array[k])*pion_gauss(sigmav_rad/ccc,(E_array[n]-E_array[k])/E_array[k])*tau_edge[n];
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
  pion_free_dvector(LOWE_EGRID,1,LOWE_GRIDNUM);  
  pion_free_dvector(LOWE_PIGRID,1,LOWE_GRIDNUM); 
  pion_free_dvector(ABUND,1,30);
  pion_free_dvector(oshe,1,30);
  pion_free_dvector(E_array,1,SPECBINS);       
  pion_free_dvector(tau,1,SPECBINS);           
  pion_free_dvector(tau_exc,1,SPECBINS);           
  pion_free_dvector(tau_edge,1,SPECBINS);           
  pion_free_dvector(ionizsigmatemp,1,SPECBINS);
  pion_free_dmatrix(Nion,1,28,1,28);
  pion_free_dmatrix(H_excite,1,28,1,6);
  pion_free_dmatrix(He_excite,1,28,1,9);
  pion_free_ivector(list,1,ELEMENTS);
  pion_free_dvector(EGRID,1,GRIDNUM);  
  pion_free_dvector(PIGRID,1,GRIDNUM); 
  if (verbose) printf("...done!\n");
  
  /* A data file could not be opened: everything above ran on zeros. Report it
   * the same way as the other input errors, after the normal cleanup. */
  if (pion_ad_failed()) {
    for (i=0;i<ne;++i) photar[i]=0.;
    return 1;
  }
  return 0;
}


/* For convolution: Normalized gaussian profile centered at 0. with s=sigma */

/* Photoexcitation cross section [cm^2] from ground.
   global: E0, DELTANUD, OSCILLATOR, ALPHA */

/* Differential oscillator strength in atomic units: 1/Hartree = 1/(2.*13.6 eV) */

/* Photoionization cross-section in atomic units */


/* ground-state photoionization cross-section from Verner */

/* ground-state photoionization cross-section from Verner */


