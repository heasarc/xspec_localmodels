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

int siabs
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(siabs,SIABS,siabs,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))

/* START */

/* Forward declarations (converted from K&R style). */


/* FINISH */


char* FGMSTR(char* name);

/* Subroutines from PHOTOION */

/* SIABS Subroutines */
void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_edge_p[]) ;


int siabs
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

  
  

  double redshift,v_rad,sigmav_rad,N_e,**Nion,COLNORM;
  int verbose;

  int LOG;

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



  int *list;

  char name[]="PHOTOION_DIR";
  char *DATADIR;

  int ion_A,ion_z;

  double **H_excite, **He_excite;

  double *ABUND,*oshe;

  int klo,khi,jlow,jhigh;

  int DIST,lines=0;


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
      printf("MIABS: Failed to open %s\n", probepath);
      return 1;
    }
    fclose(probe);
  }

  /* Validate the ion before allocating. These returns used to sit ~40 lines
   * below the allocation block, so a bad ion_A or ion_z leaked everything
   * taken so far on every evaluation -- and a fit re-enters this path each
   * iteration. They need only param[], so they belong up here. */
  ion_A = (int) param[0];
  ion_z = (int) param[1];
  if (!(ion_A == 1 || ion_A == 2 || ion_A == 6 || ion_A == 7 || ion_A == 8 || ion_A == 10 || ion_A == 12 || ion_A == 13 || ion_A == 14 || ion_A == 16 || ion_A == 18 || ion_A == 20 || ion_A == 26 || ion_A == 28)) {
    printf("Selected element A = %2d not supported.\n",ion_A);
    return 0.;
  }
  if (ion_z>ion_A) {
    printf("NOTE: z > A not allowed.\n");
    return 0.;
  }

  /* FILE NAMES */



  ABUND=pion_dvector(1,30);
  pion_ad_abundance(ABUND);

  oshe=pion_dvector(1,30);
  pion_ad_oscillator_he(oshe);

  Nion=pion_dmatrix(1,28,1,28);
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) Nion[i][j]=0.;

  ion_A = (int) param[0];
  ion_z = (int) param[1];
  COLNORM=param[2];
  redshift = param[3];     /* cosmological redshift */
  v_rad = 1.e5*param[4];   /* velocity shift: convert [km/s] to [cm/s] */
  sigmav_rad = 1.e5*param[5];   /* gaussian (1-sigma) velocity width: convert [km/s] to [cm/s] */
  if (sigmav_rad>0.) lines=1; else if (sigmav_rad<=0.) lines=0;
  EMIN=1000.*param[6]; /* from keV to eV */
  EMAX=1000.*param[7]; /* from keV to eV */
  SPECBINS=(int) param[8];

  LOG=0;
  verbose=0;

  Nion[ion_A][ion_z]=COLNORM;


  if (LOG!=0) {
    for (i=1;i<=28;++i) {
      for (j=1;j<=28;++j) { 
	if (Nion[i][j]) {
	  Nion[i][j]=pow(10.,Nion[i][j]);
	}
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
  pion_free_dmatrix(Nion,1,28,1,28);
  pion_free_dvector(ionizsigmatemp,1,SPECBINS);
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
  return 0.;
}


/* For convolution: Normalized gaussian profile centered at 0. with s=sigma */

/* Photoexcitation cross section [cm^2] from ground.
   global: E0, DELTANUD, OSCILLATOR, ALPHA */

/* Differential oscillator strength in atomic units: 1/Hartree = 1/(2.*13.6 eV) */

/* Photoionization cross-section in atomic units */


/* ground-state photoionization cross-section from Verner */

/* ground-state photoionization cross-section from Verner */


