#include "cfortran.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "photoion_phys.h"
#include "photoion_atomdata.h"
#include "photoion_const.h"

#include "photoion_state.h"

#include "photoion_alloc.h"

int neutral
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(neutral,NEUTRAL,neutral,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))

/* START */

/* Forward declarations (converted from K&R style). */


/* FINISH */


char* FGMSTR(char* name);

/* Subroutines from PHOTOION */

/* NEUTRAL Subroutines */
void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_edge_p[]) ;

int neutral
(float *ear,int ne,float *param,int ifl,float *photar,float *photer)
{
  /*ear[0] -> ear[ne]   in keV
    photar[0] -> photar[ne-1]  
    photar[i] = [photons/cm^2/s] for bin ear[i] -> ear[i+1]
    param[0] -> param[TOTAL-1]   gives XSPEC parameters from 1 to TOTAL*/


  const struct PION_LINE_ROW *line_rows;   /* line.dat, cached */
  int nline_rows;
  int pathlen;     /* buffer size for any data path under DATADIR (open issue 2) */

  
  

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

  char *root,*element_name,*temp;
  const char *ext;   /* always a string literal, never owned */
  int i,j,k,n;
  /* *.en file */
  double energy;


  double earlo,earhi;
  
  int HIGHN=50;

  int ELEMENTS=12;

  double E0=0.,DELTANUD,ALPHA,OSCILLATOR;

  int element,electron;

  int *list;

  char name[]="PHOTOION_DIR";
  char *DATADIR;

  double **H_excite, **He_excite;

  double *ABUND,*oshe;

  int klo,khi,jlow,jhigh;

  int DIST,lines;

  FILE *input;

  /* initialize photar array */
  for (i=0;i<ne;++i) photar[i]=0.;

  DATADIR=FGMSTR(name);
  pion_ad_begin(DATADIR);
  pathlen=pion_ad_pathlen();

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
      printf("NEUTRAL: Failed to open %s\n", probepath);
      return 1;
    }
    fclose(probe);
  }

  /* FILE NAMES */

  root=malloc(pathlen);
  element_name=malloc(3);  /* 2-char symbols ("Ne") need 3 bytes with the NUL */
  temp=malloc(pathlen);
  sjunk=malloc(50);
  sjunk1=malloc(50);
  sjunk2=malloc(50);
  sjunk3=malloc(50);
  sjunk4=malloc(50);
  line=malloc(400);

  Nion=pion_dmatrix(1,28,1,28);
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) Nion[i][j]=0.;


  ABUND=pion_dvector(1,30);
  pion_ad_abundance(ABUND);

  oshe=pion_dvector(1,30);
  pion_ad_oscillator_he(oshe);

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

  doppler_rad=pion_doppler(v_rad);

  EGRID=pion_dvector(1,GRIDNUM);  
  PIGRID=pion_dvector(1,GRIDNUM); 
  PIGRID_2=pion_dvector(1,GRIDNUM);

  LOWE_EGRID=pion_dvector(1,LOWE_GRIDNUM);  
  LOWE_PIGRID=pion_dvector(1,LOWE_GRIDNUM); 
  LOWE_PIGRID_2=pion_dvector(1,LOWE_GRIDNUM);

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


  /* Total electron Thomson depth */
  for (k=1;k<=SPECBINS;++k) {
    tau[k]+=N_e*sigmaT;
  }

  if (lines) {
    if (verbose) printf("Determining LOW-n Photoexcitation Cross Sections & Opacity for C to Fe...\n");
    line_rows=pion_ad_line_rows(&nline_rows);
    pion_lown_line_opacity(Nion,sigmav_rad,line_rows,nline_rows,tau_exc,1.0,1.0,0);

    if (verbose) printf("Determining HIGH-n Photoexcitation Cross Sections & Opacity for C to Fe\n");
    pion_highn_line_opacity(Nion,sigmav_rad,pion_ad_highn(),list,ELEMENTS,HIGHN,oshe,tau_exc,1.0,1.0,0);
    
    
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
	  snprintf(root,pathlen,"%s/photoion_dat/L_shell/%s",DATADIR,element_name);
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  if (electron>=3) {
	    ext="tr_short"; 
	    if (electron<10) snprintf(temp,pathlen,"%s0%da.%s",root,electron,ext);
	    else snprintf(temp,pathlen,"%s%da.%s",root,electron,ext);
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
		pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc,1.0,1.0,0);
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
	  snprintf(root,pathlen,"%s/photoion_dat/L_shell/%s",DATADIR,element_name);
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  ext="tr_short"; 
	  if (electron<10) snprintf(temp,pathlen,"%s0%da.%s",root,electron,ext);
	  else snprintf(temp,pathlen,"%s%da.%s",root,electron,ext);
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
	      pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc,1.0,1.0,0);
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
	  snprintf(root,pathlen,"%s/photoion_dat/M_shell/%s",DATADIR,element_name);
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  
	  /* Transitions */
	  ext="tr_short"; 
	  snprintf(temp,pathlen,"%s%2da.%s",root,electron,ext);
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
	      pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc,1.0,1.0,0);
	    }
	  }
	  fclose(input);
	}      
      }
    }
  }
  
  /* Read in edge opacity */
  snprintf(temp,pathlen,"%s/photoion_dat/neutral.tau",DATADIR);
  input=fopen(temp,"r");
  j=1;

  E_bin=pion_dvector(1,SPECBINS);
  /* Open issue 11: the 3 header lines below are parsed as data, sscanf converts
   * nothing, and tau_edge[1..3] take whatever djunk last held. That used to be
   * the last value read from oscillator_he.dat, by an inline loop now in
   * photoion_atomdata.c. Keep it, so this refactor stays bit-for-bit; drop it
   * when issue 11 is fixed. */
  djunk=pion_ad_oscillator_he_lastraw();
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
  free(root);
  free(temp);
  free(sjunk);
  free(sjunk1);
  free(sjunk2);
  free(sjunk3);
  free(sjunk4);
  free(line);
  free(element_name);
  pion_free_dvector(LOWE_EGRID,1,LOWE_GRIDNUM);  
  pion_free_dvector(LOWE_PIGRID,1,LOWE_GRIDNUM); 
  pion_free_dvector(LOWE_PIGRID_2,1,LOWE_GRIDNUM);
  pion_free_dvector(ABUND,1,30);
  pion_free_dvector(oshe,1,30);
  pion_free_dvector(E_array,1,SPECBINS);       
  pion_free_dvector(tau,1,SPECBINS);           
  pion_free_dvector(tau_exc,1,SPECBINS);           
  pion_free_dvector(tau_edge,1,SPECBINS);           
  pion_free_dvector(ionizsigmatemp,1,SPECBINS);
  pion_free_dvector(E_bin,1,SPECBINS);       
  pion_free_dmatrix(Nion,1,28,1,28);
  pion_free_dmatrix(H_excite,1,28,1,6);
  pion_free_dmatrix(He_excite,1,28,1,9);
  pion_free_ivector(list,1,ELEMENTS);
  pion_free_dvector(EGRID,1,GRIDNUM);  
  pion_free_dvector(PIGRID,1,GRIDNUM); 
  pion_free_dvector(PIGRID_2,1,GRIDNUM);
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


