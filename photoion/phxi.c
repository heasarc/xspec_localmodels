#include "cfortran.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "photoion_phys.h"
#include "photoion_atomdata.h"
#include "photoion_emission.h"
#include "photoion_const.h"

#include "photoion_integrate.h"
#include "photoion_spline.h"
#include "photoion_state.h"

#include "photoion_alloc.h"

int phxi
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(phxi,PHXI,phxi,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))
#define TEMPERATURES (11)

/* START */

/* Forward declarations (converted from K&R style). */


/* FINISH */


char* FGMSTR(char* name);

/* PHXI Subroutines */

void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_edge_p[]) ;

double REC_spline_px(double temp);


int phxi
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
  const struct HIGHER_ORDER_STRUCT (*highn)[3];
  int pathlen;     /* buffer size for any data path under DATADIR (open issue 2) */

  
  
  
  

  double *type1_spectrum,*type2_spectrum,*type3_spectrum,*type4_spectrum,*type5_spectrum,*type6_spectrum,*type7_spectrum,*type8_spectrum;

  double int_ans;

  int LENGTH=400;

  double Elo,Ehi; /* Voigt function param.'s */
  double Ewidth,earBIN;

  int first=1;
  

  /* junk values for strings, ints, and floats */
  char *line;
  int ijunk;
  double djunk;


  double ELO=1.e-3,EHI=1.e6;  /* MAKE SURE THIS RANGE IS OK!!! CHECK HERE!!! */


  char *temp;
  int i,j,k,n;
  /* *.rr file */

  double *Labsorb,*specRR,*specDR,*spec,*convert;
  double *rec_spectrum_temp,*specRR_temp,*specDR_temp;
  double **rateRR,**rateDR,**ratePI;



  int verbose;

  double earlo,earhi;
  
  int HIGHN=50;

  int ELEMENTS=12;


  int element,electron;
  int type;


  int *list;

  char name[]="PHOTOION_DIR";
  char *DATADIR;

  double **H_excite, **He_excite;

  double redshift,**ion,N_H;
  double H0_cosmo=0.,Omega_m_cosmo=0.; /* from XSPEC, when D=0 */
  int fileincr;


  double FLUXAVE;

  double *ABUND,*oshe;

  int klo,khi,jlow,jhigh;
  double spectemp;

  int DIST,lines=0;

  FILE *E_specfile=NULL,*l_specfile=NULL;
  FILE *input,*E_output=NULL,*l_output=NULL;
  char *specfile_name;

  /* Initialize the model array */
  for (i=0;i<ne;++i) photar[i] = 0.;


  /* input.qdp is not read until ~900 lines below, by which point ~90 arrays are
   * live, around 40 of them SPECBINS-sized. The open there is checked, but it
   * returns without freeing any of them, so a missing file leaked tens of
   * megabytes per evaluation. Probe for it here, before anything is allocated.
   * The message and the return value are unchanged; the later check stays as a
   * guard. */
  if (param[20] > 0) {
    FILE *probe = fopen("input.qdp","r");
    if (probe == NULL) {
      printf("The file 'input.qdp' must exist in this directory.\n");
      return 0;
    }
    fclose(probe);
  }

  /* With D=0 and a non-zero redshift, the distance comes from the Hubble law,
   * using XSPEC's cosmology. Check it here, before anything is allocated. */
  if (param[28] == 0. && param[15] != 0.) {
    if (pion_cosmology(&H0_cosmo,&Omega_m_cosmo)) return 1;
  }

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
      printf("PHXI: Failed to open %s\n", probepath);
      return 1;
    }
    fclose(probe);
  }

  /* xi.dat is read a few hundred lines below, after the allocation block; the
   * open there is checked but returns without freeing. Probe here instead. */
  {
    FILE *probe = fopen("xi.dat","r");
    if (probe == NULL) {
      printf("The file 'xi.dat' must exist in this directory.\n");
      printf("See $DATADIR/photoion_dat/xi.dat for an example.\n");
      return 0;
    }
    fclose(probe);
  }

  /* FILE NAMES */
  specfile_name=malloc(200);
  temp=malloc(pathlen);

  line=malloc(400);

  Nion=pion_dmatrix(1,28,1,28);
  Tion=pion_dmatrix(1,28,1,28);
  EMion=pion_dmatrix(1,28,1,28);
  EM=pion_dmatrix(1,28,1,28);
  ratePI=pion_dmatrix(1,28,1,28);
  rateRR=pion_dmatrix(1,28,1,28);
  rateDR=pion_dmatrix(1,28,1,28);
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) {Nion[i][j]=0.;Tion[i][j]=10.;EMion[i][j]=0.;EM[i][j]=0.;ratePI[i][j]=0.;rateRR[i][j]=0.;rateDR[i][j]=0.;}


  pion_ad_temperature(Tion);

  oshe=pion_dvector(1,30);
  pion_ad_oscillator_he(oshe);

  ABUND=pion_dvector(1,30);
  pion_ad_abundance(ABUND);

  type = (int) param[0];
  N_H=param[1];
  djunk=param[2];  ABUND[2]*=djunk;
  djunk=param[3];  ABUND[6]*=djunk;
  djunk=param[4];  ABUND[7]*=djunk;
  djunk=param[5];  ABUND[8]*=djunk;
  djunk=param[6];  ABUND[10]*=djunk;
  djunk=param[7];  ABUND[12]*=djunk;
  djunk=param[8];  ABUND[13]*=djunk;
  djunk=param[9];  ABUND[14]*=djunk;
  djunk=param[10];  ABUND[16]*=djunk;
  djunk=param[11];  ABUND[18]*=djunk;
  djunk=param[12];  ABUND[20]*=djunk;
  djunk=param[13];  ABUND[26]*=djunk;
  djunk=param[14];  ABUND[28]*=djunk;
  redshift = param[15];
  v_rad = 1.e5*param[16];   /* convert [km/s] to [cm/s] */
  v_trans = 1.e5*param[17];   /* convert [km/s] to [cm/s] */
  sigmav_rad = 1.e5*param[18];   /* convert [km/s] to [cm/s] */
    if (sigmav_rad>0.) lines=1; else if (sigmav_rad<=0.) lines=0;
  sigmav_trans = 1.e5*param[19]; /* convert [km/s] to [cm/s] */
  INPUT = param[20];
  INPUT_SHIFT = param[21];
  GAMMA = param[22];
  L_EMIN = 1000.*param[23]; /* from keV to eV */
  L_EMAX = 1000.*param[24]; /* from keV to eV */
  L_X = 1.e30*param[25];  /* convert [1e30 ergs/s] to [ergs/s] */
  FLUXAVE=param[26];
  f_COVERING = param[27];
  D = parsectocm*param[28];    /* convert [pc] to [cm] */
  EMIN=1000.*param[29]; /* from keV to eV */
  EMAX=1000.*param[30]; /* from keV to eV */
  SPECBINS=(int) param[31];
  fileincr = (int) param[32];
  verbose =(int) param[33];

  snprintf(temp,pathlen,"xi.dat");
  input=fopen(temp,"r");
  if (input==NULL) {
    printf("The file 'xi.dat' must exist in this directory.\n");
    printf("See $DATADIR/photoion_dat/xi.dat for an example.\n");
    return 0;
  }
  FRACXINUM=0;
  while (fscanf(input,"%lf%lf",&djunk,&djunk)!=EOF) ++FRACXINUM;
  fclose(input);
  xi_frac_grid=pion_dvector(1,FRACXINUM);  
  frac_grid=pion_dvector(1,FRACXINUM);  
  
  input=fopen(temp,"r");
  for (i=1;i<=FRACXINUM;++i) fscanf(input,"%lf%lf",&(xi_frac_grid[i]),&(frac_grid[i]));
  fclose(input);

  snprintf(temp,pathlen,"%s/photoion_dat/xi_ions.dat",DATADIR);
  input=fopen(temp,"r");
  fscanf(input,"%d",&FIONXINUM);
  xi_fion_grid=pion_dvector(1,FIONXINUM);  
  fion_grid=pion_dvector(1,FIONXINUM);  
  fion_grid_2=pion_dvector(1,FIONXINUM);  
  ion=pion_dmatrix(1,29,1,FIONXINUM);

  xi_array=pion_dvector(1,XINUM);
  fion_array=pion_dvector(1,XINUM);
  fion_array_2=pion_dvector(1,XINUM);
  for (i=1;i<=XINUM;++i) xi_array[i]=XIMIN+(XIMAX-XIMIN)*((double) i-1)/((double) (XINUM-1));

  /* read in each element - where subscript for ion[2][i] = ROMAN numeral*/
  /* hydrogen */
  element=1;
  N_e=0.;
  HNORM=0.;
  for (i=1;i<=FIONXINUM;++i) {
    fscanf(input,"%d%lf",&ijunk,&(xi_fion_grid[i]));
    for (k=1;k<=element+1;++k) {
      fscanf(input,"%lf",&(ion[k][i]));
      ion[k][i]=fabs(ion[k][i]);
    }
  }
  for (k=1;k<=element+1;++k) {
    for (i=1;i<=FIONXINUM;++i) fion_grid[i]=ion[k][i];
    pion_spline(xi_fion_grid,fion_grid,FIONXINUM,1.e40,1.e40,fion_grid_2);
    for (i=1;i<=XINUM;++i) fion_array[i]=pion_frac(xi_array[i])*pion_fion(xi_array[i]);
    pion_spline(xi_array,fion_array,XINUM,1.e40,1.e40,fion_array_2);

    electron=element-k+1;
    if (electron>=1) Nion[element][electron]=pion_integrate(pion_fion_integrand,XIMIN,XIMAX);
    if (electron==0) N_e+=pion_integrate(pion_fion_integrand,XIMIN,XIMAX);
    HNORM+=pion_integrate(pion_fion_integrand,XIMIN,XIMAX);
  }
  Nion[1][1]=Nion[1][1]*N_H/HNORM;
  N_e=N_e*N_H/HNORM;
  if (verbose) printf("Nion[%2d][%2d]=%e\n",element,element,Nion[element][element]);
  
  /* All other elements: He, C, N, O, Ne, Mg, Si, S, Ar, Ca, Fe */
  /* Aluminum and Nickel not calculated by xstar */
  list=pion_ivector(1,12);
  list[1]=1;list[2]=2; list[3]=6; list[4]=7; list[5]=8; list[6]=10; list[7]=12; list[8]=14; list[9]=16; list[10]=18; list[11]=20; list[12]=26;
  for (j=2;j<=12;++j) {
    element=list[j];
    for (i=1;i<=FIONXINUM;++i) {
      fscanf(input,"%d%lf",&ijunk,&(xi_fion_grid[i]));
      for (k=1;k<=element+1;++k) {
	fscanf(input,"%lf",&(ion[k][i]));
	ion[k][i]=fabs(ion[k][i]);
      }
    }
    for (k=1;k<=element+1;++k) {
      fion_integrate=0;
      for (i=1;i<=FIONXINUM;++i) fion_grid[i]=ion[k][i];
      pion_spline(xi_fion_grid,fion_grid,FIONXINUM,1.e40,1.e40,fion_grid_2);
	/*      printf("%e  %e  %e  %e\n",xi_fion_grid[i],pion_frac(xi_fion_grid[i]),fion_grid[i],ion[k][i]);*/
      
      fion_integrate=0;
      for (i=1;i<=XINUM;++i) {
	fion_array[i]=pion_frac(xi_array[i])*pion_fion(xi_array[i]);
	if (fion_array[i]!=0.) fion_integrate=1;
      }
      if (fion_integrate) {
	pion_spline(xi_array,fion_array,XINUM,1.e40,1.e40,fion_array_2);
	electron=element-k+1;
	if (electron>=1) {
	  Nion[element][electron]=ABUND[element]*N_H/HNORM*pion_integrate(pion_fion_integrand,XIMIN,XIMAX);
	}
	if (electron>=1) N_e+=((double) (element-electron))*Nion[element][electron];
	else N_e+=((double) (element-electron))*ABUND[element]*N_H/HNORM*pion_integrate(pion_fion_integrand,XIMIN,XIMAX);
      }
    }
  }

  for (element=2;element<=28;++element) {
    for (electron=1;electron<=28;++electron) {
      if (Nion[element][electron]) {
        Nion[element][electron]=fabs(Nion[element][electron]);
	if (verbose) printf("Nion[%2d][%2d]=%e\n",element,electron,Nion[element][electron]);
      }
    }
  }
  if (verbose) printf("N_e=%e\n",N_e);

  fclose(input);

  pion_free_ivector(list,1,10);
  pion_free_dvector(xi_frac_grid,1,FRACXINUM);
  pion_free_dvector(frac_grid,1,FRACXINUM);
  pion_free_dmatrix(ion,1,29,1,FIONXINUM);
  pion_free_dvector(xi_fion_grid,1,FIONXINUM);
  pion_free_dvector(fion_grid,1,FIONXINUM);
  pion_free_dvector(fion_grid_2,1,FIONXINUM);
  pion_free_dvector(xi_array,1,XINUM);
  pion_free_dvector(fion_array,1,XINUM);
  pion_free_dvector(fion_array_2,1,XINUM);

  if (type == 4 || type == 5) {v_trans=v_rad; sigmav_trans=sigmav_rad;}


  /* SAME AS "PHOTOION" FROM HERE ON OUT (aside from pion_fion_integrand, pion_frac,pion_fion) */

  /* To put EMIN and EMAX in source rest-frame energy units
     We convert back at the very end */  
  EMIN=EMIN*(1.+redshift);
  EMAX=EMAX*(1.+redshift);
 
  voigt_lim=1.e-4;
  tau_lim=1.e-4;
  pi_rate_lim=1.e-3;

  /*  if (type<=0) tau_lim=0.;*/

  if (D==0.) {
    if (redshift!=0.) {
      /* hubble equation integral */
      z_array=pion_dvector(1,HUBBLE_BINS);
      hubble_array=pion_dvector(1,HUBBLE_BINS);
      hubble_array_2=pion_dvector(1,HUBBLE_BINS);
      for (k=1;k<=HUBBLE_BINS;++k) {
	z_array[k]=((double) k-1.)/((double) HUBBLE_BINS-1.)*redshift;
	hubble_array[k]=pion_hubble_integrand(z_array[k],H0_cosmo,Omega_m_cosmo);
      }
      pion_spline(z_array,hubble_array,HUBBLE_BINS,1.e40,1.e40,hubble_array_2);
      D=pion_integrate(pion_hubble_integrate,0.,redshift); /* if D=0, use Hubble law */
      if (verbose) printf("redshift = %e   D = %e Mpc  (using H_0=%g km/s/Mpc, Omega_m=%g, Omega_lambda=%g, from XSPEC's cosmo setting)\n",redshift,D/parsectocm/1.e6,H0_cosmo,Omega_m_cosmo,1.-Omega_m_cosmo);
      pion_free_dvector(z_array,1,HUBBLE_BINS);
      pion_free_dvector(hubble_array,1,HUBBLE_BINS);
      pion_free_dvector(hubble_array_2,1,HUBBLE_BINS);
    } else {
      D=1.e6*parsectocm;
      if (verbose) printf("Since redshift=0 and D=0, taking D=1 Mpc\n");
    }
  }

  doppler_rad=pion_doppler(v_rad);
  doppler_trans=pion_doppler(v_trans);

  EGRID=pion_dvector(1,GRIDNUM);  
  PIGRID=pion_dvector(1,GRIDNUM); 
  PIGRID_2=pion_dvector(1,GRIDNUM);

  LOWE_EGRID=pion_dvector(1,LOWE_GRIDNUM);  
  LOWE_PIGRID=pion_dvector(1,LOWE_GRIDNUM); 
  LOWE_PIGRID_2=pion_dvector(1,LOWE_GRIDNUM);

  specRR=pion_dvector(1,SPECBINS);  
  specDR=pion_dvector(1,SPECBINS);  
  specRR_temp=pion_dvector(1,SPECBINS);  
  specDR_temp=pion_dvector(1,SPECBINS);  
  spec=pion_dvector(1,SPECBINS);  

  L_kT=pion_dvector(1,RECNUM);  
  L_RR=pion_dvector(1,RECNUM);  
  L_RR_2=pion_dvector(1,RECNUM);  
  L_DR=pion_dvector(1,RECNUM);  
  L_DR_2=pion_dvector(1,RECNUM);  
  L_REC=pion_dvector(1,RECNUM);  
  L_REC_2=pion_dvector(1,RECNUM);  
  RR_line=pion_dvector(1,RECNUM);  
  RR_line_2=pion_dvector(1,RECNUM);  
  DR_line=pion_dvector(1,RECNUM);  
  DR_line_2=pion_dvector(1,RECNUM);  

  E_array=pion_dvector(1,SPECBINS);       /* energy axis */
  E_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type1_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type2_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type3_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type4_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type5_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type6_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type7_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  type8_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (energy units) */
  l_array=pion_dvector(1,SPECBINS);       /* wavelength axis */
  l_spectrum=pion_dvector(1,SPECBINS);    /* spectrum (wavelength units) */
  convert=pion_dvector(1,SPECBINS);    /* spectrum (wavelength units) */
  abs_spectrum=pion_dvector(1,SPECBINS);    /* pure absorption spectrum */
  exc_spectrum=pion_dvector(1,SPECBINS);    /* photoexcitation spectrum */
  rec_spectrum=pion_dvector(1,SPECBINS);    /* recombination spectrum */
  rec_spectrum_temp=pion_dvector(1,SPECBINS);    /* recombination spectrum */
  tau=pion_dvector(1,SPECBINS);           /* total opacity in all ions */
  tau_edge=pion_dvector(1,SPECBINS);           /* total opacity in all ions */
  tau_exc=pion_dvector(1,SPECBINS);           /* total opacity in all ions */
  /* Zero the opacity arrays: pion_dvector() uses malloc, and every opacity
     routine accumulates onto these with +=, so stale values from a previous
     model evaluation would otherwise carry over. */
  for (k=1;k<=SPECBINS;++k) {tau[k]=0.; tau_edge[k]=0.; tau_exc[k]=0.;}
  Labsorb=pion_dvector(1,SPECBINS);    
  int_array=pion_dvector(1,SPECBINS);  
  int_array_2=pion_dvector(1,SPECBINS); 
  Hionizsigmaconv=pion_dmatrix(1,28,1,SPECBINS);
  Heionizsigmaconv=pion_dmatrix(1,28,1,SPECBINS);
  ionizsigmatemp=pion_dvector(1,SPECBINS);
  spectrumtemp=pion_dvector(1,SPECBINS);
  H_excite=pion_dmatrix(1,28,1,6);
  He_excite=pion_dmatrix(1,28,1,9);
  list=pion_ivector(1,ELEMENTS);
  Tvec=pion_dvector(1,TEMPERATURES);
  Yvec=pion_dvector(1,TEMPERATURES);
  Yvec2=pion_dvector(1,TEMPERATURES);

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
    l_array[k]=HC_KEV_ANGSTROM/(E_array[k]/1000.);  /* [Angstrom] */
  }

  /* for velocity convolution out to >= 4 sigma */
  DIST=(int) (4.*sigmav_rad/ccc*EMAX/EBIN); /*  +/-4 sigma */
  /*  if (verbose) printf("4-sigma_v^rad = %d bins\n",DIST);*/
  
  /* logarithmic energy grid for Fe-L RR and PI cross sections*/
  for (k=1;k<=GRIDNUM;++k) {
    /*EGRID=log10(energy)*/
    EGRID[k]=((double) k-1.)/((double) GRIDNUM)*log10(EHI/ELO)+log10(ELO);
  }

  /* Verner total and partial photoionization fits, parsed once per process */
  vernerionizsigma=pion_ad_verner_full();
  partialsigma=pion_ad_verner_partial();
  npartial=pion_ad_verner_partial_count();


  /* H- and He-like cross sections for C, N, and O */
  if (verbose) printf("H- and He-like edge cross sections for H,He,C to Ni...\n");
  pion_verner_edges(Nion,vernerionizsigma,partialsigma,npartial,28,tau_edge);
  if (verbose) printf("Determining LOW-n Photoexcitation Cross Sections & Opacity for C to Ni...\n");
  highn=pion_ad_highn();
  line_rows=pion_ad_line_rows(&nline_rows);
  if (lines) pion_lown_line_opacity(Nion,sigmav_rad,line_rows,nline_rows,tau_exc,0.9,1.1,1);

  if (lines) {
    if (verbose) printf("Determining HIGH-n Photoexcitation Cross Sections & Opacity for C to Fe\n");
    pion_highn_line_opacity(Nion,sigmav_rad,highn,list,ELEMENTS,HIGHN,oshe,tau_exc,0.9,1.1,1);
  }

  if (verbose) printf("L-shell ions:  Ne through Ni\n");
  pion_fac_shell_opacity(PION_FAC_L_NE_NI,Nion,sigmav_rad,1,lines,verbose,tau_edge,tau_exc,0.9,1.1,1);


  if (verbose) printf("L-shell ions:  C through O\n");
  pion_fac_shell_opacity(PION_FAC_L_C_O,Nion,sigmav_rad,1,lines,verbose,tau_edge,tau_exc,0.9,1.1,1);


  /* M-shell ions */
  if (verbose) printf("M-shell ions:  Mg through Ni\n");
  pion_fac_shell_opacity(PION_FAC_M,Nion,sigmav_rad,1,lines,verbose,tau_edge,tau_exc,0.9,1.1,1);
  

  if (0 && sigmav_rad) {
    if (verbose) printf("Convolving spectrum with appropriate velocity distribution.\n");
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
  
  /**************************************************************************/
  /* Determine integrand term for rate integrals, Labsorb, and abs_spectrum */

  /* Get rid of line excitation depth for high radial velocity width limit */
  if (type==7 || type==8) for (k=1;k<=SPECBINS;++k) tau[k]=tau[k]-tau_exc[k];

  /* Total electron Thomson depth */
  for (k=1;k<=SPECBINS;++k) {
    tau[k]+=N_e*sigmaT;
    tau_exc[k]+=N_e*sigmaT;
  }

  if (INPUT==0) {
    /* Normalize intrinsic power law */
    if (GAMMA != 2.) LNORM=L_X*ergstoeV*(2.-GAMMA)/(pow(L_EMAX,2.-GAMMA)-pow(L_EMIN,2.-GAMMA)); 
    else LNORM=L_X*ergstoeV/log(L_EMAX/L_EMIN);
    for (k=1;k<=SPECBINS;++k) abs_spectrum[k]=pion_L(E_array[k])*exp(-tau[k]);
    if (type==4) { /* Filled cone - rec. cont. (lower limit) */
      for (k=1;k<=SPECBINS;++k)	{
	if (tau[k]>1.e-5) {
	  Labsorb[k]=pion_L(E_array[k])*exp(-tau[k]);
	  abs_spectrum[k]=pion_L(E_array[k])*exp(-tau[k]);
	} else {
	  Labsorb[k]=pion_L(E_array[k]);
	  abs_spectrum[k]=pion_L(E_array[k]);
	}
      }
    } else if (type==5) { /* Filled cone - rec. cont. (upper limit) */
      for (k=1;k<=SPECBINS;++k) {
	if (tau[k]>1.e-5) {
	  Labsorb[k]=pion_L(E_array[k]);
	  abs_spectrum[k]=pion_L(E_array[k])*exp(-tau[k]);
	} else {
	  Labsorb[k]=pion_L(E_array[k]);
	  abs_spectrum[k]=pion_L(E_array[k]);
	}
      }
    } else {
      for (k=1;k<=SPECBINS;++k) {
	if (tau[k]>1.e-5) {
	  Labsorb[k]=pion_L(E_array[k])*(1.-exp(-tau[k]))/tau[k];
	  abs_spectrum[k]=pion_L(E_array[k])*exp(-tau[k]);
	} else {
	  Labsorb[k]=pion_L(E_array[k]);
	  abs_spectrum[k]=pion_L(E_array[k]);
	}
      }
    }
  } else if (INPUT>0) {
    snprintf(temp,pathlen,"input.qdp");
    input=fopen(temp,"r");
    if (input==NULL) {
      printf("The file 'input.qdp' must exist in this directory.\n");
      return 0;
    }
    for (k=1;k<=3;++k) fgets(line,LENGTH,input);
    /* Get # of data lines in input file */
    INPUT_SIZE=0;
    while (fgets(line,LENGTH,input) != NULL) ++INPUT_SIZE;
    fclose(input);
    E_input=pion_dvector(1,INPUT_SIZE);
    L_input=pion_dvector(1,INPUT_SIZE);
    L_input_2=pion_dvector(1,INPUT_SIZE);
    EtimesL_input=pion_dvector(1,INPUT_SIZE);
    EtimesL_input_2=pion_dvector(1,INPUT_SIZE);
    input=fopen(temp,"r");
    for (k=1;k<=3;++k) fgets(line,LENGTH,input);
    k=1;
    while (fgets(line,LENGTH,input) != NULL) {
      sscanf(line,"%lf%lf%lf",&(E_input[k])/* keV */,&HALFBIN_SIZE/* keV */,&(L_input[k])/* Flux units!!! [ph/cm^2/s/keV] */);
      if (INPUT_SHIFT) {E_input[k]=1000.*(1.+redshift)*E_input[k]; /* to convert from keV to eV */
      } else {E_input[k]=1000.*E_input[k];}
      HALFBIN_SIZE=1000.*HALFBIN_SIZE;
      L_input[k]=4.*PI*sqr(D)*L_input[k]/1000.;
      EtimesL_input[k]=E_input[k]*L_input[k];
      if (k==1) L_EMIN=E_input[1]-HALFBIN_SIZE;
      ++k;
    }
    fclose(input);
    L_EMAX=E_input[INPUT_SIZE]+HALFBIN_SIZE;
    pion_spline(E_input,EtimesL_input,INPUT_SIZE,1.e40,1.e40,EtimesL_input_2);
    LinterpNORM=1.;
    if (L_X>0.) {
      djunk=pion_integrate(pion_EtimesL,1.001*L_EMIN,0.999*L_EMAX);
      LinterpNORM=L_X*ergstoeV/djunk;
      if (verbose) printf("LinterpNORM = %e\n",LinterpNORM);
    }
    pion_spline(E_input,L_input,INPUT_SIZE,1.e40,1.e40,L_input_2);
    for (k=1;k<=SPECBINS;++k) abs_spectrum[k]=pion_Linterp(E_array[k])*exp(-tau[k]);
    if (type==4) { /* Filled cone - rec. cont. (lower limit) */
      for (k=1;k<=SPECBINS;++k)	{
	if (tau[k]>1.e-5) {
	  Labsorb[k]=pion_Linterp(E_array[k])*exp(-tau[k]);
	  abs_spectrum[k]=pion_Linterp(E_array[k])*exp(-tau[k]);
	} else {
	  Labsorb[k]=pion_Linterp(E_array[k]);
	  abs_spectrum[k]=pion_Linterp(E_array[k]);
	}
      }
    } else if (type==5) { /* Filled cone - rec. cont. (upper limit) */
      for (k=1;k<=SPECBINS;++k) {
	if (tau[k]>1.e-5) {
	  Labsorb[k]=pion_Linterp(E_array[k]);
	  abs_spectrum[k]=pion_Linterp(E_array[k])*exp(-tau[k]);
	} else {
	  Labsorb[k]=pion_Linterp(E_array[k]);
	  abs_spectrum[k]=pion_Linterp(E_array[k]);
	}
      }
    } else {
      for (k=1;k<=SPECBINS;++k) {
	if (tau[k]>1.e-5) {
	  Labsorb[k]=pion_Linterp(E_array[k])*(1.-exp(-tau[k]))/tau[k];
	  abs_spectrum[k]=pion_Linterp(E_array[k])*exp(-tau[k]);
	} else {
	  Labsorb[k]=pion_Linterp(E_array[k]);
	  abs_spectrum[k]=pion_Linterp(E_array[k]);
	}
      }
    }
  }

  if (type>1) {
    pion_emission_reemission(Nion,N_e,sigmav_rad,lines,verbose,list,ELEMENTS,HIGHN,oshe,ABUND,Labsorb,
                             H_excite,He_excite,ratePI,rateRR,rateDR,specRR,specDR);
  }

  if (type>1 && sigmav_trans) {
    if (verbose) printf("Convolving spectrum with appropriate velocity distribution.\n");
    DIST=(int) (4.*sigmav_trans/ccc*EMAX/EBIN); /*  +/-4 sigma */
    for (k=1;k<=SPECBINS;++k) {
    /* radiative decay after photoexcitation */
      spectrumtemp[k]=exc_spectrum[k]; 
      exc_spectrum[k]=0.;
    /* recombination */
      rec_spectrum_temp[k]=rec_spectrum[k];
      rec_spectrum[k]=0.;
      specRR_temp[k]=specRR[k];
      specRR[k]=0.;
      specDR_temp[k]=specDR[k];
      specDR[k]=0.;
    }
    klo=1;
    khi=SPECBINS;
    for (k=klo;k<=khi;++k) {
      if (k<=DIST) {jlow = 1;      jhigh = k+DIST;} 
      else         {jlow = k-DIST; jhigh = k+DIST;}
      if (jlow<1) jlow=1;
      if (jhigh>SPECBINS) jhigh=SPECBINS;
      for (n=jlow;n<=jhigh;++n) {
	/* radiative decay after photoexcitation */
	spectemp=EBIN/(E_array[k])*pion_gauss(sigmav_trans/ccc,(E_array[n]-E_array[k])/E_array[k])*spectrumtemp[n];
	exc_spectrum[k]+=spectemp;
	/* recombination */
	djunk=EBIN/(E_array[k])*pion_gauss(sigmav_trans/ccc,(E_array[n]-E_array[k])/E_array[k]);
	rec_spectrum[k]+=djunk*rec_spectrum_temp[n];
	specRR[k]+=djunk*specRR_temp[n];
	specDR[k]+=djunk*specDR_temp[n];
      }
    }
  }

  if (verbose) printf("Calculating final spectrum...\n");
  /* Output different spectral "types" */
  /* FLUXAVE is the factor times the observed continuum flux that yields the "average flux"; this modifies Sy2-echo contribution accordingly */
  for (k=1;k<=SPECBINS;++k) {
    /* Seyfert 1 - pure absorption */
    type1_spectrum[k]=abs_spectrum[k];
    /* Seyfert 2 */
    type2_spectrum[k]=FLUXAVE*rec_spectrum[k]+FLUXAVE*exc_spectrum[k];
    /* Pure recombination */
    type3_spectrum[k]=FLUXAVE*rec_spectrum[k];
    /* Seyfert 1 + reemission (lower limit = 4, upper limit = 5) */
    if (tau[k]>1.e-6) {
      type4_spectrum[k]=abs_spectrum[k]+f_COVERING*FLUXAVE*abs_spectrum[k]*tau_exc[k]+FLUXAVE*rec_spectrum[k]*(1.-exp(-tau[k]))/tau[k];
      type5_spectrum[k]=abs_spectrum[k]+f_COVERING*FLUXAVE*abs_spectrum[k]*tau_exc[k]+FLUXAVE*rec_spectrum[k]*(1.-exp(-tau[k]))/tau[k];
    } else {
      type4_spectrum[k]=abs_spectrum[k]+f_COVERING*FLUXAVE*abs_spectrum[k]*tau_exc[k]+FLUXAVE*rec_spectrum[k];
      type5_spectrum[k]=abs_spectrum[k]+f_COVERING*FLUXAVE*abs_spectrum[k]*tau_exc[k]+FLUXAVE*rec_spectrum[k];
    }
    if (INPUT==0) {
      type6_spectrum[k]=pion_L(E_array[k])+FLUXAVE*type2_spectrum[k];
    }
    else if (INPUT>0) {
      type6_spectrum[k]=pion_Linterp(E_array[k])+FLUXAVE*type2_spectrum[k];
    }
  }
  /* Sy1 - Pure Absorption */
  if (type<=1) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type1_spectrum[k];
  /* Sy2 - Photoionization and Photoexcitation */
  if (type==2) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type2_spectrum[k];
  /* Sy2 - Pure Recombination */
  if (type==3) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type3_spectrum[k];
  /* Sy1 - Patchy cone */
  if (type==4) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type4_spectrum[k];
  /* Sy1 - Filled cone */
  if (type==5) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type5_spectrum[k];
  /* CV - Sy2, except with unobscured intrinsic continuum */
  if (type==6) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type6_spectrum[k];
  /* CV - Sy2, completely unsaturated lines + unobscured intrinsic continuum */
  if (type==7) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type6_spectrum[k];
  /* CV - Sy2, completely unsaturated lines + no intrinsic continuum */
  if (type==8) for (k=1;k<=SPECBINS;++k) E_spectrum[k]=type2_spectrum[k];
  if (fileincr >= 0) {
    sprintf(specfile_name,"l_output_%d.qdp",fileincr);
    l_output=fopen(specfile_name,"w");
    sprintf(specfile_name,"E_output_%d.qdp",fileincr);
    E_output=fopen(specfile_name,"w");
    sprintf(specfile_name,"E_spectrum_%d.qdp",fileincr);
    E_specfile=fopen(specfile_name,"w");
    sprintf(specfile_name,"l_spectrum_%d.qdp",fileincr);
    l_specfile=fopen(specfile_name,"w");
    fprintf(E_specfile,"\n\n\n");
    fprintf(l_specfile,"\n\n\n");
  }
  for (k=1;k<=SPECBINS;++k)  {
    /* Convert object rest-frame luminosity to observed flux using redshift */
    E_spectrum[k]=E_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type1_spectrum[k]=type1_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type2_spectrum[k]=type2_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type3_spectrum[k]=type3_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type4_spectrum[k]=type4_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type5_spectrum[k]=type5_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type6_spectrum[k]=type6_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type7_spectrum[k]=type7_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    type8_spectrum[k]=type8_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    exc_spectrum[k]=exc_spectrum[k]/4./PI/sqr(D)/(1.+redshift);
    specRR[k]=specRR[k]/4./PI/sqr(D)/(1.+redshift);
    specDR[k]=specDR[k]/4./PI/sqr(D)/(1.+redshift);
    convert[k]=(1.e3*HC_KEV_ANGSTROM)/(sqr(l_array[k]));
    l_spectrum[k]=E_spectrum[k]*convert[k];
    if (fileincr >= 0) {
      if (INPUT==0) {
	fprintf(E_output,"%e   %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",E_array[k]/(1.+redshift),tau[k],pion_L(E_array[k])/4./PI/sqr(D),type1_spectrum[k],type2_spectrum[k],type3_spectrum[k],type4_spectrum[k],type5_spectrum[k],type6_spectrum[k],exc_spectrum[k],specRR[k],specDR[k]);
	fprintf(l_output,"%e   %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",(1.+redshift)*l_array[k],tau[k],pion_L(E_array[k])*convert[k]/4./PI/sqr(D),type1_spectrum[k]*convert[k],type2_spectrum[k]*convert[k],type3_spectrum[k]*convert[k],type4_spectrum[k]*convert[k],type5_spectrum[k]*convert[k],type6_spectrum[k]*convert[k],exc_spectrum[k]*convert[k],specRR[k]*convert[k],specDR[k]*convert[k]);
      } else if (INPUT>0) {
	fprintf(E_output,"%e   %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",E_array[k]/(1.+redshift),tau[k],pion_Linterp(E_array[k])/4./PI/sqr(D),type1_spectrum[k],type2_spectrum[k],type3_spectrum[k],type4_spectrum[k],type5_spectrum[k],type6_spectrum[k],exc_spectrum[k]/4./PI/sqr(D),specRR[k],specDR[k]);
	fprintf(l_output,"%e   %e  %e  %e  %e  %e  %e  %e  %e  %e  %e  %e\n",(1.+redshift)*l_array[k],tau[k],pion_Linterp(E_array[k])*convert[k]/4./PI/sqr(D),type1_spectrum[k]*convert[k],type2_spectrum[k]*convert[k],type3_spectrum[k]*convert[k],type4_spectrum[k]*convert[k],type5_spectrum[k]*convert[k],type6_spectrum[k]*convert[k],exc_spectrum[k]*convert[k]/4./PI/sqr(D),specRR[k]*convert[k],specDR[k]*convert[k]);
      }
      fprintf(E_specfile,"%e %e %e\n",E_array[k]/(1.+redshift)/1000.,EBIN/2./1000.,1000.*E_spectrum[k]);
      if (k==1) {    /* to avoid nonsense for first bin size */
	fprintf(l_specfile,"%e %e %e\n",(1.+redshift)*l_array[k],(l_array[k]-l_array[k+1])/2.,l_spectrum[k]);
      } else {
	fprintf(l_specfile,"%e %e %e\n",(1.+redshift)*l_array[k],(l_array[k-1]-l_array[k])/2.,l_spectrum[k]);
      }
    }
  }
  if (fileincr >= 0) {
    fclose(E_output);
    fclose(l_output);
    fclose(E_specfile);
    fclose(l_specfile);
  }    


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
    if (type==-1) photar[i]+=Ewidth/earBIN*tau[j]*EBIN/1000.*sqr(l_array[j])/HC_KEV_ANGSTROM;/* ph/cm^2/s in bin */
    else if (type==0) photar[i]+=Ewidth/earBIN*tau[j]*EBIN/1000.;/* ph/cm^2/s in bin */
    else photar[i]+=Ewidth*E_spectrum[j]*(1.+redshift);/* ph/cm^2/s in bin */
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
  
  first=1;
  if (verbose && type>1) {
    first=1;
    for (element=1;element<=28;++element) {
      for (electron=1;electron<=28;++electron) {
	if (EM[element][electron]) {
	  if (first) { 
	    printf("******************************************************************\n");
	    printf("*  Z  z  T[eV] PI=RR[1/s] EMion[cm^-3] Abundance f_ion*EM[cm^-3] *\n");
	    first=0;
	  }
	  printf("* %2d %2d %6.2lf  %4.2e    %4.2e    %4.2e     %4.2e    *\n",element,electron,Tion[element][electron],ratePI[element][electron],EMion[element][electron],ABUND[element],EM[element][electron]);
	}
      }
    }
    if (first==0) {
      printf("******************************************************************\n");
    }
  }

  if (INPUT==0) {
    for (i=1;i<=SPECBINS;++i) {
      int_array[i]=E_array[i]*eVtoergs/ccc*pion_L(E_array[i])*(1.-exp(-tau[i]));
    }
  } else {
    for (i=1;i<=SPECBINS;++i) {
      int_array[i]=E_array[i]*eVtoergs/ccc*pion_Linterp(E_array[i])*(1.-exp(-tau[i]));
    }
  }
  pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
  int_ans=f_COVERING*pion_integrate(pion_integrand,1.0001*EMIN,0.9999*EMAX);
  if (first==1) printf("******************************************************************\n");
  printf("* Radiation Pressure = %4.2e dyne (%7.3lf - %7.3lf keV, rf) *\n",int_ans,EMIN/1000.,EMAX/1000.);
  if (INPUT==0 && verbose) printf("* Power-law Norm at 1 keV: %e [photons/cm^2/s/keV]     *\n", LNORM/(1.+redshift)/1000./(4.*PI*sqr(D))*pow(1000.,2.-GAMMA));
  printf("******************************************************************\n");
    
  
  /* Free all the memory */
  free(temp);
  free(line);
  /* These were freed inside the conditional block that used them -- `if
   * (lines)`, `if (type>1)`, `if (fileincr >= 0)` -- while the malloc at the
   * top is unconditional, so each leaked whenever its branch was skipped.
   * `lines` is 0 at default parameters. Freed here instead. */
  free(specfile_name);
  /* Guard must match the allocating branch, which is `else if (INPUT > 0)`.
   * With `!= 0` a negative INPUT -- inside the declared parameter range -- frees
   * five arrays it never allocated. The first evaluation survives on the zeroed
   * globals; a later one, after an INPUT > 0 call left them dangling, aborts
   * XSPEC on a double free. */
  if (INPUT>0) {
    pion_free_dvector(E_input,1,INPUT_SIZE);
    pion_free_dvector(L_input,1,INPUT_SIZE);
    pion_free_dvector(L_input_2,1,INPUT_SIZE);
    pion_free_dvector(EtimesL_input,1,INPUT_SIZE);
    pion_free_dvector(EtimesL_input_2,1,INPUT_SIZE);
  }
  pion_free_dmatrix(ratePI,1,28,1,28);
  pion_free_dmatrix(rateRR,1,28,1,28);
  pion_free_dmatrix(rateDR,1,28,1,28);
  pion_free_dvector(LOWE_EGRID,1,LOWE_GRIDNUM);  
  pion_free_dvector(LOWE_PIGRID,1,LOWE_GRIDNUM); 
  pion_free_dvector(LOWE_PIGRID_2,1,LOWE_GRIDNUM);
  pion_free_dvector(specRR,1,SPECBINS);  
  pion_free_dvector(specDR,1,SPECBINS);  
  pion_free_dvector(specRR_temp,1,SPECBINS);  
  pion_free_dvector(specDR_temp,1,SPECBINS);  
  pion_free_dvector(spec,1,SPECBINS);  
  pion_free_dvector(L_kT,1,RECNUM);  
  pion_free_dvector(L_RR,1,RECNUM);  
  pion_free_dvector(L_RR_2,1,RECNUM);  
  pion_free_dvector(L_DR,1,RECNUM);  
  pion_free_dvector(L_DR_2,1,RECNUM);  
  pion_free_dvector(L_REC,1,RECNUM);  
  pion_free_dvector(L_REC_2,1,RECNUM);  
  pion_free_dvector(RR_line,1,RECNUM);  
  pion_free_dvector(RR_line_2,1,RECNUM);  
  pion_free_dvector(DR_line,1,RECNUM);  
  pion_free_dvector(DR_line_2,1,RECNUM);  
  pion_free_dvector(ABUND,1,30);
  pion_free_dvector(oshe,1,30);
  pion_free_dvector(E_array,1,SPECBINS);       
  pion_free_dvector(E_spectrum,1,SPECBINS);
  pion_free_dvector(type1_spectrum,1,SPECBINS);
  pion_free_dvector(type2_spectrum,1,SPECBINS);
  pion_free_dvector(type3_spectrum,1,SPECBINS);
  pion_free_dvector(type4_spectrum,1,SPECBINS);
  pion_free_dvector(type5_spectrum,1,SPECBINS);
  pion_free_dvector(type6_spectrum,1,SPECBINS);
  pion_free_dvector(type7_spectrum,1,SPECBINS);
  pion_free_dvector(type8_spectrum,1,SPECBINS);
  pion_free_dvector(l_array,1,SPECBINS);
  pion_free_dvector(l_spectrum,1,SPECBINS);
  pion_free_dvector(convert,1,SPECBINS);
  pion_free_dvector(abs_spectrum,1,SPECBINS);    
  pion_free_dvector(exc_spectrum,1,SPECBINS);    
  pion_free_dvector(rec_spectrum,1,SPECBINS);    
  pion_free_dvector(rec_spectrum_temp,1,SPECBINS);    
  pion_free_dvector(tau,1,SPECBINS);           
  pion_free_dvector(tau_exc,1,SPECBINS);           
  pion_free_dvector(tau_edge,1,SPECBINS);           
  pion_free_dvector(Labsorb,1,SPECBINS);    
  pion_free_dvector(int_array,1,SPECBINS);  
  pion_free_dvector(int_array_2,1,SPECBINS); 
  pion_free_dmatrix(Nion,1,28,1,28);
  pion_free_dmatrix(Tion,1,28,1,28);
  pion_free_dmatrix(EMion,1,28,1,28);
  pion_free_dmatrix(EM,1,28,1,28);
  pion_free_dmatrix(Hionizsigmaconv,1,28,1,SPECBINS);
  pion_free_dmatrix(Heionizsigmaconv,1,28,1,SPECBINS);
  pion_free_dvector(ionizsigmatemp,1,SPECBINS);
  pion_free_dvector(spectrumtemp,1,SPECBINS);
  pion_free_dmatrix(H_excite,1,28,1,6);
  pion_free_dmatrix(He_excite,1,28,1,9);
  pion_free_ivector(list,1,ELEMENTS);
  pion_free_dvector(Tvec,1,TEMPERATURES);
  pion_free_dvector(Yvec,1,TEMPERATURES);
  pion_free_dvector(Yvec2,1,TEMPERATURES);
  pion_free_dvector(EGRID,1,GRIDNUM);  
  pion_free_dvector(PIGRID,1,GRIDNUM); 
  pion_free_dvector(PIGRID_2,1,GRIDNUM);
  
  if (verbose) printf("Done!\n");
  
  /* A data file could not be opened: everything above ran on zeros. Report it
   * the same way as the other input errors, after the normal cleanup. */
  if (pion_ad_failed()) {
    for (i=0;i<ne;++i) photar[i]=0.;
    return 1;
  }
  return 0.;
}


/* integrand for photoionization rate or photoexcitation rate */


/* Maxwellian electron distribution */

/* Recombination-to-ground distribution (needs to be normalized) */

/* Recombination-to-ground distribution (needs to be normalized) */

/* For calculating recombination line strengths */

/* integrand for photoionization rate or photoexcitation rate */

/* Recombination cross-section in atomic units (a0^2 units) */


/* For convolution: Normalized gaussian profile centered at 0. with s=sigma */

/* Photoexcitation cross section [cm^2] from ground.
   global: E0, DELTANUD, OSCILLATOR, ALPHA */

/* Differential oscillator strength in atomic units: 1/Hartree = 1/(2.*13.6 eV) */

/* Photoionization cross-section in atomic units */


/* ground-state photoionization cross-section from Verner */

/* ground-state photoionization cross-section from Verner */


