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

int photoion
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(photoion,PHOTOION,photoion,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))
#define TEMPERATURES (11)

/* START */

/* Forward declarations (converted from K&R style). */


/* FINISH */


char* FGMSTR(char* name);

/* PHOTOION Subroutines */
/* Subroutines from PHOTOION */

/* MIABS Subroutines */
void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_edge_p[]) ;

double REC_spline_pi(double temp);


int photoion
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

  
  
  
  

  double *type1_spectrum,*type2_spectrum,*type3_spectrum,*type4_spectrum,*type5_spectrum,*type6_spectrum,*type7_spectrum,*type8_spectrum;

  double int_ans;


  double Elo,Ehi; /* Voigt function param.'s */
  double Ewidth,earBIN;

  int first=1;


  /* junk values for strings, ints, and floats */
  double djunk;


  double ELO=1.e-3,EHI=1.e6;  /* MAKE SURE THIS RANGE IS OK!!! CHECK HERE!!! */


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

  double redshift;
  double H0_cosmo=0.,Omega_m_cosmo=0.; /* from XSPEC, when D=0 */
  int fileincr;

  int LOG;


  double COLNORM,FLUXAVE;

  double *ABUND,*oshe;

  int klo,khi,jlow,jhigh;
  double spectemp;

  int DIST,lines=0;

  FILE *E_specfile=NULL,*l_specfile=NULL;
  FILE *E_output=NULL,*l_output=NULL;
  char *specfile_name;

  /* initialize photar array */
  for (i=0;i<ne;++i) photar[i]=0.;


  /* input.qdp is not read until ~900 lines below, by which point ~90 arrays are
   * live, around 40 of them SPECBINS-sized. The open there is checked, but it
   * returns without freeing any of them, so a missing file leaked tens of
   * megabytes per evaluation. Probe for it here, before anything is allocated.
   * The message and the return value are unchanged; the later check stays as a
   * guard. */
  if (param[6] > 0) {
    FILE *probe = fopen("input.qdp","r");
    if (probe == NULL) {
      printf("The file 'input.qdp' must exist in this directory.\n");
      return 0;
    }
    fclose(probe);
  }

  /* With D=0 and a non-zero redshift, the distance comes from the Hubble law,
   * using XSPEC's cosmology. Check it here, before anything is allocated. */
  if (param[14] == 0. && param[1] != 0.) {
    if (pion_cosmology(&H0_cosmo,&Omega_m_cosmo)) return 1;
  }

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
      printf("PHOTOION: Failed to open %s\n", probepath);
      return 1;
    }
    fclose(probe);
  }

  /* FILE NAMES */
  specfile_name=malloc(200);



  ABUND=pion_dvector(1,30);
  pion_ad_abundance(ABUND);

  oshe=pion_dvector(1,30);
  pion_ad_oscillator_he(oshe);

  /* Initialize the model array */
  for (i=0;i<ne;++i) photar[i] = 0.;

  Nion=pion_dmatrix(1,28,1,28);
  Tion=pion_dmatrix(1,28,1,28);
  EMion=pion_dmatrix(1,28,1,28);
  EM=pion_dmatrix(1,28,1,28);
  ratePI=pion_dmatrix(1,28,1,28);
  rateRR=pion_dmatrix(1,28,1,28);
  rateDR=pion_dmatrix(1,28,1,28);
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) {Nion[i][j]=0.;Tion[i][j]=10.;EMion[i][j]=0.;EM[i][j]=0.;ratePI[i][j]=0.;rateRR[i][j]=0.;rateDR[i][j]=0.;}

  type = (int) param[0];
  redshift = param[1];
  v_rad = 1.e5*param[2];   /* convert [km/s] to [cm/s] */
  v_trans = 1.e5*param[3];   /* convert [km/s] to [cm/s] */
  sigmav_rad = 1.e5*param[4];   /* convert [km/s] to [cm/s] */
    if (sigmav_rad>0.) lines=1; else if (sigmav_rad<=0.) lines=0;
  sigmav_trans = 1.e5*param[5]; /* convert [km/s] to [cm/s] */
  INPUT = param[6];
  INPUT_SHIFT = param[7];
  GAMMA = param[8];
  L_EMIN = 1000.*param[9]; /* from keV to eV */
  L_EMAX = 1000.*param[10]; /* from keV to eV */
  L_X = 1.e30*param[11];  /* convert [1e30 ergs/s] to [ergs/s] */
  FLUXAVE=param[12];
  f_COVERING = param[13];
  D = parsectocm*param[14];    /* convert [pc] to [cm] */
  EMIN=1000.*param[15]; /* from keV to eV */
  EMAX=1000.*param[16]; /* from keV to eV */
  SPECBINS=(int) param[17];
  fileincr = (int) param[18];
  verbose =(int) param[19];
  COLNORM=param[20];

  /* Total electron column density */
  N_e = param[21];

  /* H */
  Nion[1][1] = param[22];

  /* He */
  Nion[2][1] = param[23];
  Nion[2][2] = param[24];

  /* C */
  Nion[6][1] = param[25];
  Nion[6][2] = param[26];
  Nion[6][3] = param[27];
  Nion[6][4] = param[28];
  Nion[6][5] = param[29];
  Nion[6][6] = param[30];

  /* N */
  Nion[7][1] = param[31];
  Nion[7][2] = param[32];
  Nion[7][3] = param[33];
  Nion[7][4] = param[34];
  Nion[7][5] = param[35];
  Nion[7][6] = param[36];
  Nion[7][7] = param[37];

  /* O */
  Nion[8][1] = param[38];
  Nion[8][2] = param[39];
  Nion[8][3] = param[40];
  Nion[8][4] = param[41];
  Nion[8][5] = param[42];
  Nion[8][6] = param[43];
  Nion[8][7] = param[44];
  Nion[8][8] = param[45];

  /* Ne */
  Nion[10][1] = param[46];
  Nion[10][2] = param[47];
  Nion[10][3] = param[48];
  Nion[10][4] = param[49];
  Nion[10][5] = param[50];
  Nion[10][6] = param[51];
  Nion[10][7] = param[52];
  Nion[10][8] = param[53];
  Nion[10][9] = param[54];
  Nion[10][10] = param[55];

  /* Mg */
  Nion[12][1] = param[56];
  Nion[12][2] = param[57];
  Nion[12][3] = param[58];
  Nion[12][4] = param[59];
  Nion[12][5] = param[60];
  Nion[12][6] = param[61];
  Nion[12][7] = param[62];
  Nion[12][8] = param[63];
  Nion[12][9] = param[64];
  Nion[12][10] = param[65];
  Nion[12][11] = param[66];
  Nion[12][12] = param[67];

  /* Al */
  Nion[13][1] = param[68];
  Nion[13][2] = param[69];
  Nion[13][3] = param[70];
  Nion[13][4] = param[71];
  Nion[13][5] = param[72];
  Nion[13][6] = param[73];
  Nion[13][7] = param[74];
  Nion[13][8] = param[75];
  Nion[13][9] = param[76];
  Nion[13][10] = param[77];
  Nion[13][11] = param[78];
  Nion[13][12] = param[79];
  Nion[13][13] = param[80];

  /* Si */
  Nion[14][1] = param[81];
  Nion[14][2] = param[82];
  Nion[14][3] = param[83];
  Nion[14][4] = param[84];
  Nion[14][5] = param[85];
  Nion[14][6] = param[86];
  Nion[14][7] = param[87];
  Nion[14][8] = param[88];
  Nion[14][9] = param[89];
  Nion[14][10] = param[90];
  Nion[14][11] = param[91];
  Nion[14][12] = param[92];
  Nion[14][13] = param[93];
  Nion[14][14] = param[94];

  /* S */
  Nion[16][1] = param[95];
  Nion[16][2] = param[96];
  Nion[16][3] = param[97];
  Nion[16][4] = param[98];
  Nion[16][5] = param[99];
  Nion[16][6] = param[100];
  Nion[16][7] = param[101];
  Nion[16][8] = param[102];
  Nion[16][9] = param[103];
  Nion[16][10] = param[104];
  Nion[16][11] = param[105];
  Nion[16][12] = param[106];
  Nion[16][13] = param[107];
  Nion[16][14] = param[108];
  Nion[16][15] = param[109];
  Nion[16][16] = param[110];

  /* Ar */
  Nion[18][1] = param[111];
  Nion[18][2] = param[112];
  Nion[18][3] = param[113];
  Nion[18][4] = param[114];
  Nion[18][5] = param[115];
  Nion[18][6] = param[116];
  Nion[18][7] = param[117];
  Nion[18][8] = param[118];
  Nion[18][9] = param[119];
  Nion[18][10] = param[120];
  Nion[18][11] = param[121];
  Nion[18][12] = param[122];
  Nion[18][13] = param[123];
  Nion[18][14] = param[124];
  Nion[18][15] = param[125];
  Nion[18][16] = param[126];
  Nion[18][17] = param[127];
  Nion[18][18] = param[128];

  /* Ca */
  Nion[20][1] = param[129];
  Nion[20][2] = param[130];
  Nion[20][3] = param[131];
  Nion[20][4] = param[132];
  Nion[20][5] = param[133];
  Nion[20][6] = param[134];
  Nion[20][7] = param[135];
  Nion[20][8] = param[136];
  Nion[20][9] = param[137];
  Nion[20][10] = param[138];
  Nion[20][11] = param[139];
  Nion[20][12] = param[140];
  Nion[20][13] = param[141];
  Nion[20][14] = param[142];
  Nion[20][15] = param[143];
  Nion[20][16] = param[144];
  Nion[20][17] = param[145];
  Nion[20][18] = param[146];
  Nion[20][19] = param[147];
  Nion[20][20] = param[148];

  /* Fe */
  Nion[26][1] = param[149];
  Nion[26][2] = param[150];
  Nion[26][3] = param[151];
  Nion[26][4] = param[152];
  Nion[26][5] = param[153];
  Nion[26][6] = param[154];
  Nion[26][7] = param[155];
  Nion[26][8] = param[156];
  Nion[26][9] = param[157];
  Nion[26][10] = param[158];
  Nion[26][11] = param[159];
  Nion[26][12] = param[160];
  Nion[26][13] = param[161];
  Nion[26][14] = param[162];
  Nion[26][15] = param[163];
  Nion[26][16] = param[164];
  Nion[26][17] = param[165];
  Nion[26][18] = param[166];
  Nion[26][19] = param[167];
  Nion[26][20] = param[168];
  Nion[26][21] = param[169];
  Nion[26][22] = param[170];
  Nion[26][23] = param[171];
  Nion[26][24] = param[172];
  Nion[26][25] = param[173];
  Nion[26][26] = param[174];

  /* Ni */
  Nion[28][1] = param[175];
  Nion[28][2] = param[176];
  Nion[28][3] = param[177];
  Nion[28][4] = param[178];
  Nion[28][5] = param[179];
  Nion[28][6] = param[180];
  Nion[28][7] = param[181];
  Nion[28][8] = param[182];
  Nion[28][9] = param[183];
  Nion[28][10] = param[184];
  Nion[28][11] = param[185];
  Nion[28][12] = param[186];
  Nion[28][13] = param[187];
  Nion[28][14] = param[188];
  Nion[28][15] = param[189];
  Nion[28][16] = param[190];
  Nion[28][17] = param[191];
  Nion[28][18] = param[192];
  Nion[28][19] = param[193];
  Nion[28][20] = param[194];
  Nion[28][21] = param[195];
  Nion[28][22] = param[196];
  Nion[28][23] = param[197];
  Nion[28][24] = param[198];
  Nion[28][25] = param[199];
  Nion[28][26] = param[200];
  Nion[28][27] = param[201];
  Nion[28][28] = param[202];
  Tion[6][1] = param[203];
  Tion[6][2] = param[204];
  Tion[7][1] = param[205];
  Tion[7][2] = param[206];
  Tion[8][1] = param[207];
  Tion[8][2] = param[208];
  Tion[10][1] = param[209];
  Tion[10][2] = param[210];
  Tion[10][3] = param[211];
  Tion[10][4] = param[212];
  Tion[10][5] = param[213];
  Tion[10][6] = param[214];
  Tion[10][7] = param[215];
  Tion[10][8] = param[216];
  Tion[10][9] = param[217];
  Tion[10][10] = param[218];
  Tion[12][1] = param[219];
  Tion[12][2] = param[220];
  Tion[12][3] = param[221];
  Tion[12][4] = param[222];
  Tion[12][5] = param[223];
  Tion[12][6] = param[224];
  Tion[12][7] = param[225];
  Tion[12][8] = param[226];
  Tion[12][9] = param[227];
  Tion[12][10] = param[228];
  Tion[13][1] = param[229];
  Tion[13][2] = param[230];
  Tion[13][3] = param[231];
  Tion[13][4] = param[232];
  Tion[13][5] = param[233];
  Tion[13][6] = param[234];
  Tion[13][7] = param[235];
  Tion[13][8] = param[236];
  Tion[13][9] = param[237];
  Tion[13][10] = param[238];
  Tion[14][1] = param[239];
  Tion[14][2] = param[240];
  Tion[14][3] = param[241];
  Tion[14][4] = param[242];
  Tion[14][5] = param[243];
  Tion[14][6] = param[244];
  Tion[14][7] = param[245];
  Tion[14][8] = param[246];
  Tion[14][9] = param[247];
  Tion[14][10] = param[248];
  Tion[16][1] = param[249];
  Tion[16][2] = param[250];
  Tion[16][3] = param[251];
  Tion[16][4] = param[252];
  Tion[16][5] = param[253];
  Tion[16][6] = param[254];
  Tion[16][7] = param[255];
  Tion[16][8] = param[256];
  Tion[16][9] = param[257];
  Tion[16][10] = param[258];
  Tion[18][1] = param[259];
  Tion[18][2] = param[260];
  Tion[18][3] = param[261];
  Tion[18][4] = param[262];
  Tion[18][5] = param[263];
  Tion[18][6] = param[264];
  Tion[18][7] = param[265];
  Tion[18][8] = param[266];
  Tion[18][9] = param[267];
  Tion[18][10] = param[268];
  Tion[20][1] = param[269];
  Tion[20][2] = param[270];
  Tion[20][3] = param[271];
  Tion[20][4] = param[272];
  Tion[20][5] = param[273];
  Tion[20][6] = param[274];
  Tion[20][7] = param[275];
  Tion[20][8] = param[276];
  Tion[20][9] = param[277];
  Tion[20][10] = param[278];
  Tion[26][1] = param[279];
  Tion[26][2] = param[280];
  Tion[26][3] = param[281];
  Tion[26][4] = param[282];
  Tion[26][5] = param[283];
  Tion[26][6] = param[284];
  Tion[26][7] = param[285];
  Tion[26][8] = param[286];
  Tion[26][9] = param[287];
  Tion[26][10] = param[288];
  Tion[28][1] = param[289];
  Tion[28][2] = param[290];
  Tion[28][3] = param[291];
  Tion[28][4] = param[292];
  Tion[28][5] = param[293];
  Tion[28][6] = param[294];
  Tion[28][7] = param[295];
  Tion[28][8] = param[296];
  Tion[28][9] = param[297];
  Tion[28][10] = param[298];

  if (type == 4 || type ==5) {v_trans=v_rad; sigmav_trans=sigmav_rad;}

  LOG = 0;
  if (LOG!=0) {
    for (i=1;i<=28;++i) {
      for (j=1;j<=28;++j) { 
	if (Nion[i][j]) {
	  Nion[i][j]=pow(10.,Nion[i][j]);
	}
      }
    }
    if (N_e) N_e=pow(10.,N_e);
  }

  if (COLNORM!=1.) {
    for (i=1;i<=28;++i) for (j=1;j<=28;++j) Nion[i][j]*=COLNORM;
    N_e*=COLNORM;
  }

  /* SAME AS "PHSI" and "PHXI" FROM HERE ON OUT */

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
    if (pion_input_continuum(redshift,type,verbose,Labsorb)) return 0;   /* no input.qdp: as before; the entry probe catches it first */
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


