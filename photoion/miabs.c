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

int miabs
(float *ear,int ne,float *param,int ifl,float *photar,float *photer);

FCALLSCSUB6(miabs,MIABS,miabs,FLOATV,INT,FLOATV,INT,FLOATV,FLOATV)

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))


/* START */

/* Forward declarations (converted from K&R style). */


/* FINISH */


char* FGMSTR(char* name);

/* Subroutines from PHOTOION */

/* MIABS Subroutines */
void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_edge_p[]) ;


int miabs
(float *ear,int ne,float *param,int ifl,float *photar,float *photer)
{
  /*ear[0] -> ear[ne]   in keV
    photar[0] -> photar[ne-1]  
    photar[i] = [photons/cm^2/s] for bin ear[i] -> ear[i+1]
    param[0] -> param[TOTAL-1]   gives XSPEC parameters from 1 to TOTAL*/

  const struct VERNER_STRUCT (*vernerionizsigma)[31];  /* cached, see photoion_atomdata.h */

  /*  struct VERNER_PARTIAL_STRUCT partialsigma[29][PARTIAL_NUM];*/
  const struct VERNER_PARTIAL_STRUCT (*partialsigma)[125];
  const int *npartial;
  int pathlen;     /* buffer size for any data path under DATADIR (open issue 2) */

  struct HYDROGEN_STRUCT {
    double lambda[8];
    double f[8];
    double A[8];
    double b[8]; /* branching ratios? */
  };
  struct HYDROGEN_STRUCT hydrogen[29];
  
  struct HELIUM_STRUCT {
    double lambda[11];
    double f[11];
    double A[11];
    double b[11]; /* branching ratios? */
  };
  struct HELIUM_STRUCT helium[29];
  
  struct HIGHER_ORDER_STRUCT {
    double lambda[101];
    double f[101];
  };
  struct HIGHER_ORDER_STRUCT highn[27][3];

  double redshift,v_rad,sigmav_rad,N_e,**Nion,COLNORM;
  int verbose;

  int LOG;

  double Elo,Ehi; /* Voigt function param.'s */
  double earBIN,Ewidth;

  /* junk values for strings, ints, and floats */
  char *sjunk,*line,*sjunk1,*sjunk2,*sjunk3,*sjunk4;
  double djunk;

  double g_i,g_j;

  double ELO=1.e-3,EHI=1.e6;  /* MAKE SURE THIS RANGE IS OK!!! CHECK HERE!!! */

  double en,Atemp,ftemp,ga;

  char *root,*element_name,*temp;
  const char *ext;   /* always a string literal, never owned */
  int i,j,k,n;
  /* *.rr file */
  double p0,p1,p2,p3;

  double grounddeg[3],leveldeg[3][11],freedeg[3];

  double earlo,earhi;
  
  int HIGHN=50;

  int ELEMENTS=12;

  double E0=0.,DELTANUD,ALPHA,OSCILLATOR;

  int element,electron;
  int LINE;
  double WAVE,AMtemp,fMtemp;
  double Etemp;

  int *list;

  char name[]="PHOTOION_DIR";
  char *DATADIR;

  double **H_excite, **He_excite;

  double *ABUND,*oshe,oscillatornorm=0;

  int klo,khi,jlow,jhigh;

  int DIST,lines=0;

  FILE *linedat,*highnfile;
  char *linedat_name,*highnfile_name;
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
      printf("MIABS: Failed to open %s\n", probepath);
      return 1;
    }
    fclose(probe);
  }

  /* FILE NAMES */
  linedat_name=malloc(pathlen);
  highnfile_name=malloc(pathlen);

  root=malloc(pathlen);
  element_name=malloc(3);  /* 2-char symbols ("Ne") need 3 bytes with the NUL */
  temp=malloc(pathlen);
  sjunk=malloc(50);
  sjunk1=malloc(50);
  sjunk2=malloc(50);
  sjunk3=malloc(50);
  sjunk4=malloc(50);
  line=malloc(400);


  ABUND=pion_dvector(1,30);
  pion_ad_abundance(ABUND);

  oshe=pion_dvector(1,30);
  pion_ad_oscillator_he(oshe);

  Nion=pion_dmatrix(1,28,1,28);
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) Nion[i][j]=0.;

  redshift = param[0];     /* cosmological redshift */
  v_rad = 1.e5*param[1];   /* velocity shift: convert [km/s] to [cm/s] */
  sigmav_rad = 1.e5*param[2];   /* gaussian (1-sigma) velocity width: convert [km/s] to [cm/s] */
  if (sigmav_rad>0.) lines=1; else if (sigmav_rad<=0.) lines=0;
  EMIN=1000.*param[3]; /* from keV to eV */
  EMAX=1000.*param[4]; /* from keV to eV */
  SPECBINS=(int) param[5];
  verbose =(int) param[6];
  LOG=(int) param[7];  /* LOG=0 (normal units for N_ion), LOG=1 (log units)*/
  COLNORM=param[8];

  /* Total electron column density [cm^-2] */
  N_e = param[9];

  /* H */
  Nion[1][1] = param[10];

  /* He */
  Nion[2][1] = param[11];
  Nion[2][2] = param[12];

  /* C */
  Nion[6][1] = param[13];
  Nion[6][2] = param[14];
  Nion[6][3] = param[15];
  Nion[6][4] = param[16];
  Nion[6][5] = param[17];
  Nion[6][6] = param[18];

  /* N */
  Nion[7][1] = param[19];
  Nion[7][2] = param[20];
  Nion[7][3] = param[21];
  Nion[7][4] = param[22];
  Nion[7][5] = param[23];
  Nion[7][6] = param[24];
  Nion[7][7] = param[25];

  /* O */
  Nion[8][1] = param[26];
  Nion[8][2] = param[27];
  Nion[8][3] = param[28];
  Nion[8][4] = param[29];
  Nion[8][5] = param[30];
  Nion[8][6] = param[31];
  Nion[8][7] = param[32];
  Nion[8][8] = param[33];

  /* Ne */
  Nion[10][1] = param[34];
  Nion[10][2] = param[35];
  Nion[10][3] = param[36];
  Nion[10][4] = param[37];
  Nion[10][5] = param[38];
  Nion[10][6] = param[39];
  Nion[10][7] = param[40];
  Nion[10][8] = param[41];
  Nion[10][9] = param[42];
  Nion[10][10] = param[43];

  /* Mg */
  Nion[12][1] = param[44];
  Nion[12][2] = param[45];
  Nion[12][3] = param[46];
  Nion[12][4] = param[47];
  Nion[12][5] = param[48];
  Nion[12][6] = param[49];
  Nion[12][7] = param[50];
  Nion[12][8] = param[51];
  Nion[12][9] = param[52];
  Nion[12][10] = param[53];
  Nion[12][11] = param[54];
  Nion[12][12] = param[55];

  /* Al */
  Nion[13][1] = param[56];
  Nion[13][2] = param[57];
  Nion[13][3] = param[58];
  Nion[13][4] = param[59];
  Nion[13][5] = param[60];
  Nion[13][6] = param[61];
  Nion[13][7] = param[62];
  Nion[13][8] = param[63];
  Nion[13][9] = param[64];
  Nion[13][10] = param[65];
  Nion[13][11] = param[66];
  Nion[13][12] = param[67];
  Nion[13][13] = param[68];

  /* Si */
  Nion[14][1] = param[69];
  Nion[14][2] = param[70];
  Nion[14][3] = param[71];
  Nion[14][4] = param[72];
  Nion[14][5] = param[73];
  Nion[14][6] = param[74];
  Nion[14][7] = param[75];
  Nion[14][8] = param[76];
  Nion[14][9] = param[77];
  Nion[14][10] = param[78];
  Nion[14][11] = param[79];
  Nion[14][12] = param[80];
  Nion[14][13] = param[81];
  Nion[14][14] = param[82];

  /* S */
  Nion[16][1] = param[83];
  Nion[16][2] = param[84];
  Nion[16][3] = param[85];
  Nion[16][4] = param[86];
  Nion[16][5] = param[87];
  Nion[16][6] = param[88];
  Nion[16][7] = param[89];
  Nion[16][8] = param[90];
  Nion[16][9] = param[91];
  Nion[16][10] = param[92];
  Nion[16][11] = param[93];
  Nion[16][12] = param[94];
  Nion[16][13] = param[95];
  Nion[16][14] = param[96];
  Nion[16][15] = param[97];
  Nion[16][16] = param[98];

  /* Ar */
  Nion[18][1] = param[99];
  Nion[18][2] = param[100];
  Nion[18][3] = param[101];
  Nion[18][4] = param[102];
  Nion[18][5] = param[103];
  Nion[18][6] = param[104];
  Nion[18][7] = param[105];
  Nion[18][8] = param[106];
  Nion[18][9] = param[107];
  Nion[18][10] = param[108];
  Nion[18][11] = param[109];
  Nion[18][12] = param[110];
  Nion[18][13] = param[111];
  Nion[18][14] = param[112];
  Nion[18][15] = param[113];
  Nion[18][16] = param[114];
  Nion[18][17] = param[115];
  Nion[18][18] = param[116];

  /* Ca */
  Nion[20][1] = param[117];
  Nion[20][2] = param[118];
  Nion[20][3] = param[119];
  Nion[20][4] = param[120];
  Nion[20][5] = param[121];
  Nion[20][6] = param[122];
  Nion[20][7] = param[123];
  Nion[20][8] = param[124];
  Nion[20][9] = param[125];
  Nion[20][10] = param[126];
  Nion[20][11] = param[127];
  Nion[20][12] = param[128];
  Nion[20][13] = param[129];
  Nion[20][14] = param[130];
  Nion[20][15] = param[131];
  Nion[20][16] = param[132];
  Nion[20][17] = param[133];
  Nion[20][18] = param[134];
  Nion[20][19] = param[135];
  Nion[20][20] = param[136];

  /* Fe */
  Nion[26][1] = param[137];
  Nion[26][2] = param[138];
  Nion[26][3] = param[139];
  Nion[26][4] = param[140];
  Nion[26][5] = param[141];
  Nion[26][6] = param[142];
  Nion[26][7] = param[143];
  Nion[26][8] = param[144];
  Nion[26][9] = param[145];
  Nion[26][10] = param[146];
  Nion[26][11] = param[147];
  Nion[26][12] = param[148];
  Nion[26][13] = param[149];
  Nion[26][14] = param[150];
  Nion[26][15] = param[151];
  Nion[26][16] = param[152];
  Nion[26][17] = param[153];
  Nion[26][18] = param[154];
  Nion[26][19] = param[155];
  Nion[26][20] = param[156];
  Nion[26][21] = param[157];
  Nion[26][22] = param[158];
  Nion[26][23] = param[159];
  Nion[26][24] = param[160];
  Nion[26][25] = param[161];
  Nion[26][26] = param[162];

  /* Ni */
  Nion[28][1] = param[163];
  Nion[28][2] = param[164];
  Nion[28][3] = param[165];
  Nion[28][4] = param[166];
  Nion[28][5] = param[167];
  Nion[28][6] = param[168];
  Nion[28][7] = param[169];
  Nion[28][8] = param[170];
  Nion[28][9] = param[171];
  Nion[28][10] = param[172];
  Nion[28][11] = param[173];
  Nion[28][12] = param[174];
  Nion[28][13] = param[175];
  Nion[28][14] = param[176];
  Nion[28][15] = param[177];
  Nion[28][16] = param[178];
  Nion[28][17] = param[179];
  Nion[28][18] = param[180];
  Nion[28][19] = param[181];
  Nion[28][20] = param[182];
  Nion[28][21] = param[183];
  Nion[28][22] = param[184];
  Nion[28][23] = param[185];
  Nion[28][24] = param[186];
  Nion[28][25] = param[187];
  Nion[28][26] = param[188];
  Nion[28][27] = param[189];
  Nion[28][28] = param[190];

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
  if (verbose) printf("4*sigmav_rad = %d bins\n",DIST);

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
  pion_verner_edges(Nion,vernerionizsigma,partialsigma,npartial,28,tau_edge);

  if (lines) {
    if (verbose) printf("Determining LOW-n Photoexcitation Cross Sections & Opacity for C to Fe...\n");
    grounddeg[1]=2.;
    freedeg[1]=1.;
    for (i=1;i<=6;++i) leveldeg[1][i]=6.;
    grounddeg[2]=1.;
    freedeg[2]=2.;
    for (i=3;i<=8;++i) leveldeg[2][i]=3.;
    snprintf(linedat_name,pathlen,"%s/photoion_dat/line.dat",DATADIR);
    linedat=fopen(linedat_name,"r");
    while(fscanf(linedat,"%s%d%d%d%s%lf%s%lf%lf",sjunk,&element,&electron,&LINE,sjunk,&WAVE,sjunk,&AMtemp,&fMtemp) != EOF) {
      if (electron == 1) {
	hydrogen[element].lambda[LINE]=WAVE;
	hydrogen[element].A[LINE]=AMtemp;
	hydrogen[element].f[LINE]=fMtemp;
	if (Nion[element][electron] && fMtemp && LINE <= 4) {
	  E0=1000.*HC_KEV_ANGSTROM/hydrogen[element].lambda[LINE];
	  E0=E0*doppler_rad;
	  OSCILLATOR=hydrogen[element].f[LINE];
	  g_j=leveldeg[electron][LINE];
	  g_i=grounddeg[electron];
	  Atemp=hydrogen[element].A[LINE];
	  ga=Atemp/*g_j*Atemp/(3.*g_i*OSCILLATOR)*/;
	  DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	  ALPHA=ga/(4.*PI*DELTANUD);
	  pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc,1.0,1.0,0);
	}
      } else if (electron == 2) {
	helium[element].lambda[LINE]=WAVE;
	helium[element].A[LINE]=AMtemp;
	helium[element].f[LINE]=fMtemp;
	if (Nion[element][electron] && fMtemp  && LINE <= 6) {
	  E0=1000.*HC_KEV_ANGSTROM/helium[element].lambda[LINE];
	  E0=E0*doppler_rad;
	  OSCILLATOR=helium[element].f[LINE];
	  g_j=leveldeg[electron][LINE];
	  g_i=grounddeg[electron];
	  Atemp=helium[element].A[LINE];
	  ga=Atemp/*g_j*Atemp/(3.*g_i*OSCILLATOR)*/;
	  DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	  ALPHA=ga/(4.*PI*DELTANUD);
	  pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc,1.0,1.0,0);
	}
      }
    }
    fclose(linedat);

    if (verbose) printf("Determining HIGH-n Photoexcitation Cross Sections & Opacity for C to Fe\n");
    snprintf(highnfile_name,pathlen,"%s/photoion_dat/highn.dat",DATADIR);
    highnfile=fopen(highnfile_name,"r");
    while (fscanf(highnfile,"%d%d%d%lf%lf",&element,&electron,&n,&Etemp,&ftemp)!=EOF) {
      highn[element][electron].lambda[n]=HC_KEV_ANGSTROM/Etemp*1000.;
      highn[element][electron].f[n]=ftemp;
    }
    fclose(highnfile);
    for (i=1;i<=ELEMENTS;++i) {
      element=list[i];
      for (electron=1;electron<=2;++electron) {
	if (Nion[element][electron]) {
	  if (electron==1) oscillatornorm=1.6; /* Bethe-Salpeter p. 265 */
	  if (electron==2) oscillatornorm=oshe[element]; /* defined above */
	  for (n=6;n<=HIGHN;++n) {
	    OSCILLATOR=oscillatornorm/cube((double) n);
	    if (element!=28) {
	      E0=HC_KEV_ANGSTROM/highn[element][electron].lambda[n]*1000.;
	    } else if (electron==1) {/* Use Fe numbers for Ni */
	      E0=HC_KEV_ANGSTROM/(highn[26][electron].lambda[n]/1.1614)*1000.;
	    } else if (electron==2) {/* Use Fe numbers for Ni */
	      E0=HC_KEV_ANGSTROM/(highn[26][electron].lambda[n]/1.165)*1000.;
	    }
	    E0=E0*doppler_rad;
	    DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	    /* taking classical value for "ga=gamma" from p. 112-114 B+D */
	    ga=DAMPING_CLASSICAL*sqr(E0*eVtoergs/hhh);
	    ALPHA=ga/(4.*PI*DELTANUD);
	    pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc,1.0,1.0,0);
	  }
	}
      }
    }
  }    


  /* File root for L-shell input files for Ne through Ca */
  if (verbose) printf("L-shell ions:  Ne through Ni\n");
  for (element=10;element<=28;++element) {
    for (electron=1;electron<=10;++electron) {
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
	/* for photoionization cross-sections */
	ext="pi_short"; 
	if (electron<10) snprintf(temp,pathlen,"%s0%da.%s",root,electron,ext);
	else snprintf(temp,pathlen,"%s%da.%s",root,electron,ext);
	input=fopen(temp,"r");
	while (fscanf(input,"%d%lf%d%lf%lf%lf",&i,&g_i,&j,&g_j,&THRESHOLD,&ANGULAR) != EOF) {
	  i=i+1; j=j+1; /* lowest level is '1' not '0'!!! */
	  g_i=g_i+1.;
	  g_j=g_j+1.;
	  fscanf(input,"%lf%lf%lf%lf",&p0,&p1,&p2,&p3);
	  for (k=1;k<=LOWE_GRIDNUM;++k) {
	    fscanf(input,"%lf%lf%lf%lf",&(LOWE_EGRID[k]),&djunk,&(LOWE_PIGRID[k]),&djunk);
	    LOWE_EGRID[k]=log10(LOWE_EGRID[k]+THRESHOLD);
	    LOWE_PIGRID[k]=log10(1.e-20*LOWE_PIGRID[k]);
	  }
	  pion_spline(LOWE_EGRID,LOWE_PIGRID,LOWE_GRIDNUM,1.e40,1.e40,LOWE_PIGRID_2);	  
	  if (Nion[element][electron]*pion_pisigma(g_i,p0,p1,p2,p3,THRESHOLD) >= tau_lim) {
	    /* PHOTOIONIZATION OUT OF GRD STATE OF n+1 ION */
	    for (k=1;k<=GRIDNUM;++k) {
	      if (EGRID[k]<log10(THRESHOLD)) PIGRID[k]=1.e-90;
	      else if (EGRID[k]>=log10(THRESHOLD) && EGRID[k]<=LOWE_EGRID[LOWE_GRIDNUM]) PIGRID[k]=pion_lowEpispline(pow(10.,EGRID[k]));
	      else PIGRID[k]=pion_pisigma(g_i,p0,p1,p2,p3,pow(10.,EGRID[k]));
	    }
	    for (k=1;k<=GRIDNUM;++k) PIGRID[k]=log10(PIGRID[k]);
	    pion_spline(EGRID,PIGRID,GRIDNUM,1.e40,1.e40,PIGRID_2);
	    
	    /* calculate opacity of given edge and modify "tau" accordingly */
	    pion_fac_edge_opacity(Nion[element][electron],THRESHOLD,tau_edge);
	  }
	}
	fclose(input);
	
	if (lines) {
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
	/* for photoionization cross-sections */
	ext="pi_short"; 
	if (electron<10) snprintf(temp,pathlen,"%s0%da.%s",root,electron,ext);
	else snprintf(temp,pathlen,"%s%da.%s",root,electron,ext);
	input=fopen(temp,"r");
	while (fscanf(input,"%d%lf%d%lf%lf%lf",&i,&g_i,&j,&g_j,&THRESHOLD,&ANGULAR) != EOF) {
	  i=i+1; j=j+1; /* lowest level is '1' not '0'!!! */
	  g_i=g_i+1.;
	  g_j=g_j+1.;
	  fscanf(input,"%lf%lf%lf%lf",&p0,&p1,&p2,&p3);
	  for (k=1;k<=LOWE_GRIDNUM;++k) {
	    fscanf(input,"%lf%lf%lf%lf",&(LOWE_EGRID[k]),&djunk,&(LOWE_PIGRID[k]),&djunk);
	    LOWE_EGRID[k]=log10(LOWE_EGRID[k]+THRESHOLD);
	    LOWE_PIGRID[k]=log10(1.e-20*LOWE_PIGRID[k]);
	  }
	  pion_spline(LOWE_EGRID,LOWE_PIGRID,LOWE_GRIDNUM,1.e40,1.e40,LOWE_PIGRID_2);	  
	  if (Nion[element][electron]*pion_pisigma(g_i,p0,p1,p2,p3,THRESHOLD) >= tau_lim) {
	    /* PHOTOIONIZATION OUT OF GRD STATE OF n+1 ION */
	    for (k=1;k<=GRIDNUM;++k) {
	      if (EGRID[k]<log10(THRESHOLD)) PIGRID[k]=1.e-90;
	      else if (EGRID[k]>=log10(THRESHOLD) && EGRID[k]<=LOWE_EGRID[LOWE_GRIDNUM]) PIGRID[k]=pion_lowEpispline(pow(10.,EGRID[k]));
	      else PIGRID[k]=pion_pisigma(g_i,p0,p1,p2,p3,pow(10.,EGRID[k]));
	    }
	    for (k=1;k<=GRIDNUM;++k) PIGRID[k]=log10(PIGRID[k]);
	    pion_spline(EGRID,PIGRID,GRIDNUM,1.e40,1.e40,PIGRID_2);
	    
	    /* calculate opacity of given edge and modify "tau" accordingly */
	    pion_fac_edge_opacity(Nion[element][electron],THRESHOLD,tau_edge);
	  }
	}
	fclose(input);
	
	if (lines) {
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

	if (lines) {
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

	/* Edges */
	ext="pi_short"; 
	snprintf(temp,pathlen,"%s%2da.%s",root,electron,ext);
	input=fopen(temp,"r");
	while (fscanf(input,"%d%lf%d%lf%lf%lf",&i,&g_i,&j,&g_j,&THRESHOLD,&ANGULAR) != EOF) {
	  i=i+1; j=j+1; /* lowest level is '1' not '0'!!! */
	  g_i=g_i+1.;
	  g_j=g_j+1.;
	  fscanf(input,"%lf%lf%lf%lf",&p0,&p1,&p2,&p3);
	  for (k=1;k<=LOWE_GRIDNUM;++k) {
	    fscanf(input,"%lf%lf%lf%lf",&(LOWE_EGRID[k]),&djunk,&(LOWE_PIGRID[k]),&djunk);
	    LOWE_EGRID[k]=log10(LOWE_EGRID[k]+THRESHOLD);
	    LOWE_PIGRID[k]=log10(1.e-20*LOWE_PIGRID[k]);
	  }
	  pion_spline(LOWE_EGRID,LOWE_PIGRID,LOWE_GRIDNUM,1.e40,1.e40,LOWE_PIGRID_2);	  
	  if (Nion[element][electron]*pion_pisigma(g_i,p0,p1,p2,p3,1.01*THRESHOLD) >= tau_lim) {
	    /* PHOTOIONIZATION OUT OF GRD STATE OF n+1 ION */
	    for (k=1;k<=GRIDNUM;++k) {
	      if (EGRID[k]<log10(THRESHOLD)) PIGRID[k]=1.e-90;
	      else if (EGRID[k]>=log10(THRESHOLD) && EGRID[k]<=LOWE_EGRID[LOWE_GRIDNUM]) PIGRID[k]=pion_lowEpispline(pow(10.,EGRID[k]));
	      else PIGRID[k]=pion_pisigma(g_i,p0,p1,p2,p3,pow(10.,EGRID[k]));
	    }
	    for (k=1;k<=GRIDNUM;++k) PIGRID[k]=log10(PIGRID[k]);
	    pion_spline(EGRID,PIGRID,GRIDNUM,1.e40,1.e40,PIGRID_2);
	    
	    /* calculate opacity of given edge and modify "tau" accordingly */
	    pion_fac_edge_opacity(Nion[element][electron],THRESHOLD,tau_edge);
	  }
	}
	fclose(input);
      }
    }
  }
  

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
  /* Both were freed inside the `if (lines)` block that used them, while the
   * malloc at the top is unconditional, so each leaked whenever that branch was
   * skipped -- `lines` is 0 whenever sigma_v is 0, which this model's default
   * of 100 km/s avoids, but a user setting it to 0 does not. Freed here instead. */
  free(linedat_name);
  free(highnfile_name);
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
  pion_free_dmatrix(Nion,1,28,1,28);
  pion_free_dvector(ionizsigmatemp,1,SPECBINS);
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
  return 0.;
}


/* For convolution: Normalized gaussian profile centered at 0. with s=sigma */

/* Photoexcitation cross section [cm^2] from ground.
   global: E0, DELTANUD, OSCILLATOR, ALPHA */

/* Differential oscillator strength in atomic units: 1/Hartree = 1/(2.*13.6 eV) */

/* Photoionization cross-section in atomic units */


/* ground-state photoionization cross-section from Verner */

/* ground-state photoionization cross-section from Verner */


