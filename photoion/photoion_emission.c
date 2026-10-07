/* The reemission half of the emission models. See the header.
 *
 * This is the type>1 block of photoion.c moved verbatim, with one kind of
 * change: the data it read on every evaluation (trates, tr_shorter, pi_short,
 * .dat, rr_short, H_ and He_recombination) now comes from photoion_atomdata,
 * parsed once. The arithmetic, and the order in which it accumulates into the
 * spectra, are unchanged. */
#include <stdio.h>
#include <math.h>

#include "photoion_alloc.h"
#include "photoion_spline.h"
#include "photoion_state.h"
#include "photoion_phys.h"
#include "photoion_const.h"
#include "photoion_atomdata.h"
#include "photoion_integrate.h"
#include "photoion_emission.h"

#define sqr(X) ((X)*(X))
#define cube(X) ((X)*(X)*(X))

/* The H/He-like rates and recombination spectra need the ion's Verner
 * photoionization fit. Without one (Emax 0) they would be 0/0 = NaN in
 * every bin -- what Ni gave before verner_photo.dat had Ni rows. Skip such
 * an ion instead, and say so once. */
static int have_verner(const struct VERNER_STRUCT *v, int element, int electron)
{
  static int warned[31][3];
  if (v->Emax > 0.) return 1;
  if (element >= 0 && element <= 30 && electron >= 1 && electron <= 2 && !warned[element][electron]) {
    printf("PHOTOION: no Verner photoionization fit for Z=%d with %d electron%s; its recombination emission is left out\n",
           element, electron, electron == 1 ? "" : "s");
    warned[element][electron] = 1;
  }
  return 0;
}

void pion_emission_reemission(double **Nion, double N_e, double sigmav_rad,
                              int lines, int verbose,
                              const int list[], int ELEMENTS, int HIGHN,
                              const double oshe[], const double ABUND[], const double Labsorb[],
                              double **H_excite, double **He_excite,
                              double **ratePI, double **rateRR, double **rateDR,
                              double specRR[], double specDR[])
{
  int i,j,k,n,r,g,element,electron,LINE,typenum,itemp,jtemp,SUMlo,SUMhi;
  int ntrsh,npirec,nrrl,nrrs;
  double E0,OSCILLATOR,ga,DELTANUD,ALPHA,Elo,Ehi,int_junk,int_ans,strength,oscillatornorm;
  double ratePE=0.;   /* read only if a photoexcitation integral comes out negative, which a sum of non-negative terms cannot */
  double Ttemp,RRtemp,DRtemp,g_i,g_j,ftemp,Atemp,en,p0,p1,p2,p3,djunk,RECNORM,R;
  double intMIN,intMAX;
  const struct PION_TRATES_ROW *trates;
  const struct PION_TRSHORTER_ROW *trsh;
  const struct PION_FAC_PIREC *pirec, *rrs;
  const struct PION_RRLINE_ROW *rrl;
  const struct H_REC_STRUCT *H_rec;
  const struct HE_REC_STRUCT *He_rec;
  const struct HYDROGEN_STRUCT *hydrogen=pion_ad_hydrogen();
  const struct HELIUM_STRUCT *helium=pion_ad_helium();
  const struct HIGHER_ORDER_STRUCT (*highn)[3]=pion_ad_highn();
  const struct VERNER_STRUCT (*vernerionizsigma)[31]=pion_ad_verner_full();

    /**************************************************************************/
    /*                  Determine Emission Line spectrum                      */
    /**************************************************************************/
    if (lines) {
      if (verbose) printf("Determining LOW-n Photoexcitation Rates (n<=5) \n");
      /* hydrogenic */      
      for (i=1;i<=ELEMENTS;++i) {
	electron=1;
	element=list[i];
	if (Nion[element][electron]) {
	  for (LINE=1;LINE<=4;++LINE) {
	    OSCILLATOR=hydrogen[element].f[LINE];
	    if (OSCILLATOR) {
	      E0=1000.*HC_KEV_ANGSTROM/hydrogen[element].lambda[LINE];
	      E0=E0*doppler_rad;
	      ga=hydrogen[element].A[LINE]/*leveldeg[electron][LINE]/(3.*grounddeg[electron]*hydrogen[element].f[LINE])*/;
	      DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	      ALPHA=ga/(4.*PI*DELTANUD);
	      pion_line_limits(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,&Elo,&Ehi,&SUMlo,&SUMhi,0.9,1.1,1);
	      int_junk=0.;
	      for (k=SUMlo;k<=SUMhi;++k) {
		int_array[k]=pion_excitsigma(E0,OSCILLATOR,DELTANUD,ALPHA,E_array[k])*Labsorb[k];
		int_junk+=EBIN*int_array[k];
	      }
	      /*	    pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
			    int_ans=pion_integrate(pion_integrand,1.001*Elo,0.999*Ehi);*/
	      if (1 || int_ans<0.) int_ans=int_junk;
	      H_excite[element][LINE]=f_COVERING*Nion[element][electron]*int_ans;
	      if (H_excite[element][LINE]<0.) H_excite[element][LINE]=f_COVERING*Nion[element][electron]*ratePE;
	      /* add line to Seyfert 2 spectrum */
	      E0=doppler_trans/doppler_rad*E0;
	      k=(int) ((E0-EMIN+0.5*EBIN)/EBIN);
	      if ((E0-EMIN+0.5*EBIN)/EBIN-(double) k >= 0.5) ++k;
	      if (k <= SPECBINS && k >= 1) {
		exc_spectrum[k]+=H_excite[element][LINE]/EBIN;
	      }
	    }
	  }
	}
      }
      /* He-like photoexcitation - low n */
      for (i=1;i<=ELEMENTS;++i) {
	element=list[i];
	electron=2;
	if (Nion[element][electron]) {
	  for (LINE=3;LINE<=6;++LINE) {
	    OSCILLATOR=helium[element].f[LINE];
	    if (OSCILLATOR) {
	      E0=1000.*HC_KEV_ANGSTROM/helium[element].lambda[LINE];
	      E0=E0*doppler_rad;
	      ga=helium[element].A[LINE]/*leveldeg[electron][LINE]/(3.*grounddeg[electron]*helium[element].f[LINE])*/;
	      DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	      ALPHA=ga/(4.*PI*DELTANUD);
	      pion_line_limits(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,&Elo,&Ehi,&SUMlo,&SUMhi,0.9,1.1,1);
	      int_junk=0.;
	      for (k=SUMlo;k<=SUMhi;++k) {
		int_array[k]=pion_excitsigma(E0,OSCILLATOR,DELTANUD,ALPHA,E_array[k])*Labsorb[k];
		int_junk+=EBIN*int_array[k];
	      }
	      /*	    pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
			    int_ans=pion_integrate(pion_integrand,1.001*Elo,0.999*Ehi);*/
	      if (1 || int_ans<0.) int_ans=int_junk;
	      He_excite[element][LINE]=f_COVERING*Nion[element][electron]*int_ans;
	      /* Add line to Seyfert 2 spectrum */
	      E0=doppler_trans/doppler_rad*E0;
	      k=(int) ((E0-EMIN+0.5*EBIN)/EBIN);
	      if ((E0-EMIN+0.5*EBIN)/EBIN-(double) k >= 0.5) ++k;
	      if (k <= SPECBINS && k >= 1) {
		exc_spectrum[k]+=He_excite[element][LINE]/EBIN;
	      }
	    }
	  }
	}
      }  
      
      if (verbose) printf("Determining HIGH-n Photoexcitation Rates (n>5)\n");
      /* Hydrogenic and Heliumlike */
      for (i=1;i<=ELEMENTS;++i) {
	element=list[i];
	for (electron=1;electron<=2;++electron) {
	  if (Nion[element][electron]) {
	    if (electron==1) oscillatornorm=1.6; /* Bethe-Salpeter p. 265 */
	    if (electron==2) oscillatornorm=oshe[element]; /* defined above */
	    for (n=6;n<=HIGHN;++n) {
	      OSCILLATOR=oscillatornorm/cube((double) n);
	      E0=HC_KEV_ANGSTROM/highn[element][electron].lambda[n]*1000.;
	      E0=E0*doppler_rad;
	      if (E0>=EMIN && E0<=EMAX) {
		DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
		/* taking classical value for "ga=gamma" from p. 112-114 B+D */
		ga=DAMPING_CLASSICAL*sqr(E0*eVtoergs/hhh);
		ALPHA=ga/(4.*PI*DELTANUD);
		pion_line_limits(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,&Elo,&Ehi,&SUMlo,&SUMhi,0.9,1.1,1);
		int_junk=0.;
		for (k=SUMlo;k<=SUMhi;++k) {
		  int_array[k]=pion_excitsigma(E0,OSCILLATOR,DELTANUD,ALPHA,E_array[k])*Labsorb[k];
		  int_junk+=EBIN*int_array[k];
		}
		/*	      pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
			      int_ans=pion_integrate(pion_integrand,1.001*Elo,0.999*Ehi);*/
		if (1 || int_ans<0.) int_ans=int_junk;
		strength=f_COVERING*Nion[element][electron]*int_ans;
		/* add line to seyfert 2 like spectrum */
		E0=doppler_trans/doppler_rad*E0;
		k=(int) ((E0-EMIN+0.5*EBIN)/EBIN);
		if ((E0-EMIN+0.5*EBIN)/EBIN-(double) k >= 0.5) ++k;
		if (k <= SPECBINS && k >= 1) {
		  exc_spectrum[k]+=strength/EBIN;
		}
	      }
	    }
	  }
	}
      }
    }    

    /* Calculate Thomson scattering of spectrum */
    if (N_e) {
      for (k=1;k<=SPECBINS;++k) exc_spectrum[k]+=f_COVERING*N_e*sigmaT*Labsorb[k];
    }
    
    
    /* File root for all L-shell input files */
    for (element=10;element<=28;++element) {
      for (electron=3;electron<=10;++electron) {
	/* line contributions */
	if (Nion[element][electron]) {
	  if (verbose) printf("%2d %2d Recombination ...\n",element,electron);
	  /* Read in RR and DR contributions */
	  if (verbose) printf("Read in total RR and DR rates for %2d %2d\n",element,electron);
	  trates=pion_ad_trates(element,electron);
	  for (k=1;k<=RECNUM && trates;++k) {
	    Ttemp=trates[k-1].T;
	    RRtemp=trates[k-1].RR;
	    DRtemp=trates[k-1].DR;
	    L_kT[k]=Ttemp;
	    L_RR[k]=1.e-10*RRtemp;
	    L_DR[k]=1.e-10*DRtemp;
	    L_REC[k]=L_RR[k]+L_DR[k];
	  }

	  pion_spline(L_kT,L_RR,RECNUM,1.e40,1.e40,L_RR_2);
	  pion_spline(L_kT,L_DR,RECNUM,1.e40,1.e40,L_DR_2);
	  pion_spline(L_kT,L_REC,RECNUM,1.e40,1.e40,L_REC_2);
	  
	  kT=Tion[element][electron];
	  L_RR_kT=pion_L_RR_spline(kT);
	  L_DR_kT=pion_L_DR_spline(kT);
	  L_REC_kT=pion_L_REC_spline(kT);
	  
	  if (verbose) printf("PE rates %2d %2d ...\n",element,electron);
	  trsh=pion_ad_tr_shorter(element,electron,&ntrsh);
	  for (r=0;r<ntrsh;++r) {
	    j=trsh[r].j;
	    i=trsh[r].i;
	    en=trsh[r].en;
	    ftemp=trsh[r].f;
	    Atemp=trsh[r].A;
	    /*    sscanf(line,"%d%d%d%d%lf%lf%lf%lf%lf",&j,&jjunk,&i,&ijunk,&en,&ftemp,&Atemp,&DECAYRATE,&AIRATE);*/
	    j=j+1;
	    i=i+1;
	    g_i=((double) i)+1.;  /* deg=2J */
	    g_j=((double) j)+1.;  /* deg=2J */
	    ftemp=ftemp/g_i;     /* CHECK THIS - VERY IMPORTANT!!! */
	    E0=en;
	    E0=E0*doppler_rad;
	    DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
	    if ((E0>=EMIN && E0<=EMAX) && ftemp>tau_lim*FACTOR*DELTANUD/Nion[element][electron]) {
	      /*	    printf("%5d %2d %5d %2d %e %e %e\n",j,jjunk,i,ijunk,en,ftemp,Atemp);*/
	      OSCILLATOR=ftemp;
	      ga=Atemp/*g_j*Atemp/(3.*g_i*ftemp)*/;
	      ALPHA=ga/(4.*PI*DELTANUD);
	      pion_line_limits(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,&Elo,&Ehi,&SUMlo,&SUMhi,0.9,1.1,1);
	      int_junk=0.;
	      for (k=SUMlo;k<=SUMhi;++k) {
		int_array[k]=f_COVERING*Nion[element][electron]*pion_excitsigma(E0,OSCILLATOR,DELTANUD,ALPHA,E_array[k])*Labsorb[k];
		int_junk+=EBIN*int_array[k];
	      }      
	      /*	      pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
			      int_ans=pion_integrate(pion_integrand,1.001*Elo,0.999*Ehi);*/
	      if (1 || int_ans<0.) int_ans=int_junk;
	      ratePE=int_ans;
	      /* *Atemp/(DECAYRATE+AIRATE);
	       ratePI[element][electron]+=djunk*AIRATE/(DECAYRATE+AIRATE);*/
	      /*	      printf("line: PI[%2d][%2d] = %e    %e\n",element,electron,djunk*AIRATE/(DECAYRATE+AIRATE),E0);*/
	      /* Need to add line emission to overall spectrum */
	      E0=doppler_trans/doppler_rad*E0;
	      k=(int) ((E0-EMIN+0.5*EBIN)/EBIN);
	      if ((E0-EMIN+0.5*EBIN)/EBIN-(double) k >= 0.5) ++k;
	      if (k <= SPECBINS && k >= 1)
		exc_spectrum[k]+=ratePE/EBIN;
	    }
	  }
	  
	  if (verbose) printf("PI rates %2d %2d ... \n",element,electron);
	  pirec=pion_ad_fac_pi(PION_L_SHELL,element,electron,&npirec);
	  for (r=0;r<npirec;++r) {
	    /* PHOTOIONIZATION OUT OF GRD STATE OF n+1 ION */
	    pion_fac_load_record(&pirec[r],&g_i,&g_j);
	    p0=pirec[r].p[0]; p1=pirec[r].p[1]; p2=pirec[r].p[2]; p3=pirec[r].p[3];
	    if ((THRESHOLD>=EMIN/doppler_rad && THRESHOLD<=EMAX/doppler_rad) && pion_pisigma(g_i,p0,p1,p2,p3,THRESHOLD) >= 0. /* This one should be left at zero????*/) {
	      pion_fac_build_table(&pirec[r],g_i);

	      /* calculate PI rate and modify ratePI[element][electron] */
	      djunk=Nion[element][electron]*pion_fac_PI_rate_integral(THRESHOLD,Labsorb);
	      ratePI[element][electron]+=djunk;
	      /*printf("edge: PI[%2d][%2d]=%e    %e\n",element,electron,djunk,THRESHOLD*doppler_rad);*/
	    }
	  }
	  /*	  if (verbose) printf("%e (PI)\n",ratePI[element][electron]);*/

	  EMion[element][electron]=ratePI[element][electron]/L_REC_kT;
	  EM[element][electron]=1.2/ABUND[element]*EMion[element][electron];

	  rrl=pion_ad_rrlines(element,electron,&nrrl);
	  for (g=0;g<nrrl;g+=10) {
	    /* One group of 10 rows: kT, RR and DR per row; type, levels, energy
	       and wavelength as the last row of the group leaves them. */
	    for (k=1;k<=10;++k) {
	      L_kT[k]=rrl[g+k-1].kT;
	      RR_line[k]=rrl[g+k-1].RR;
	      DR_line[k]=rrl[g+k-1].DR;
	    }
	    typenum=rrl[g+9].typenum;
	    itemp=rrl[g+9].itemp;
	    jtemp=rrl[g+9].jtemp;
	    en=rrl[g+9].en;
	    /*    itemp=itemp+1;
		  jtemp=jtemp+1;*/
	    for (k=1;k<=10;++k) {
	      RR_line[k]*=1.e-10;
	      DR_line[k]*=1.e-10;
	    }
	    pion_spline(L_kT,DR_line,RECNUM,1.e40,1.e40,DR_line_2);	
	    pion_spline(L_kT,RR_line,RECNUM,1.e40,1.e40,RR_line_2);	
	    if (typenum<100) {
	      /*read in RRC contributions */
	      rrs=pion_ad_rr_short(element,electron,&nrrs);
	      for (r=0;r<nrrs;++r) {
		i=rrs[r].i;
		j=rrs[r].j;
		THRESHOLD=rrs[r].THRESHOLD;
		ANGULAR=rrs[r].ANGULAR;
		g_i=rrs[r].g_i+1.;
		g_j=rrs[r].g_j+1.;
		p0=rrs[r].p[0]; p1=rrs[r].p[1]; p2=rrs[r].p[2]; p3=rrs[r].p[3];
		for (k=1;k<=LOWE_GRIDNUM;++k) {
		  LOWE_EGRID[k]=rrs[r].grid[k-1][0];
		  LOWE_EGRID[k]=log10(LOWE_EGRID[k]/doppler_trans);
		}
		
		if (i==itemp && j==jtemp) {
		  /*		  printf("%5d %5d  %5d %5d  %e\n",i,itemp,j,jtemp,pion_rrsigma(g_i,g_j,p0,p1,p2,p3,THRESHOLD)); */
		  if (THRESHOLD>=EMIN/doppler_trans && THRESHOLD<=EMAX/doppler_trans) {
		    /* A 20000-point RRGRID and a 6-point LOWE_RRGRID spline used to be built
		       here. Nothing read either (pion_rrspline and pion_lowErrspline had no
		       callers), so both were removed. The LOWE_EGRID read above is kept: it
		       still sets that shared array. */
		    int_junk=0.;
		    for (k=1;k<=SPECBINS;++k) {
		      int_array[k]=pion_fac_recombination(g_i,g_j,p0,p1,p2,p3,kT,E_array[k]/doppler_trans-THRESHOLD);
		      int_junk+=EBIN*int_array[k];
		    }
		    /* CHANGE HERE */
		    intMIN=1.00001*THRESHOLD*doppler_trans;
		    intMAX=100.*THRESHOLD*doppler_trans;
		    if (intMAX>EMAX) intMAX=EMAX;
		    /*		    pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
				    int_ans=pion_integrate(pion_integrand,intMIN,intMAX);*/
		    if (1 || int_ans<0.) int_ans=int_junk;
		    RECNORM=int_ans;
		    rateRR[element][electron]=ratePI[element][electron]*pion_RR_line_spline(kT)/L_REC_kT;		  
		    for (k=1;k<=SPECBINS;++k) {
		      djunk=rateRR[element][electron]*pion_fac_recombination(g_i,g_j,p0,p1,p2,p3,kT,E_array[k]/doppler_trans-THRESHOLD)/RECNORM;
		      specRR[k]+=djunk;
		      rec_spectrum[k]+=djunk;
		    }
		  }
		}
	      }
	    } else if (typenum>=100 && en>=EMIN/doppler_trans && en<=EMAX/doppler_trans) {
	      /*read in line contributions */
	      E0=en;
	      E0=doppler_trans*E0;
	      rateRR[element][electron]=ratePI[element][electron]*pion_RR_line_spline(kT)/L_REC_kT;
	      rateDR[element][electron]=ratePI[element][electron]*pion_DR_line_spline(kT)/L_REC_kT;
	      k=(int) ((E0-EMIN+0.5*EBIN)/EBIN);
	      if ((E0-EMIN+0.5*EBIN)/EBIN-(double) k >= 0.5) ++k;
	      if (k <= SPECBINS && k >= 1) {
		specRR[k]+=rateRR[element][electron]/EBIN;
		specDR[k]+=rateDR[element][electron]/EBIN;
		rec_spectrum[k]+=specRR[k]+specDR[k];
	      }
	    }
	  }
	}
      }
    }

    if (verbose) printf("Determining Photoionization Rates for H- and He-like...\n");
    /* Photoionization Rates */
    for (i=1;i<=ELEMENTS;++i) {
      element=list[i];
      for (electron=1;electron<=2;++electron) {
	if (Nion[element][electron]) {
	  if (!have_verner(&vernerionizsigma[element][electron],element,electron)) continue;
	  THRESHOLD=vernerionizsigma[element][electron].Eth;
	  int_junk=0.;
	  for (k=1;k<=SPECBINS;++k) {
	    int_array[k]=pion_vernerph(vernerionizsigma[element][electron],E_array[k])*Labsorb[k];
	    int_junk+=EBIN*int_array[k];
	  }
	  /*	  pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
		  int_ans=pion_integrate(pion_integrand,1.001*THRESHOLD*doppler_rad,0.999*EMAX);*/
	  if (1 || int_ans<0.) int_ans=int_junk;
	  ratePI[element][electron]=f_COVERING*Nion[element][electron]*int_ans;
	}
      }
    }
    
    if (verbose) printf("Determining Recombination Line/RRC Strengths & EM's ...\n");
    if (verbose) printf("H-like...\n");
    H_rec=pion_ad_h_rec();
    for (i=1;i<=ELEMENTS;++i) {
      element=list[i];
      electron=1;
      if (Nion[element][electron]) {
	if (!have_verner(&vernerionizsigma[element][electron],element,electron)) continue;
	R=ratePI[element][electron];
	for (k=1;k<=PION_REC_TEMPERATURES;++k) {
	  Tvec[k]=log10(H_rec[element].T[k]);
	}
	/* Now calculate each line strength for given ion kT */
	for (LINE=1;LINE<=5;++LINE) {
	  for (k=1;k<=PION_REC_TEMPERATURES;++k) {
	    Yvec[k]=log10(H_rec[element].lines[LINE][k]);
	  }
	  pion_spline(Tvec,Yvec,PION_REC_TEMPERATURES,1.e40,1.e40,Yvec2);
	  strength=pow(10.,pion_loglinestrength(log10(Tion[element][electron]),PION_REC_TEMPERATURES));
	  E0=HC_KEV_ANGSTROM/hydrogen[element].lambda[LINE]*1000.;
	  E0=doppler_trans*E0;
	  k=(int) ((E0-EMIN+0.5*EBIN)/EBIN);
	  if ((E0-EMIN+0.5*EBIN)/EBIN-(double) k >= 0.5) ++k;
	  if (k >= 1 && k <= SPECBINS) rec_spectrum[k]+=R*strength/EBIN;
	}
	/* RRC strength */
	kT=Tion[element][electron];
	for (k=1;k<=PION_REC_TEMPERATURES;++k) Yvec[k]=log10(H_rec[element].rrc[k]);
	pion_spline(Tvec,Yvec,PION_REC_TEMPERATURES,1.e40,1.e40,Yvec2);
	strength=pow(10.,pion_loglinestrength(log10(Tion[element][electron]),PION_REC_TEMPERATURES));
	E0=HC_KEV_ANGSTROM/hydrogen[element].lambda[6/*rrc*/]*1000.;
	E0=doppler_trans*E0;
	/* normalize "recombination" */
	int_junk=0.;
	for (k=1;k<=SPECBINS;++k) {
	  int_array[k]=pion_verner_recombination(vernerionizsigma[element][electron],kT,E_array[k]/doppler_trans);
	  int_junk+=EBIN*int_array[k];
	}
	THRESHOLD=vernerionizsigma[element][electron].Eth;
	/* CHANGE HERE */
	intMIN=1.00001*THRESHOLD*doppler_trans;
	intMAX=100.*THRESHOLD*doppler_trans;
	if (intMAX>EMAX) intMAX=EMAX;
	/*pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
		int_ans=pion_integrate(pion_integrand,1.001*EMIN,0.999*EMAX);*/
	if (1 || int_ans<0.) int_ans=int_junk;
	RECNORM=int_ans;
	for (k=1;k<=SPECBINS;++k) {
	  rec_spectrum[k]+=R*strength/RECNORM*pion_verner_recombination(vernerionizsigma[element][electron],kT,E_array[k]/doppler_trans);
	}
	/* C-coefficient magnitude */
	for (k=1;k<=PION_REC_TEMPERATURES;++k) {
	  Yvec[k]=log10(H_rec[element].C[k]);
	}
	pion_spline(Tvec,Yvec,PION_REC_TEMPERATURES,1.e40,1.e40,Yvec2);
	strength=pow(10.,pion_loglinestrength(log10(Tion[element][electron]),PION_REC_TEMPERATURES));

	EMion[element][electron]=R/strength;
	EM[element][electron]=1.2/ABUND[element]*EMion[element][electron];
      }
    }
    
    /* Heliumlike */
    if (verbose) printf("He-like...\n");
    He_rec=pion_ad_he_rec();
    for (i=1;i<=ELEMENTS;++i) {
      element=list[i];
      electron=2;
      if (Nion[element][electron]) {
	if (!have_verner(&vernerionizsigma[element][electron],element,electron)) continue;
	R=ratePI[element][electron];
	for (k=1;k<=PION_REC_TEMPERATURES;++k) {
	  Tvec[k]=log10(He_rec[element].T[k]);
	}
	/* Now calculate each line strength for given ion kT */
	for (LINE=1;LINE<=7;++LINE) {
	  for (k=1;k<=PION_REC_TEMPERATURES;++k) {
	    Yvec[k]=log10(He_rec[element].lines[LINE][k]);
	  }
	  pion_spline(Tvec,Yvec,PION_REC_TEMPERATURES,1.e40,1.e40,Yvec2);
	  strength=pow(10.,pion_loglinestrength(log10(Tion[element][electron]),PION_REC_TEMPERATURES));
	  /* line.dat has no wavelength (0) for He-like LINE 7 of Ca, Fe and Ni.
	   * The energy would be infinite and its bin index an int overflow
	   * (undefined behaviour); the line was never added anyway, since that
	   * index fell outside the spectrum. */
	  if (helium[element].lambda[LINE] <= 0.) continue;
	  E0=HC_KEV_ANGSTROM/helium[element].lambda[LINE]*1000.;
	  E0=doppler_trans*E0;
	  k=(int) ((E0-EMIN+0.5*EBIN)/EBIN);
	  if ((E0-EMIN+0.5*EBIN)/EBIN-(double) k >= 0.5) ++k;
	  if (k <= SPECBINS && k >= 1) rec_spectrum[k]+=R*strength/EBIN;
	}
	/* RRC strength */
	kT=Tion[element][electron];
	for (k=1;k<=PION_REC_TEMPERATURES;++k) {
	  Yvec[k]=log10(He_rec[element].rrc[k]);
	}
	pion_spline(Tvec,Yvec,PION_REC_TEMPERATURES,1.e40,1.e40,Yvec2);
	strength=pow(10.,pion_loglinestrength(log10(Tion[element][electron]),PION_REC_TEMPERATURES));
	E0=HC_KEV_ANGSTROM/helium[element].lambda[9/*rrc*/]*1000.;
	E0=doppler_trans*E0;
	/* normalize "recombination" */
	int_junk=0.;
	for (k=1;k<=SPECBINS;++k) {
	  int_array[k]=pion_verner_recombination(vernerionizsigma[element][electron],kT,E_array[k]/doppler_trans);
	  int_junk+=EBIN*int_array[k];
	}
	THRESHOLD=vernerionizsigma[element][electron].Eth;
	/* CHANGE HERE */
	intMIN=1.00001*THRESHOLD*doppler_trans;
	intMAX=100.*THRESHOLD*doppler_trans;
	if (intMAX>EMAX) intMAX=EMAX;
	/*	pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
		int_ans=pion_integrate(pion_integrand,intMIN,intMAX);*/
	if (1 || int_ans<0.) int_ans=int_junk;
	RECNORM=int_ans;
	for (k=1;k<=SPECBINS;++k) {
	  rec_spectrum[k]+=R*strength/RECNORM*pion_verner_recombination(vernerionizsigma[element][electron],kT,E_array[k]/doppler_trans);
	}
	/* C-coefficient magnitude */
	for (k=1;k<=PION_REC_TEMPERATURES;++k) {
	  Yvec[k]=log10(He_rec[element].C[k]);
	}
	pion_spline(Tvec,Yvec,PION_REC_TEMPERATURES,1.e40,1.e40,Yvec2);
	strength=pow(10.,pion_loglinestrength(log10(Tion[element][electron]),PION_REC_TEMPERATURES));

	EMion[element][electron]=R/strength;
	EM[element][electron]=1.2/ABUND[element]*EMion[element][electron];
      }
    }
}

/* The tabulated input continuum (INPUT > 0): read input.qdp from the cwd
 * (a user file, so every evaluation), spline E*L and L, normalize to L_X, and
 * set Labsorb and abs_spectrum for the given type. The INPUT>0 branch of
 * photoion, phsi and phxi, identical in all three, moved verbatim. Returns 1
 * if input.qdp is missing (the models return 0 then, as before; their entry
 * probe catches it first), else 0. The arrays it allocates are the globals
 * E_input, L_input, ... which the models free. */
int pion_input_continuum(double redshift, int type, int verbose, double Labsorb[])
{
  const char *temp="input.qdp";
  char line[400];
  const int LENGTH=400;
  int k;
  double djunk;
  FILE *input;

    input=fopen(temp,"r");
    if (input==NULL) {
      printf("The file 'input.qdp' must exist in this directory.\n");
      return 1;
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
  return 0;
}

/* The report photoion, phsi and phxi print at the end of an evaluation when
 * verbose is set (the caller tests it): the ion table (type > 1 only), then
 * the radiation force on the absorbing gas and, for the power-law continuum
 * (INPUT == 0), its normalization at 1 keV. For information only; nothing
 * here enters the model output. The block was identical in all three models.
 * Reads the state globals (EM, Tion, EMion, tau, ...) and overwrites the
 * scratch array int_array. */
void pion_emission_report(int type, double redshift,
                          double **ratePI, const double ABUND[])
{
  int element,electron,i;
  int first=1;
  double int_ans;

  if (type>1) {
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

  /* Radiation force on the absorbing gas */
  if (INPUT==0) {
    for (i=1;i<=SPECBINS;++i) {
      int_array[i]=E_array[i]*eVtoergs/ccc*pion_L(E_array[i])*(1.-exp(-tau[i]));
    }
  } else {
    for (i=1;i<=SPECBINS;++i) {
      int_array[i]=E_array[i]*eVtoergs/ccc*pion_Linterp(E_array[i])*(1.-exp(-tau[i]));
    }
  }
  int_ans=0.;
  for (i=1;i<=SPECBINS;++i) int_ans+=EBIN*int_array[i];
  int_ans=f_COVERING*int_ans;
  if (first==1) printf("******************************************************************\n");
  printf("* Radiation Pressure = %4.2e dyne (%7.3lf - %7.3lf keV, rf) *\n",int_ans,EMIN/1000.,EMAX/1000.);
  if (INPUT==0) printf("* Power-law Norm at 1 keV: %e [photons/cm^2/s/keV]     *\n", LNORM/(1.+redshift)/1000./(4.*PI*sqr(D))*pow(1000.,2.-GAMMA));
  printf("******************************************************************\n");
}
