/* Atomic physics, recombination and ion-fraction routines shared by the
 * photoion models.
 *
 * Every routine here was duplicated verbatim across the model files that use
 * it -- up to eight copies -- differing only by a name suffix. The one
 * exception is frac(), where phxi.c spelled a condition as two branches and
 * xiabs.c as one; they are equivalent and xiabs.c's form is used.
 *
 * NOT shared, deliberately: line_limits() and line_opacity() exist in two
 * genuinely different versions. The emission models (photoion, phsi, phxi)
 * carry a reworked pair that moves the bin-range calculation into
 * line_limits(), widens the integration range by 10% and clamps both ends;
 * the absorption models (miabs, neutral, siabs, vneutral, xiabs) still have
 * the original. That looks like an improvement that was never propagated, but
 * unifying them changes numbers, so it is a physics decision rather than a
 * refactoring one and both copies are left in place.
 *
 * These routines read and write the shared state declared in photoion_state.h
 * rather than taking it as arguments, as the originals did. See that header
 * for the re-entrancy caveat.
 *
 * loglinestrength() has the one changed signature: it takes the number of
 * temperature points as an argument instead of reading the TEMPERATURES macro,
 * which is 11 in the emission models and 30 in the absorption models.
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "photoion_alloc.h"
#include "photoion_spline.h"
#include "photoion_state.h"
#include "photoion_phys.h"
#include "photoion_const.h"
#include "photoion_atomdata.h"

#define sqr(X) ((X)*(X))
#define SMALL (1.e-6)
#define cube(X) ((X)*(X)*(X))

double pion_DR_line_spline(double temp)
{
  double answer;
  
  pion_splint(L_kT,DR_line,DR_line_2,RECNUM,temp,&answer);
  return answer;
}

double pion_EtimesL(double E /*[eV]*/)   /* [photon/s] */
{
  double answer;

  if (E>=L_EMIN && E<=L_EMAX) pion_splint(E_input,EtimesL_input,EtimesL_input_2,INPUT_SIZE,E,&answer);
  else answer=0.;
  return LinterpNORM*answer;
}

double pion_HeI_edge(double E) /* Yan, Sadeghpour, Dalgarno (1998) */
{
  double answer=0.,x;

  /* The fit (their eq. 14) holds at and above threshold only. Below it the
   * polynomial in x^(-1/2) keeps rising, and pion_HeI_edge_opacity starts at
   * 0.95*threshold, so without this test He I opacity began 1.2 eV under the
   * edge. */
  if (E < 24.58) return 0.;

  x=E/24.58;
  answer+=1.;
  answer+=-4.7416/pow(x,1./2.);
  answer+=14.8200/pow(x,2./2.);
  answer+=-30.8678/pow(x,3./2.);
  answer+=37.3584/pow(x,4./2.);
  answer+=-23.4585/pow(x,5./2.);
  answer+=5.9133/pow(x,6./2.);
  answer*=733.0/pow(E/1000.,3.5)*1.e-24;

  return answer;
}

void pion_HeI_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_p[])
{
  int k,klo,kth;
  double taujunk;
  
  klo=(int) ((0.95*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.05*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  taujunk=Nion_column_density*pion_HeI_edge(E_array[k]/doppler_rad);
  while ((taujunk > tau_lim || k <= kth) && k < SPECBINS) {
    k=k+1;
    taujunk=Nion_column_density*pion_HeI_edge(E_array[k]/doppler_rad);
    tau_p[k]+=taujunk;
  }
}

double pion_L(double E /*[eV]*/) /* [photons/s/eV] */
{
  double answer;

  if (E<=L_EMAX && E>=L_EMIN) answer=LNORM*pow(E,-GAMMA);
  else answer=0.;
  return answer;
}

double pion_L_DR_spline(double temp)
{
  double answer;
  
  pion_splint(L_kT,L_DR,L_DR_2,RECNUM,temp,&answer);
  return answer;
}

double pion_L_REC_spline(double temp)
{
  double answer;
  
  pion_splint(L_kT,L_REC,L_REC_2,RECNUM,temp,&answer);
  return answer;
}

double pion_L_RR_spline(double temp)
{
  double answer;
  
  pion_splint(L_kT,L_RR,L_RR_2,RECNUM,temp,&answer);
  return answer;
}

double pion_Linterp(double E /*[eV]*/)   /* [photons/s/eV] */
{
  double answer;

  if (E>=L_EMIN && E<=L_EMAX) pion_splint(E_input,L_input,L_input_2,INPUT_SIZE,E,&answer);
  else answer=0.;
  return LinterpNORM*answer;
}

double pion_RR_line_spline(double temp)
{
  double answer;
  
  pion_splint(L_kT,RR_line,RR_line_2,RECNUM,temp,&answer);
  return answer;
}

double pion_dfdE(double g_i,double p0, double p1, double p2, double p3, double Te)
{
  double x,y,answer,E;

  E=Te+THRESHOLD;
  x=(Te+p3)/p3;
  y=(1.+p2)/(sqrt(x)+p2);
  answer=E/(Te+p3)*p0*pow(x,-3.5-ANGULAR+0.5*p1)*pow(y,p1)/g_i;
  return answer;
}

double pion_doppler(double v /*[cm/s]*/) /*[unitless] red: v>0.,d(v)>1, blue: v<0.,d(v)<1*/
{
  double answer;

  answer=sqrt((1.-v/ccc)/(1.+v/ccc));
  return answer;
}

double pion_excitsigma(double E0, double OSCILLATOR,double DELTANUD,double ALPHA, double E)
{
  double answer;
  double V;

  V=(E-E0)/(hhh*DELTANUD*ergstoeV);  
  answer=PI*re*ccc*OSCILLATOR/sqrt(PI)/DELTANUD*pion_voigt(ALPHA,V)*doppler_rad;
  return answer;
}

double pion_fac_PI_rate_integral(double THRESHOLD,const double Labsorb[])
{
  int k,klo,khi,kth;
  double strength,pitemp=0.;
  double int_junk,int_ans;

  klo=(int) ((0.5*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.1*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  pi_rate_lim=pion_fac_ionizsigma(THRESHOLD,E_array[kth]/doppler_rad)*Labsorb[kth];
  /*  k=1;
  while (k<SPECBINS) {*/
  int_junk=0.;
  while ((k<=kth || pitemp>=1.e-4*pi_rate_lim) && k < SPECBINS) {
    k=k+1;
    pitemp=pion_fac_ionizsigma(THRESHOLD,E_array[k]/doppler_rad)*Labsorb[k];
    int_array[k]=pitemp;
    int_junk+=EBIN*pitemp;
  }
  khi=(int) ((double) k/2.);
  if (khi>SPECBINS) khi=SPECBINS;
  /*  pion_spline(E_array,int_array,SPECBINS,1.e40,1.e40,int_array_2);
      int_ans=pion_integrate(pion_integrand,(1.001)*THRESHOLD*doppler_rad,E_array[khi]);*/
  if (1 || int_ans<0.) int_ans=int_junk;
  strength=f_COVERING*int_ans;
  return strength;
}

void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_p[]) 
{
  int k,klo,kth;
  double taujunk;
  
  klo=(int) ((0.95*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.05*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  taujunk=Nion_column_density*pion_fac_ionizsigma(THRESHOLD,E_array[k]/doppler_rad);
  while ((taujunk > tau_lim || k <= kth) && k < SPECBINS) {
    k=k+1;
    taujunk=Nion_column_density*pion_fac_ionizsigma(THRESHOLD,E_array[k]/doppler_rad);
    tau_p[k]+=taujunk;
  }
}

double pion_fac_ionizsigma(double THRESHOLD, double E)
{
  double answer=0.;
  
  if (E>=THRESHOLD) {
    if (E<=pow(10.,LOWE_EGRID[1])) answer=pion_lowEpispline(E);
    else answer=pion_pispline(E);
  }
  return answer;
}

double pion_fac_recombination(double g_i,double g_j,double p0,double p1,double p2,double p3,double kT, double Te)
{
  double E,answer;
  
  E=Te-THRESHOLD;
  if (Te > 0.)   answer=sqrt(2.*Te/meeV)*pion_maxwell(Te,kT,1.)*sqr(E)/(meeV*Te)*pion_rrsigma(g_i,g_j,p0,p1,p2,p3,Te);
  else answer=0.;
  return answer;
}

double pion_fion(double xi)
{
  double answer;
  
  pion_splint(xi_fion_grid,fion_grid,fion_grid_2,FIONXINUM,xi,&answer);
  return answer;
}

double pion_fion_integrand(double temp)
{
  double answer;

  pion_splint(xi_array,fion_array,fion_array_2,XINUM,temp,&answer);
  return answer;
}

double pion_frac(double xi)
{
  int i;
  double slope,answer;
  
  i=1;
  while (xi_frac_grid[i]<xi) ++i;
  if (i==1 || i>FRACXINUM) answer=0.;
  else {
    slope=(frac_grid[i]-frac_grid[i-1])/(xi_frac_grid[i]-xi_frac_grid[i-1]);
    answer=slope*(xi-xi_frac_grid[i-1])+frac_grid[i-1];
  }

  return answer;
}

double pion_gauss(double s,double x) 
{
  double answer=0.;

  answer=1./sqrt(2.*PI)/s*exp(-sqr((x)/s)/2.);  
  return answer;
}

/* XSPEC's cosmology, as set by its "cosmo" command (libXSFunctions, the
 * library that also provides FGMSTR). H0 is in km/s/Mpc. */
float csmgh0(void);
float csmgl0(void);

/* Read XSPEC's cosmology for the Hubble-law distance. The universe is taken
 * to be flat, Omega_m = 1 - lambda0, which is also how XSPEC itself treats a
 * non-zero lambda0; q0 is not used. lambda0 = 0 is therefore Einstein-de
 * Sitter. lambda0 >= 1 would leave Omega_m <= 0 and the integrand undefined,
 * so it is refused, as is H0 <= 0. Returns 0 on success. */
int pion_cosmology(double *H0, double *Omega_m)
{
  double h0=csmgh0(), lambda0=csmgl0();
  if (!(h0 > 0.) || !(lambda0 < 1.)) {
    printf("PHOTOION: cannot use the cosmology H0 = %g, lambda0 = %g for the "
           "Hubble-law distance (D = 0). This model needs H0 > 0 and "
           "lambda0 < 1, because it takes the universe to be flat with "
           "Omega_m = 1 - lambda0. Set it with XSPEC's cosmo command, or give "
           "a distance D.\n", h0, lambda0);
    return 1;
  }
  *H0=h0;
  *Omega_m=1.-lambda0;
  return 0;
}

/* Comoving distance integrand c/H(z) [cm] for a flat universe; H0 in
 * km/s/Mpc. Peacock, Eq. 3.39, p. 76 */
double pion_hubble_integrand(double z /* redshift */, double H0, double Omega_m)
{
  double answer;
  double H0_per_s=H0*1.e5/(1.e6*parsectocm);
  answer=ccc/H0_per_s/sqrt(Omega_m*cube(1.+z)+(1.-Omega_m));
  return answer;
}

double pion_hubble_integrate(double z /* redshift */)
{
  double answer;
  pion_splint(z_array,hubble_array,hubble_array_2,HUBBLE_BINS,z,&answer);
  return answer;
}

double pion_integrand(double temp)
{
  double answer;

  pion_splint(E_array,int_array,int_array_2,SPECBINS,temp,&answer);
  return answer;
}

double pion_loglinestrength(double logkT, int ntemps)
{
  double answer;

  pion_splint(Tvec,Yvec,Yvec2,ntemps,logkT,&answer);
  return answer;
}

double pion_lowEpispline(double E)
{
  double answer;
  
  E=log10(E);
  pion_splint(LOWE_EGRID,LOWE_PIGRID,LOWE_PIGRID_2,LOWE_GRIDNUM,E,&answer);
  return pow(10.,answer);
}

double pion_maxwell(double Te, double kT, double NORM)
{
  double answer;

  answer=NORM*sqrt(Te)*exp(-Te/kT);
  return answer;
}

double pion_pisigma(double g_i, double p0, double p1, double p2, double p3, double E)
{
  double answer;

  answer=sqr(a0)*2.*PI*FINE_STRUCTURE*pion_dfdE(g_i,p0,p1,p2,p3,E-THRESHOLD);
  return answer;
}

double pion_pispline(double E)
{
  double answer;
  
  E=log10(E);
  pion_splint(EGRID,PIGRID,PIGRID_2,GRIDNUM,E,&answer);
  return pow(10.,answer);
}

double pion_rrsigma(double g_i, double g_j, double p0, double p1, double p2, double p3, double Te /*electron energy?*/)
{
  double answer,E;

  E=Te+THRESHOLD;
  answer=sqr(FINE_STRUCTURE)/2.*g_i/g_j*(sqr(E)/Te*(eVtoHartree))*pion_pisigma(g_i,p0,p1,p2,p3,E);
  return answer;
}

void pion_verner_full_edge_opacity(double Nion_column_density, double THRESHOLD, struct VERNER_STRUCT verner, double tau_p[]) 
{
  int k,klo,kth;
  double taujunk;
  
  klo=(int) ((0.95*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.05*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  taujunk=Nion_column_density*pion_vernerph(verner, E_array[k]/doppler_rad);
  while ((taujunk > tau_lim || k <= kth) && k < SPECBINS) {
    k=k+1;
    taujunk=Nion_column_density*pion_vernerph(verner, E_array[k]/doppler_rad);
    tau_p[k]+=taujunk;
  }
}

void pion_verner_partial_edge_opacity(double Nion_column_density, double THRESHOLD, struct VERNER_PARTIAL_STRUCT verner, double tau_p[]) 
{
  int k,klo,kth;
  double taujunk;
  
  klo=(int) ((0.95*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  kth=(int) ((1.05*THRESHOLD*doppler_rad-EMIN+0.5*EBIN)/EBIN);
  if (klo < 1) klo=1;
  if (klo > SPECBINS) klo=SPECBINS;
  k=klo;
  taujunk=Nion_column_density*pion_vernerpartialph(verner, E_array[k]/doppler_rad);
  while ((taujunk > tau_lim || k <= kth) && k < SPECBINS) {
    k=k+1;
    taujunk=Nion_column_density*pion_vernerpartialph(verner, E_array[k]/doppler_rad);
    tau_p[k]+=taujunk;
  }
}

double pion_verner_recombination(struct VERNER_STRUCT vioniz, double kT, double E)
{
  double Te,THRESHOLD,answer;
  
  THRESHOLD=vioniz.Eth;
  Te=E-THRESHOLD;
  if (Te > 0.)   answer=/*grounddeg[electron]/(ge*freedeg[electron])*ccc*/sqrt(2.*Te/meeV)*pion_maxwell(Te,kT,1.)*sqr(E)/(meeV*Te)*pion_vernerph(vioniz,E);
  else answer=0.;
  return answer;
}

double pion_vernerpartialph(struct VERNER_PARTIAL_STRUCT verner, double E)
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

double pion_vernerph(struct VERNER_STRUCT verner, double E)
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

double pion_voigt(double alpha,double v)
{
  int i;
  /* Unit-offset coefficient tables; element 0 is unused. These were three
   * pion_dvector(1,7) allocated and freed on every call, which made calloc/free
   * a measurable share of line-opacity time. Fixed tables give the same values. */
  static const double a[8]={0.,
    122.607931777104326, 214.382388694706425, 181.928533092181549,
    93.155580458134410, 30.180142196210589, 5.912626209773153,
    0.564189583562615};
  static const double b[8]={0.,
    122.607931773875350, 352.730625110963558, 457.334478783897737,
    348.703917719495792, 170.354001821091472, 53.992906912940207,
    10.479857114260399};
  static const double c[8]={0.,
    0.5641641, 0.8718681, 1.474395, -19.57862, 802.4513, -4850.316,
    8031.468};
  double v2,v3,fac1,fac2;
  double p1,p2,p3,p4,p5,p6,p7;
  double o1,o2,o3,o4,o5,o6,o7;
  double q1,q2;
  double r1,r2;
  double H;
 
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

  return H;
}

/* Line integration limits, shared by all models.
 *
 * Solves for the detuning v_lim at which the line's optical depth falls below
 * tau_lim (using the Lorentz-wing limit of the Voigt function), or, if the
 * profile is still above voigt_lim there, for where it falls to voigt_lim of
 * line centre. Returns that interval both as energies [Elo,Ehi] and as the
 * spectrum bin range [SUMlo,SUMhi] the callers loop over. tau_lim and
 * voigt_lim are per-model state; see photoion_state.h.
 *
 * pad_lo and pad_hi scale Elo and Ehi before converting to bin indices. The
 * absorption models pass 1.0 and 1.0 and get the interval unchanged; the
 * emission models pass 0.9 and 1.1, widening the range so that their
 * fixed-grid sum does not truncate the wings of an integrated quantity. Note
 * the scaling is of absolute energy, not of line width -- at typical X-ray
 * line energies 0.9/1.1 is a far wider margin than the line itself.
 *
 * clamp_both selects which clamping the caller had before this was shared, and
 * the two differ only for a line lying entirely outside [EMIN,EMAX]. With
 * clamp_both false (absorption models) SUMlo may exceed SUMhi, so the caller's
 * loop runs zero times, which is correct: the line does not reach the grid.
 * With it true (emission models) both indices are pinned to the nearest edge,
 * so the caller deposits one bin of far-wing opacity at the boundary for a line
 * that is not in the band at all. Only the two clamps applied in both branches
 * are needed to keep the index in range; the extra pair is not a safety
 * measure. Preserved as-is so that sharing these routines changes no results;
 * which behavior is wanted is a separate question.
 */
void pion_line_limits(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double *Elo,double *Ehi,int *SUMlo,int *SUMhi,double pad_lo,double pad_hi,int clamp_both)
{
  double Vlim;

  Vlim=sqrt(ALPHA*Nion_column_density*re*ccc*OSCILLATOR/DELTANUD/tau_lim);
  if (pion_voigt(ALPHA,Vlim)>voigt_lim) {
    Vlim=sqrt(ALPHA/sqrt(PI)/(voigt_lim*pion_voigt(ALPHA,0.)));
  }
  if (E0-Vlim*(hhh*DELTANUD*ergstoeV)>EMIN) {
    *Elo=E0-Vlim*(hhh*DELTANUD*ergstoeV);
  } else *Elo=(1.+SMALL)*EMIN;
  if (E0+Vlim*(hhh*DELTANUD*ergstoeV)<EMAX) {
    *Ehi=E0+Vlim*(hhh*DELTANUD*ergstoeV);
  } else *Ehi=(1.-SMALL)*EMAX;

  *SUMlo=(int) ((pad_lo*(*Elo)-EMIN+0.5*EBIN)/EBIN);
  if ((pad_lo*(*Elo)-EMIN+0.5*EBIN)/EBIN-(double) *SUMlo >= 0.5) ++(*SUMlo);
  *SUMhi=(int) ((pad_hi*(*Ehi)-EMIN+0.5*EBIN)/EBIN);
  if ((pad_hi*(*Ehi)-EMIN+0.5*EBIN)/EBIN-(double) *SUMlo >= 0.5) ++(*SUMhi);

  if (clamp_both) {
    if (*SUMlo > SPECBINS) *SUMlo=SPECBINS;
    else if (*SUMlo < 1) *SUMlo=1;
    if (*SUMhi > SPECBINS) *SUMhi=SPECBINS;
    else if (*SUMhi < 1) *SUMhi=1;
  } else {
    if (*SUMlo < 1) *SUMlo=1;
    if (*SUMhi > SPECBINS) *SUMhi=SPECBINS;
  }
}

/* Accumulate one line's excitation opacity over the range pion_line_limits
 * returns. pad_lo and pad_hi are passed straight through. */
void pion_line_opacity(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double tau_exc_p[],double pad_lo,double pad_hi,int clamp_both)
{
  double Ehi,Elo;
  int k,SUMlo,SUMhi;

  pion_line_limits(Nion_column_density,E0,OSCILLATOR,ALPHA,DELTANUD,&Elo,&Ehi,&SUMlo,&SUMhi,pad_lo,pad_hi,clamp_both);

  for (k=SUMlo;k<=SUMhi;++k) {
    tau_exc_p[k]+=Nion_column_density*pion_excitsigma(E0,OSCILLATOR,DELTANUD,ALPHA,E_array[k]);
  }
}

/* Verner photoionization edges: the full cross sections for H- and He-like
 * ions (He I from Yan et al. instead), and the partial L- and M-shell cross
 * sections. This was the same block in seven models. zmax is the last element
 * in the H-/He-like loop, 28 in most models and 26 in siabs and xiabs; the two
 * agree because verner_photo.dat stops at Z=26. The partial-edge loops run to
 * the number of records each element actually has (open issue 12: they used to
 * scan all 126 slots of an uninitialized array, one past its end). Sets the
 * global THRESHOLD exactly as the inline copies did.
 *
 * Nion is a parameter, not the global: the five absorption models declare a
 * local Nion (and N_e, sigmav_rad, v_rad) that shadows the one in
 * photoion_state.h, so the global is NULL there. Shared code must never read
 * those four globals. */
void pion_verner_edges(double **Nion,
                       const struct VERNER_STRUCT (*vernerionizsigma)[31],
                       const struct VERNER_PARTIAL_STRUCT (*partialsigma)[125],
                       const int npartial[], int zmax, double tau_p[])
{
  int element, electron, j;

  /* Photoionization opacity for H- and He-like */
  for (element=1;element<=zmax;++element) {
    for (electron=1;electron<=2;++electron) {
      if (Nion[element][electron]) {
	if (element==2 && electron==2) {
	  THRESHOLD=24.58;
	  pion_HeI_edge_opacity(Nion[element][electron],THRESHOLD,tau_p);
	} else {
	  THRESHOLD=vernerionizsigma[element][electron].Eth;
	  pion_verner_full_edge_opacity(Nion[element][electron],THRESHOLD,vernerionizsigma[element][electron],tau_p);
	}
      }
    }
  }

  /* From Verner table: L-shell edges for C,N,O */
  for (element=6;element<=8;++element) {
    for (electron=3;electron<=8;++electron) {
      if (Nion[element][electron]) {
	for (j=0;j<npartial[element];++j) {
	  if (partialsigma[element][j].electron==electron && partialsigma[element][j].principal>=2) {
	    THRESHOLD=partialsigma[element][j].Eth;
	    pion_verner_partial_edge_opacity(Nion[element][electron],THRESHOLD,partialsigma[element][j],tau_p);
	  }
	}
      }
    }
  }

  /* From Verner table: Get L-shell edges for C,N,O and M-shell edges for M-shell ions */
  for (element=1;element<=28;++element) {
    for (electron=11;electron<=28;++electron) {
      if (Nion[element][electron] && !(electron <=20 && (element == 26 || element == 28))) {
	for (j=0;j<npartial[element];++j) {
	  if (partialsigma[element][j].electron==electron && partialsigma[element][j].principal>=3) {
	    THRESHOLD=partialsigma[element][j].Eth;
	    pion_verner_partial_edge_opacity(Nion[element][electron],THRESHOLD,partialsigma[element][j],tau_p);
	  }
	}
      }
    }
  }
}

/* Low-n (n <= 5) photoexcitation opacity of H- and He-like ions, from the
 * line.dat rows in file order: LINE <= 4 for H-like and LINE <= 6 for He-like,
 * where the ion has a column and the line an oscillator strength. This was the
 * same loop in all eight models, fused with the file read. Only pad_lo,
 * pad_hi and clamp_both differed (emission 0.9, 1.1, 1; absorption 1.0, 1.0, 0;
 * see pion_line_limits). The callers test `lines`. Nion and sigmav_rad are
 * parameters because the absorption models shadow both globals. */
void pion_lown_line_opacity(double **Nion, double sigmav_rad,
                            const struct PION_LINE_ROW *rows, int nrows,
                            double tau_p[], double pad_lo, double pad_hi, int clamp_both)
{
  int r, element, electron, LINE;
  double E0, OSCILLATOR, Atemp, ga, DELTANUD, ALPHA;

  for (r=0;r<nrows;++r) {
    element=rows[r].element;
    electron=rows[r].electron;
    LINE=rows[r].LINE;
    if (electron == 1) {
      if (!(Nion[element][electron] && rows[r].f && LINE <= 4)) continue;
    } else if (electron == 2) {
      if (!(Nion[element][electron] && rows[r].f && LINE <= 6)) continue;
    } else continue;
    E0=1000.*HC_KEV_ANGSTROM/rows[r].WAVE;
    E0=E0*doppler_rad;
    OSCILLATOR=rows[r].f;
    Atemp=rows[r].A;
    ga=Atemp;
    DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
    ALPHA=ga/(4.*PI*DELTANUD);
    pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_p,pad_lo,pad_hi,clamp_both);
  }
}

/* High-n (6 <= n <= HIGHN) photoexcitation opacity of H- and He-like ions, for
 * the elements in list[1..nlist]. Oscillator strengths scale as n^-3 (1.6 for
 * H-like, Bethe & Salpeter p. 265; oshe for He-like); Ni uses Fe's
 * wavelengths, scaled. Same in all eight models apart from the pad arguments. */
void pion_highn_line_opacity(double **Nion, double sigmav_rad,
                             const struct HIGHER_ORDER_STRUCT (*highn)[3],
                             const int list[], int nlist, int HIGHN, const double oshe[],
                             double tau_p[], double pad_lo, double pad_hi, int clamp_both)
{
  int i, n, element, electron;
  double E0=0., OSCILLATOR, oscillatornorm=0., ga, DELTANUD, ALPHA;

  for (i=1;i<=nlist;++i) {
    element=list[i];
    for (electron=1;electron<=2;++electron) {
      if (Nion[element][electron]) {
	if (electron==1) oscillatornorm=1.6; /* Bethe-Salpeter p. 265 */
	if (electron==2) oscillatornorm=oshe[element];
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
	  pion_line_opacity(Nion[element][electron],E0,OSCILLATOR,ALPHA,DELTANUD,tau_p,pad_lo,pad_hi,clamp_both);
	}
      }
    }
  }
}

/* Put one FAC pi_short record where the photoionization table code reads it:
 * the globals THRESHOLD and ANGULAR (pion_pisigma and pion_dfdE use both),
 * LOWE_EGRID/LOWE_PIGRID with their log10 transforms, and the 6-point spline.
 * g_i and g_j come back with the +1 every caller applied. This is the body of
 * the models' read loop, after the fscanf calls. */
void pion_fac_load_record(const struct PION_FAC_PIREC *r, double *g_i, double *g_j)
{
  int k;

  if (LOWE_GRIDNUM != PION_LOWE_N) pion_error("pion_fac_load_record: LOWE_GRIDNUM is not 6");
  THRESHOLD=r->THRESHOLD;
  ANGULAR=r->ANGULAR;
  *g_i=r->g_i+1.;
  *g_j=r->g_j+1.;
  for (k=1;k<=LOWE_GRIDNUM;++k) {
    LOWE_EGRID[k]=r->grid[k-1][0];
    LOWE_PIGRID[k]=r->grid[k-1][2];
    LOWE_EGRID[k]=log10(LOWE_EGRID[k]+THRESHOLD);
    LOWE_PIGRID[k]=log10(1.e-20*LOWE_PIGRID[k]);
  }
  pion_spline(LOWE_EGRID,LOWE_PIGRID,LOWE_GRIDNUM,1.e40,1.e40,LOWE_PIGRID_2);
}

/* The 20000-point photoionization table for the record loaded last: zero
 * (1e-90) below threshold, the 6-point spline up to its last node (open issue
 * 7 covers the region below its first node), the FAC fit above. Then log10
 * and a spline in PIGRID/PIGRID_2. */
void pion_fac_build_table(double g_i, const double p[4])
{
  int k;
  double p0=p[0], p1=p[1], p2=p[2], p3=p[3];

  /* PHOTOIONIZATION OUT OF GRD STATE OF n+1 ION */
  for (k=1;k<=GRIDNUM;++k) {
    if (EGRID[k]<log10(THRESHOLD)) PIGRID[k]=1.e-90;
    else if (EGRID[k]>=log10(THRESHOLD) && EGRID[k]<=LOWE_EGRID[LOWE_GRIDNUM]) PIGRID[k]=pion_lowEpispline(pow(10.,EGRID[k]));
    else PIGRID[k]=pion_pisigma(g_i,p0,p1,p2,p3,pow(10.,EGRID[k]));
  }
  for (k=1;k<=GRIDNUM;++k) PIGRID[k]=log10(PIGRID[k]);
  pion_spline(EGRID,PIGRID,GRIDNUM,1.e40,1.e40,PIGRID_2);
}

/* Edges of one FAC pi_short file: build each record's table if the edge
 * clears tau_lim at gate_factor*threshold, and add its opacity. */
static void fac_edges(double Nion_ion, int shell, int Z, int nelec, double gate_factor, double tau_edge_p[])
{
  int r, n;
  double g_i, g_j;
  const struct PION_FAC_PIREC *rec = pion_ad_fac_pi(shell, Z, nelec, &n);

  for (r=0;r<n;++r) {
    pion_fac_load_record(&rec[r], &g_i, &g_j);
    if (Nion_ion*pion_pisigma(g_i,rec[r].p[0],rec[r].p[1],rec[r].p[2],rec[r].p[3],gate_factor*THRESHOLD) >= tau_lim) {
      pion_fac_build_table(g_i, rec[r].p);
      /* calculate opacity of given edge and modify "tau" accordingly */
      pion_fac_edge_opacity(Nion_ion,THRESHOLD,tau_edge_p);
    }
  }
}

/* Lines of one FAC tr_short file inside [EMIN, EMAX] that clear tau_lim. */
static void fac_lines(double Nion_ion, double sigmav_rad, int shell, int Z, int nelec,
                      double tau_exc_p[], double pad_lo, double pad_hi, int clamp_both)
{
  int r, n;
  double g_i, ftemp, E0, DELTANUD, OSCILLATOR, ga, ALPHA;
  const struct PION_FAC_TRROW *row = pion_ad_fac_tr(shell, Z, nelec, &n);

  for (r=0;r<n;++r) {
    g_i=row[r].g_i+1.;
    ftemp=row[r].f/g_i;     /* CHECK THIS - VERY IMPORTANT!!! */
    E0=row[r].en;
    E0=E0*doppler_rad;
    DELTANUD=sqrt(2.)*sigmav_rad/ccc*(E0*eVtoergs/hhh);
    if ((E0>=EMIN && E0<=EMAX) && ftemp>tau_lim*FACTOR*DELTANUD/Nion_ion) {
      OSCILLATOR=ftemp;
      ga=row[r].A;
      ALPHA=ga/(4.*PI*DELTANUD);
      pion_line_opacity(Nion_ion,E0,OSCILLATOR,ALPHA,DELTANUD,tau_exc_p,pad_lo,pad_hi,clamp_both);
    }
  }
}

/* The three FAC shell blocks that were copied into every model:
 *   PION_FAC_L_NE_NI  L shell, Ne..Ni, 3-10 electrons (and Ni's 1-2)
 *   PION_FAC_L_C_O    L shell, C..O,  3-8 electrons
 *   PION_FAC_M        M shell, Mg..Ni, 11-28 electrons; lines before edges,
 *                     and the edge gate at 1.01*threshold
 * do_edges is 0 only in neutral, which takes its edges from neutral.tau.
 * do_lines is the model's `lines`. Edges add to tau_edge_p and lines to
 * tau_exc_p, each in the original order. Nion and sigmav_rad are parameters
 * because the absorption models shadow both globals. */
void pion_fac_shell_opacity(int which, double **Nion, double sigmav_rad,
                            int do_edges, int do_lines, int verbose,
                            double tau_edge_p[], double tau_exc_p[],
                            double pad_lo, double pad_hi, int clamp_both)
{
  int element, electron;

  if (which == PION_FAC_L_NE_NI) {
    for (element=10;element<=28;++element) {
      for (electron=1;electron<=10;++electron) {
	if (Nion[element][electron] && (electron>=3 || element==28)) {
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  if (do_edges) fac_edges(Nion[element][electron],PION_L_SHELL,element,electron,1.0,tau_edge_p);
	  if (do_lines && electron>=3)
	    fac_lines(Nion[element][electron],sigmav_rad,PION_L_SHELL,element,electron,tau_exc_p,pad_lo,pad_hi,clamp_both);
	}
      }
    }
  } else if (which == PION_FAC_L_C_O) {
    for (element=6;element<=8;++element) {
      for (electron=3;electron<=8;++electron) {
	if (Nion[element][electron]) {
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  if (do_edges) fac_edges(Nion[element][electron],PION_L_SHELL,element,electron,1.0,tau_edge_p);
	  if (do_lines)
	    fac_lines(Nion[element][electron],sigmav_rad,PION_L_SHELL,element,electron,tau_exc_p,pad_lo,pad_hi,clamp_both);
	}
      }
    }
  } else if (which == PION_FAC_M) {
    for (element=12;element<=28;++element) {
      for (electron=11;electron<=28;++electron) {
	if (Nion[element][electron]) {
	  if (verbose) printf("Z = %2d   z = %2d\n",element,electron);
	  if (do_lines)
	    fac_lines(Nion[element][electron],sigmav_rad,PION_M_SHELL,element,electron,tau_exc_p,pad_lo,pad_hi,clamp_both);
	  if (do_edges) fac_edges(Nion[element][electron],PION_M_SHELL,element,electron,1.01,tau_edge_p);
	}
      }
    }
  }
}
