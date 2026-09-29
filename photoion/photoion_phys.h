/* Atomic physics, recombination and ion-fraction routines shared by the
 * photoion models. See photoion_phys.c.
 */

#ifndef PHOTOION_PHYS_H
#define PHOTOION_PHYS_H

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

double pion_DR_line_spline(double temp);
double pion_EtimesL(double E );
double pion_HeI_edge(double E);
double pion_L(double E );
double pion_L_DR_spline(double temp);
double pion_L_REC_spline(double temp);
double pion_L_RR_spline(double temp);
double pion_Linterp(double E );
double pion_RR_line_spline(double temp);
double pion_dfdE(double g_i,double p0, double p1, double p2, double p3, double Te);
double pion_doppler(double v );
double pion_excitsigma(double E0, double OSCILLATOR,double DELTANUD,double ALPHA, double E);
double pion_fac_PI_rate_integral(double THRESHOLD,double Labsorb[]);
double pion_fac_ionizsigma(double THRESHOLD, double E);
double pion_fac_recombination(double g_i,double g_j,double p0,double p1,double p2,double p3,double kT, double Te);
double pion_fion(double xi);
double pion_fion_integrand(double temp);
double pion_frac(double xi);
double pion_gauss(double s,double x);
double pion_hubble_integrand(double z );
double pion_hubble_integrate(double z );
double pion_integrand(double temp);
double pion_loglinestrength(double logkT, int ntemps);
double pion_lowEpispline(double E);
double pion_lowErrspline(double E);
double pion_maxwell(double Te, double kT, double NORM);
double pion_pisigma(double g_i, double p0, double p1, double p2, double p3, double E);
double pion_pispline(double E);
double pion_rrsigma(double g_i, double g_j, double p0, double p1, double p2, double p3, double Te );
double pion_rrspline(double E);
double pion_verner_recombination(struct VERNER_STRUCT vioniz, double kT, double E);
double pion_vernerpartialph(struct VERNER_PARTIAL_STRUCT verner, double E);
double pion_vernerph(struct VERNER_STRUCT verner, double E);
double pion_voigt(double alpha,double v);
void pion_HeI_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_p[]);
void pion_fac_edge_opacity(double Nion_column_density, double THRESHOLD, double tau_p[]);
void pion_verner_full_edge_opacity(double Nion_column_density, double THRESHOLD, struct VERNER_STRUCT verner, double tau_p[]);
void pion_verner_partial_edge_opacity(double Nion_column_density, double THRESHOLD, struct VERNER_PARTIAL_STRUCT verner, double tau_p[]);

void pion_line_limits(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double *Elo,double *Ehi,int *SUMlo,int *SUMhi,double pad_lo,double pad_hi,int clamp_both);
void pion_line_opacity(double Nion_column_density,double E0,double OSCILLATOR,double ALPHA,double DELTANUD,double tau_exc_p[],double pad_lo,double pad_hi,int clamp_both);

#endif /* PHOTOION_PHYS_H */
