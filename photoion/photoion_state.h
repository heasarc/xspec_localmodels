/* Mutable module state shared by the photoion models.
 * See photoion_state.c. Not re-entrant.
 */

#ifndef PHOTOION_STATE_H
#define PHOTOION_STATE_H

#include <stdio.h>

extern double L_X;
extern double L_EMAX;
extern double L_EMIN;
extern double LNORM;
extern double LinterpNORM;
extern double GAMMA;
extern double f_COVERING;
extern double D;
extern double kT;
extern double sigmav_rad;
extern double sigmav_trans;
extern double v_rad;
extern double v_trans;
extern double *EGRID;
extern double *PIGRID;
extern double *PIGRID_2;
extern double *RRGRID;
extern double *RRGRID_2;
extern double ELECTRON_ENERGY;
extern double THRESHOLD;
extern double ANGULAR;
extern int GRIDNUM;
extern int RECNUM;
extern double *LOWE_EGRID;
extern double *LOWE_PIGRID;
extern double *LOWE_PIGRID_2;
extern double *LOWE_RRGRID;
extern double *LOWE_RRGRID_2;
extern int LOWE_GRIDNUM;
extern double *E_array;
extern double *E_spectrum;
extern double *abs_spectrum;
extern double *exc_spectrum;
extern double *rec_spectrum;
extern double *l_array;
extern double *l_spectrum;
extern double *int_array;
extern double *int_array_2;
extern double *Tvec;
extern double *Yvec;
extern double *Yvec2;
extern double *z_array;
extern double *hubble_array;
extern double *hubble_array_2;
extern int HUBBLE_BINS;
extern double EBIN;
extern double EMIN;
extern double EMAX;
extern double voigt_lim;
extern double tau_lim;
extern double pi_rate_lim;
extern int INPUT;
extern int INPUT_SIZE;
extern int INPUT_SHIFT;
extern double *E_input;
extern double *L_input;
extern double *L_input_2;
extern double *EtimesL_input;
extern double *EtimesL_input_2;
extern double HALFBIN_SIZE;
extern double N_e;
extern double **Nion;
extern double **Tion;
extern double **EMion;
extern double **EM;
extern double doppler_rad;
extern double doppler_trans;
extern double **Hionizsigmaconv;
extern double **Heionizsigmaconv;
extern double *ionizsigmatemp;
extern double *spectrumtemp;
extern double *tau;
extern double *tau_exc;
extern double *tau_edge;
extern double *RR_line;
extern double *RR_line_2;
extern double *DR_line;
extern double *DR_line_2;
extern double *L_kT;
extern double *L_RR;
extern double *L_RR_2;
extern double *L_DR;
extern double *L_DR_2;
extern double *L_REC;
extern double *L_REC_2;
extern double L_RR_kT;
extern double L_DR_kT;
extern double L_REC_kT;
extern int SPECBINS;
extern double *xi_frac_grid;
extern double *frac_grid;
extern double *xi_fion_grid;
extern double *fion_grid;
extern double *fion_grid_2;
extern double XIMIN;
extern double XIMAX;
extern double HNORM;
extern double *xi_array;
extern double *fion_array;
extern double *fion_array_2;
extern int FRACXINUM;
extern int FIONXINUM;
extern int XINUM;
extern int fion_integrate;

#endif /* PHOTOION_STATE_H */
