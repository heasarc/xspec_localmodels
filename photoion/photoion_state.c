/* Mutable module state shared by the photoion models.
 *
 * Each model file used to declare its own private copy of these as file-scope
 * statics. The declarations agreed wherever they appeared: of the 104 names,
 * only HNORM differed, and only in carrying an explicit = 0., which is what a
 * file-scope object gets anyway.
 *
 * They are gathered here so the physics routines can be shared rather than
 * copied into every model. This is genuinely shared mutable state: XSPEC
 * evaluates one model component at a time and each model sets what it needs
 * at the start of its own evaluation, so sequential use is safe, but nothing
 * here is re-entrant.
 *
 * Names are deliberately left unprefixed. Prefixing would mean rewriting every
 * reference in the model files, where a local variable shadowing one of these
 * names would be silently captured -- a change that still compiles and links
 * but means something different.
 */

#include <stdio.h>

#include "photoion_state.h"

double L_X;
double L_EMAX;
double L_EMIN;
double LNORM;
double LinterpNORM;
double GAMMA;
double f_COVERING;
double D;
double kT;
double sigmav_rad;
double sigmav_trans;
double v_rad;
double v_trans;
double *EGRID;
double *PIGRID;
double *PIGRID_2;
double ELECTRON_ENERGY;
double THRESHOLD;
double ANGULAR;
int GRIDNUM=20000;
int RECNUM=10;
double *LOWE_EGRID;
double *LOWE_PIGRID;
double *LOWE_PIGRID_2;
int LOWE_GRIDNUM=6;
double *E_array;
double *E_spectrum;
double *abs_spectrum;
double *exc_spectrum;
double *rec_spectrum;
double *l_array;
double *l_spectrum;
double *int_array;
double *int_array_2;
double *Tvec;
double *Yvec;
double *Yvec2;
double *z_array;
double *hubble_array;
double *hubble_array_2;
int HUBBLE_BINS=100000;
double EBIN;
double EMIN;
double EMAX;
double voigt_lim;
double tau_lim;
double pi_rate_lim;
int INPUT;
int INPUT_SIZE;
int INPUT_SHIFT;
double *E_input;
double *L_input;
double *L_input_2;
double *EtimesL_input;
double *EtimesL_input_2;
double HALFBIN_SIZE;
double N_e;
double **Nion;
double **Tion;
double **EMion;
double **EM;
double doppler_rad;
double doppler_trans;
double **Hionizsigmaconv;
double **Heionizsigmaconv;
double *ionizsigmatemp;
double *spectrumtemp;
double *tau;
double *tau_exc;
double *tau_edge;
double *RR_line;
double *RR_line_2;
double *DR_line;
double *DR_line_2;
double *L_kT;
double *L_RR;
double *L_RR_2;
double *L_DR;
double *L_DR_2;
double *L_REC;
double *L_REC_2;
double L_RR_kT;
double L_DR_kT;
double L_REC_kT;
int SPECBINS;
double *xi_frac_grid;
double *frac_grid;
double *xi_fion_grid;
double *fion_grid;
double *fion_grid_2;
double XIMIN=-9.999;
double XIMAX=+9.999;
double HNORM=0.;
double *xi_array;
double *fion_array;
double *fion_array_2;
int FRACXINUM;
int FIONXINUM;
int XINUM=10000;
int fion_integrate;
