/* Atomic data from photoion_dat/, parsed once per process and kept.
 *
 * Every model used to re-read these files on every evaluation, in its own
 * copy of the same code. Each accessor here loads its file on first use and
 * hands out the cached result. The cache is tied to one PHOTOION_DIR: a model
 * entry that names a different directory drops everything and reloads lazily.
 * Like the spline cache and the integration workspace, it is allocated once
 * and deliberately never freed.
 *
 * Missing files are reported, not fatal. A loader that cannot open its file
 * prints the path once, leaves its data zero and raises a flag. Each model
 * runs to its normal cleanup and then, if the flag is set, zeroes its output
 * and returns 1. There is no early return, so nothing it allocated can leak.
 * Nothing here is re-entrant, like the rest of the package.
 */
#ifndef PHOTOION_ATOMDATA_H
#define PHOTOION_ATOMDATA_H

#include "photoion_phys.h"   /* struct VERNER_STRUCT, struct VERNER_PARTIAL_STRUCT */

/* Called on entry to every model, with the PHOTOION_DIR it was given. */
void pion_ad_begin(const char *datadir);

/* Nonzero if any data file could not be opened since pion_ad_begin. */
int pion_ad_failed(void);

/* abundance.dat into ABUND[1..30], zero where the file has no entry. */
void pion_ad_abundance(double ABUND[]);

/* oscillator_he.dat, scaled by 10^3, into oshe[] at the elements listed.
 * Other entries are left as they are. */
void pion_ad_oscillator_he(double oshe[]);

/* temperature.dat into Tion[Z][electrons] at the ions listed. */
void pion_ad_temperature(double **Tion);

/* verner_photo.dat, indexed [Z][electrons]. Ions the file lacks are zero. */
const struct VERNER_STRUCT (*pion_ad_verner_full(void))[31];

/* verner_partial_PIsigmas.dat, indexed [Z][record], with the number of
 * records for each Z. Loops run to the count, not to the array size. */
const struct VERNER_PARTIAL_STRUCT (*pion_ad_verner_partial(void))[125];
const int *pion_ad_verner_partial_count(void);

/* line.dat. The rows in file order (for the low-n opacity), and the H- and
 * He-like tables they fill, indexed [0..28] by Z (for recombination lines). */
const struct PION_LINE_ROW *pion_ad_line_rows(int *nrows);
const struct HYDROGEN_STRUCT *pion_ad_hydrogen(void);
const struct HELIUM_STRUCT *pion_ad_helium(void);

/* highn.dat, indexed [Z][electrons], Z = 0..28. The file stops at Z = 26, so
 * Ni (28) reads zeros. */
const struct HIGHER_ORDER_STRUCT (*pion_ad_highn(void))[3];

/* Element symbol used in the data file names (C N O Ne Mg Al Si S Ar Ca Fe
 * Ni), or NULL for an element the package has no FAC data for. */
const char *pion_ad_symbol(int Z);

/* FAC data for one ion: shell PION_L_SHELL or PION_M_SHELL, element Z,
 * nelec electrons. The photoionization records of <El><nn>a.pi_short and
 * the transition rows of <El><nn>a.tr_short, in file order. A missing file
 * reports, raises the flag and gives 0 rows. */
const struct PION_FAC_PIREC *pion_ad_fac_pi(int shell, int Z, int nelec, int *n);
const struct PION_FAC_TRROW *pion_ad_fac_tr(int shell, int Z, int nelec, int *n);

/* Emission-model data, L shell only. trates: the 10 rows for nelec (3..10).
 * rr_short: records in the pi_short layout. .dat: rows in file order, a
 * multiple of 10. */
const struct PION_TRATES_ROW *pion_ad_trates(int Z, int nelec);
const struct PION_TRSHORTER_ROW *pion_ad_tr_shorter(int Z, int nelec, int *n);
const struct PION_RRLINE_ROW *pion_ad_rrlines(int Z, int nelec, int *n);
const struct PION_FAC_PIREC *pion_ad_rr_short(int Z, int nelec, int *n);

/* H_recombination.dat and He_recombination.dat, indexed [0..28] by Z. */
const struct H_REC_STRUCT *pion_ad_h_rec(void);
const struct HE_REC_STRUCT *pion_ad_he_rec(void);

/* xi_ions.dat: ion fractions on a grid of xi, for the 12 elements xstar
 * computes (H He C N O Ne Mg Si S Ar Ca Fe, in that file order). For element
 * j (1..12, of Z = PION_XI_Z[j]): xi[j][i] and frac[j][(k-1)*nxi + i-1], for
 * i = 1..nxi and k = 1..Z+1, raw as the file gives them. */
struct PION_XI_IONS {
  int nxi;
  double *xi[13];
  double *frac[13];
};
extern const int PION_XI_Z[13];
const struct PION_XI_IONS *pion_ad_xi_ions(void);

/* neutral.tau, one entry per line, as `sscanf(line,"%lf%lf%lf",...)` left it:
 * n fields converted (0 for the header lines), then the values. */
struct PION_NTAU_ROW { int n; double v[3]; };
const struct PION_NTAU_ROW *pion_ad_neutral_tau(int *n);

#endif
