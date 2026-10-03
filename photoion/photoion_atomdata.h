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

/* A buffer size that holds any data path under the current directory. */
int pion_ad_pathlen(void);

/* abundance.dat into ABUND[1..30], zero where the file has no entry. */
void pion_ad_abundance(double ABUND[]);

/* oscillator_he.dat, scaled by 10^3, into oshe[] at the elements listed.
 * Other entries are left as they are. */
void pion_ad_oscillator_he(double oshe[]);

/* The last raw value in oscillator_he.dat. Needed only to keep neutral's
 * neutral.tau header bins bit-for-bit while open issue 11 stands. */
double pion_ad_oscillator_he_lastraw(void);

/* temperature.dat into Tion[Z][electrons] at the ions listed. */
void pion_ad_temperature(double **Tion);

/* verner_photo.dat, indexed [Z][electrons]. Ions the file lacks are zero. */
const struct VERNER_STRUCT (*pion_ad_verner_full(void))[31];

/* verner_partial_PIsigmas.dat, indexed [Z][record], with the number of
 * records for each Z. Loops run to the count, not to the array size. */
const struct VERNER_PARTIAL_STRUCT (*pion_ad_verner_partial(void))[125];
const int *pion_ad_verner_partial_count(void);

#endif
