#ifndef PHOTOION_CONST_H
#define PHOTOION_CONST_H

/* Physical constants shared by the photoion models, in cgs units.
 *
 * Everything here is derived from one CODATA edition, generated from NIST's
 * own table by photoion_codata_gen.py into photoion_codataYYYY.h. To change
 * edition, generate the new header and change this one #include; every value
 * below follows. The exceptions are pi (from GSL) and the parsec (IAU).
 *
 * The models originally carried CODATA 1986 values, with a 1998 value of hc
 * in some of them. Moving everything to one edition changed the spectra by at
 * most 7e-7 in summed flux, and by up to 3e-4 in single bins on steep line
 * wings and edges, almost all of it from unifying hc.
 *
 * This header is included by every model file, so the names are macros in
 * every translation unit: do not reuse them as local identifiers. */

#include "photoion_codata2022.h"
#include <gsl/gsl_math.h>

/* pi and sqrt(pi) from GSL, which defines both correctly rounded in every
 * compiler mode. <math.h> is not enough: M_PI is POSIX, not ISO C, and glibc
 * hides it under a strict -std=c11/c17; and it has no sqrt(pi) constant.
 * sqrt(M_PI) would not do either: it is one ulp away from M_SQRTPI, being the
 * root of an already-rounded pi. */
#define PI (M_PI)
#define SQRT_PI (M_SQRTPI)

#define ccc (CODATA_SPEED_OF_LIGHT*1.e2)                /* speed of light [cm/s] */
#define hhh (CODATA_PLANCK*1.e7)                        /* Planck constant [erg s] */
#define eVtoergs (CODATA_ELECTRON_VOLT*1.e7)            /* electron volt [erg] */
#define ergstoeV (1./eVtoergs)                          /* convert ergs to eV */
#define a0 (CODATA_BOHR_RADIUS*1.e2)                    /* Bohr radius [cm] */
#define re (CODATA_CLASSICAL_ELECTRON_RADIUS*1.e2)      /* classical electron radius [cm] */
#define sigmaT (CODATA_THOMSON_CROSS_SECTION*1.e4)      /* Thomson cross section [cm^2] */
#define meeV (CODATA_ELECTRON_MASS_ENERGY_MEV*1.e6)     /* electron rest energy [eV] */
#define FINE_STRUCTURE (CODATA_FINE_STRUCTURE)
#define eVtoHartree (1./CODATA_HARTREE_ENERGY_EV)       /* 1/(Hartree energy in eV) */

/* hc in keV Angstrom, i.e. E[keV] = HC_KEV_ANGSTROM/lambda[Angstrom]. Since
 * the 2019 SI, h, c and e are all exact, so this is exact too. */
#define HC_KEV_ANGSTROM (hhh*ccc/eVtoergs*1.e5)

/* sqrt(pi) m_e c/(pi e^2) = 1/(sqrt(pi) re c) [s/cm^2]. Used only in a
 * tau_lim cutoff test on line-centre optical depth. */
#define FACTOR (1./(SQRT_PI*re*ccc))

/* Classical radiation damping constant 8 pi^2 e^2/(3 m_e c^3)
 * = 8 pi^2 re/(3 c) [s], so that gamma = DAMPING_CLASSICAL * nu^2. */
#define DAMPING_CLASSICAL (8.*PI*PI*re/(3.*ccc))

/* Parsec [cm]: exactly 648000/pi au (IAU 2015 Resolution B2), with
 * 1 au = 149597870700 m exactly (IAU 2012 Resolution B2). */
#define parsectocm (648000./PI*1.495978707e13)

/* There is deliberately no Hubble constant here. H0 and lambda0 are the
 * user's choice, read from XSPEC's cosmo setting by pion_cosmology(). */

#endif
