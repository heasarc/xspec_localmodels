#ifndef PHOTOION_CONST_H
#define PHOTOION_CONST_H

/* Physical constants shared by the photoion models, in cgs units.
 *
 * These were previously copied into every model file. The values here are
 * the ones those copies held, unchanged, so that this consolidation alters no
 * numbers. Most are CODATA 1986 (the thesis dates from c. 2000). Exceptions
 * are noted per line, including two different values of hc, which are kept
 * apart under honest names rather than silently merged.
 *
 * This header is included by every model file, so the names are macros in
 * every translation unit: do not reuse them as local identifiers. */

#define PI (3.141592653589793)

#define ccc (2.99792458e10)          /* speed of light [cm/s], exact */
#define hhh (6.6260755e-27)          /* Planck constant [erg s], CODATA 1986 */
#define eVtoergs (1.60217733e-12)    /* electron volt [erg], CODATA 1986 */
#define ergstoeV (1./eVtoergs)       /* convert ergs to eV */
#define a0 (5.29177249e-9)           /* Bohr radius [cm], CODATA 1986 */
#define re (2.81794092e-13)          /* classical electron radius e^2/(m c^2) [cm], CODATA 1986 */
#define sigmaT (6.6525e-25)          /* Thomson cross-section 8 pi/3 re^2 [cm^2], 1986 rounded to 5 figures */
#define meeV (5.1099906e5)           /* electron rest energy [eV], CODATA 1986 */
#define FINE_STRUCTURE (1./137.0359895)        /* CODATA 1986 */
#define eVtoHartree (1./(2.*13.6056981))       /* 1/(Hartree energy in eV), CODATA 1986 */

/* hc in keV Angstrom, i.e. E[keV] = hc/lambda[Angstrom]. Two values are in
 * use: the absorption models (neutral, vneutral, miabs, siabs, xiabs) use
 * the 1986 value, and the emission models (photoion, phsi, phxi) and the *ext
 * models use a 1998 one, whose last digit differs from the CODATA 1998 table
 * (12.39841857) by one. They differ by 4.7e-7 relative. */
#define HC_KEV_ANGSTROM_1986 (12.3984244)
#define HC_KEV_ANGSTROM_1998 (12.39841856)

/* sqrt(pi) m_e c/(pi e^2) = 1/(sqrt(pi) re c), derived from 1986 values. It
 * appears only in a tau_lim cutoff test on line-centre optical depth. */
#define FACTOR (66.784)

/* Classical radiation damping constant 8 pi^2 e^2/(3 m_e c^3) [s], so that
 * gamma = DAMPING_CLASSICAL * nu^2. Written to 3 figures; the exact 1986
 * value is 2.47389e-22. */
#define DAMPING_CLASSICAL (2.47*pow(10.,-22.))

#define parsectocm (3.085678e18)     /* parsec [cm] */

/* Hubble constant [1/s]: 71 km/s/Mpc (WMAP), 71*1e5/3.085678e18/1e6 */
#define H_0 (2.301e-18)

#endif
