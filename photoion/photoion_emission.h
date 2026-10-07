/* The reemission half of the emission models (photoion, phsi, phxi): what
 * they compute when type > 1. Line emission after photoexcitation, the
 * L-shell recombination lines and continua, and the H- and He-like
 * photoionization rates and recombination spectra. The block was identical
 * in all three models, 426 lines each; it is one function now. Also the
 * tabulated input continuum and the verbose end-of-evaluation report, which
 * were likewise identical in the three. */
#ifndef PHOTOION_EMISSION_H
#define PHOTOION_EMISSION_H

/* Nion, N_e and sigmav_rad are the model's values (the shared-code rule: never
 * read those globals). list[1..ELEMENTS] are the H/He-like elements; Labsorb
 * the absorbed continuum. Results go to the globals exc_spectrum, rec_spectrum,
 * EM and EMion, and to H_excite/He_excite, ratePI/rateRR/rateDR and
 * specRR/specDR, which the model owns. */
void pion_emission_reemission(double **Nion, double N_e, double sigmav_rad,
                              int lines, int verbose,
                              const int list[], int ELEMENTS, int HIGHN,
                              const double oshe[], const double ABUND[], const double Labsorb[],
                              double **H_excite, double **He_excite,
                              double **ratePI, double **rateRR, double **rateDR,
                              double specRR[], double specDR[]);

/* INPUT > 0: the tabulated continuum from input.qdp. 1 if the file is missing. */
int pion_input_continuum(double redshift, int type, int verbose, double Labsorb[]);

/* The end-of-evaluation report (ion table, radiation pressure, power-law
 * norm). Call only when verbose is set. */
void pion_emission_report(int type, double redshift,
                          double **ratePI, const double ABUND[]);

#endif
