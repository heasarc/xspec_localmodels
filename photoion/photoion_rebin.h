#ifndef PHOTOION_REBIN_H
#define PHOTOION_REBIN_H

/* Rebinning of a model's energy grid onto XSPEC's ear grid.

   The model grid has nbins bins, indexed 1..nbins, centred on E_array[j]
   (eV). Their widths are E_bin[j], or EBIN for every bin when E_bin is
   NULL. The XSPEC grid has ne bins with edges ear[0..ne] (keV), mapped to
   the model's frame as escale*ear[i]. Both grids must be in increasing
   energy.

   pion_rebin_next yields each overlap in turn: XSPEC bin i, model bin j,
   the width of their overlap Ewidth and the width of XSPEC bin earBIN,
   both in eV in the model's frame. The caller adds its own quantity to
   photar[i]:

     struct pion_rebin rb;
     pion_rebin_begin(&rb,ear,ne,1000.*(1.+redshift),E_array,NULL,EBIN,SPECBINS);
     while (pion_rebin_next(&rb,&i,&j,&Ewidth,&earBIN))
       photar[i]+=Ewidth/earBIN*exp(-tau[j]);

   Reads no global state, and needs nothing beyond the arguments. */

struct pion_rebin {
  /* private: set by pion_rebin_begin, advanced by pion_rebin_next */
  const float *ear;
  int ne;
  double escale;
  const double *E_array;
  const double *E_bin;
  double EBIN;
  int nbins;
  int i,j;
  double earlo,earhi,earBIN,Elo,Ehi;
};

void pion_rebin_begin(struct pion_rebin *rb, const float ear[], int ne,
                      double escale, const double E_array[],
                      const double E_bin[], double EBIN, int nbins);
int pion_rebin_next(struct pion_rebin *rb, int *i, int *j,
                    double *Ewidth, double *earBIN);

/* A grid read from a file may come in decreasing energy (a wavelength file
   in increasing wavelength). If E_array[1] > E_array[n], reverse E_array,
   E_bin and E_spectrum ([1..n] each) in place, so that the grid increases
   as pion_rebin_begin requires. A grid that already increases is left
   untouched. An unsorted grid is not handled: it is not a valid input. */
void pion_rebin_ascending(double E_array[], double E_bin[],
                          double E_spectrum[], int n);

#endif
