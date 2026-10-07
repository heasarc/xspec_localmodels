#include "photoion_rebin.h"

/* The arithmetic below is that of the loop it replaces, operation for
   operation, so that results are unchanged bit for bit. */

static double rebin_width(const struct pion_rebin *rb, int j)
{
  return rb->E_bin ? rb->E_bin[j] : rb->EBIN;
}

static void rebin_set_ear(struct pion_rebin *rb)
{
  rb->earlo=rb->escale*rb->ear[rb->i];
  rb->earhi=rb->escale*rb->ear[rb->i+1];
  rb->earBIN=rb->earhi-rb->earlo;
}

static void rebin_set_model(struct pion_rebin *rb)
{
  rb->Elo=rb->E_array[rb->j]-rebin_width(rb,rb->j)/2.;
  rb->Ehi=rb->E_array[rb->j]+rebin_width(rb,rb->j)/2.;
}

void pion_rebin_begin(struct pion_rebin *rb, const float ear[], int ne,
                      double escale, const double E_array[],
                      const double E_bin[], double EBIN, int nbins)
{
  rb->ear=ear;
  rb->ne=ne;
  rb->escale=escale;
  rb->E_array=E_array;
  rb->E_bin=E_bin;
  rb->EBIN=EBIN;
  rb->nbins=nbins;
  rb->j=1;
  rb->i=0;
  rebin_set_ear(rb);
  rebin_set_model(rb);
  /* skip the XSPEC bins below the model grid, then the model bins below
     the XSPEC grid */
  while (rb->earhi<=rb->Elo && rb->i<rb->ne) {
    ++rb->i;
    if (rb->i<rb->ne) rebin_set_ear(rb);
  }
  while (rb->Ehi<rb->earlo && rb->j<=rb->nbins) {
    ++rb->j;
    if (rb->j<=rb->nbins) rebin_set_model(rb);
  }
}

int pion_rebin_next(struct pion_rebin *rb, int *i, int *j,
                    double *Ewidth, double *earBIN)
{
  double w;

  if (!(rb->j<=rb->nbins && rb->i<rb->ne)) return 0;
  w=rebin_width(rb,rb->j);
  if (rb->Elo<rb->earlo) {
    if (rb->Ehi-rb->earlo<rb->earBIN) *Ewidth=rb->Ehi-rb->earlo;
    else *Ewidth=rb->earBIN;
  } else {
    if (rb->earhi-rb->Elo<w) *Ewidth=rb->earhi-rb->Elo;
    else *Ewidth=w;
  }
  *i=rb->i;
  *j=rb->j;
  *earBIN=rb->earBIN;
  /* advance whichever bin ends first; the edge loads are guarded so that
     neither ear[ne+1] nor E_array[nbins+1] is read */
  if (rb->earhi<rb->Ehi) {
    ++rb->i;
    if (rb->i<rb->ne) rebin_set_ear(rb);
  } else {
    ++rb->j;
    if (rb->j<=rb->nbins) rebin_set_model(rb);
  }
  return 1;
}
