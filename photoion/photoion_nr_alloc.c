/* Numerical Recipes style array allocators shared by all photoion models.
 *
 * Extracted verbatim from the per-model copies, which were identical apart
 * from a name suffix -- with one exception: these are photoion.c's variants,
 * which zero the memory they hand back. The other ten model files used
 * allocators that returned uninitialized malloc'd memory, and several of
 * them accumulate onto arrays with += without initializing them first.
 * Sharing the zeroing version fixes that class of bug library-wide.
 *
 * The unit-offset convention of Numerical Recipes is preserved: a vector
 * allocated [nl..nh] is returned offset so that v[nl] is the first element,
 * and the matching free_* routine adds the offset back before free().
 */

#include <stdlib.h>
#include <stdio.h>

#include "photoion_nr_alloc.h"

void pion_nrerror(char error_text[])
{

   /*fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to system...\n");
   exit(1);*/
}

double *pion_vector(int nl, int nh)
{
   double *v;
   int i,sz;

   sz = nh-nl+1;
   v=(double *)malloc((unsigned) sz*sizeof(double));
   if (!v) pion_nrerror("allocation failure in vector()");
   for (i=0; i<sz; ++i) v[i] = .0; 
   return v-nl;
}

int *pion_ivector(int nl, int nh)
{
   int *v, i, sz;

   sz = nh-nl+1;
   v=(int *)malloc((unsigned) sz*sizeof(int));
   if (!v) pion_nrerror("allocation failure in ivector()");
   for (i=0; i<sz; ++i) v[i] = 0; 
   return v-nl;
}

double *pion_dvector(int nl, int nh)
{
   double *v;
   int i,sz;

   sz = nh-nl+1;
   v=(double *)malloc((unsigned) sz*sizeof(double));
   if (!v) pion_nrerror("allocation failure in dvector()");
   for (i=0; i<sz; ++i) v[i] = .0; 
   return v-nl;
}

double **pion_matrix(int nrl, int nrh, int ncl, int nch)
{
   int i, j, sz;
   double **m, *m_i;

   m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
   if (!m) pion_nrerror("allocation failure 1 in matrix()");
   m -= nrl;

   sz = nch-ncl+1;
   for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned) sz*sizeof(double));
      if (!m[i]) pion_nrerror("allocation failure 2 in matrix()");
      m_i = m[i];
      for (j=0; j<sz; ++j) m_i[j] = .0;
      m[i] -= ncl;
   }
   return m;
}

double **pion_dmatrix(int nrl, int nrh, int ncl, int nch)
{
   int i ,j, sz;
   double **m, *m_i;

   m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
   if (!m) pion_nrerror("allocation failure 1 in dmatrix()");
   m -= nrl;

   sz = nch-ncl+1;
   for(i=nrl;i<=nrh;i++) {
      m[i]=(double *) malloc((unsigned) sz*sizeof(double));
      if (!m[i]) pion_nrerror("allocation failure 2 in dmatrix()");
      m_i = m[i];
      for (j=0; j<sz; ++j) m_i[j] = .0;
      m[i] -= ncl;
   }
   return m;
}

int **pion_imatrix(int nrl, int nrh, int ncl, int nch)
{
   int i,**m, *m_i, j, sz;

   m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
   if (!m) pion_nrerror("allocation failure 1 in imatrix()");
   m -= nrl;

   sz = nch-ncl+1;
   for(i=nrl;i<=nrh;i++) {
      m[i]=(int *)malloc((unsigned) sz*sizeof(int));
      if (!m[i]) pion_nrerror("allocation failure 2 in imatrix()");
      m_i = m[i];
      for (j=0; j<sz; ++j) m_i[j] = 0;
      m[i] -= ncl;
   }
   return m;
}

double **pion_submatrix(double **a, int oldrl, int oldrh, int oldcl, int oldch, int newrl, int newcl)
{
   int i,j;
   double **m;

   m=(double **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(double*));
   if (!m) pion_nrerror("allocation failure in submatrix()");
   m -= newrl;

   for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

   return m;
}

void pion_free_vector(double *v, int nl, int nh)
{
   free((char*) (v+nl));
}

void pion_free_ivector(int *v, int nl, int nh)
{
   free((char*) (v+nl));
}

void pion_free_dvector(double *v, int nl, int nh)
{
   free((char*) (v+nl));
}

void pion_free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void pion_free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void pion_free_imatrix(int **m, int nrl, int nrh, int ncl, int nch)
{
   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
   free((char*) (m+nrl));
}

void pion_free_submatrix(double **b, int nrl, int nrh, int ncl, int nch)
{
   free((char*) (b+nrl));
}

double **pion_convert_matrix(double *a, int nrl, int nrh, int ncl, int nch)
{
   int i,j,nrow,ncol;
   double **m;

   nrow=nrh-nrl+1;
   ncol=nch-ncl+1;
   m = (double **) malloc((unsigned) (nrow)*sizeof(double*));
   if (!m) pion_nrerror("allocation failure in convert_matrix()");
   m -= nrl;
   for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
   return m;
}

void pion_free_convert_matrix(double **b, int nrl, int nrh, int ncl, int nch)
{
   free((char*) (b+nrl));
}
