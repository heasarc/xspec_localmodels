/* Atomic data from photoion_dat/, parsed once per process. See the header
 * for the contract.
 *
 * Each dataset is parsed with exactly the fscanf formats and assignments the
 * models used, so the numbers that reach the physics are the same bits. Only
 * the timing changes: once per process instead of once per evaluation.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "photoion_atomdata.h"

#define cube(X) ((X)*(X)*(X))

/* Longest suffix appended to PHOTOION_DIR anywhere in the package is under
 * 50 characters ("/photoion_dat/verner_partial_PIsigmas.dat" is 41). */
#define PATH_PAD 64

enum { UNREAD = 0, LOADED, MISSING };

static char *datadir = NULL;     /* the directory the cache belongs to */
static int failed = 0;           /* a file was missing since pion_ad_begin */

static int st_abund, st_oshe, st_temp, st_vfull, st_vpart;

static double abund[31];
static double oshe_raw[31];
static int    oshe_set[31];
static double oshe_last;
static double temp_val[31][31];
static int    temp_set[31][31];
static struct VERNER_STRUCT vfull[31][31];
static struct VERNER_PARTIAL_STRUCT vpart[29][125];
static int    npart[29];

static void reset(void)
{
  st_abund = st_oshe = st_temp = st_vfull = st_vpart = UNREAD;
  memset(abund, 0, sizeof abund);
  memset(oshe_raw, 0, sizeof oshe_raw);
  memset(oshe_set, 0, sizeof oshe_set);
  oshe_last = 0.;
  memset(temp_val, 0, sizeof temp_val);
  memset(temp_set, 0, sizeof temp_set);
  memset(vfull, 0, sizeof vfull);
  memset(vpart, 0, sizeof vpart);
  memset(npart, 0, sizeof npart);
}

void pion_ad_begin(const char *dir)
{
  if (dir == NULL) dir = "";
  if (datadir == NULL || strcmp(datadir, dir) != 0) {
    size_t n = strlen(dir) + 1;
    free(datadir);
    datadir = malloc(n);
    if (datadir) memcpy(datadir, dir, n);
    reset();
  }
  failed = 0;
}

int pion_ad_failed(void)
{
  return failed;
}

int pion_ad_pathlen(void)
{
  return (int) (datadir ? strlen(datadir) : 0) + PATH_PAD;
}

/* Open photoion_dat/<name> under the current directory. On failure, report
 * once per directory, mark the dataset MISSING and raise the flag. */
static FILE *open_data(const char *name, int *state)
{
  char path[4096];
  FILE *f;

  snprintf(path, sizeof path, "%s/photoion_dat/%s", datadir ? datadir : "", name);
  f = fopen(path, "r");
  if (f == NULL) {
    printf("PHOTOION: cannot open %s\n", path);
    *state = MISSING;
    failed = 1;
  }
  return f;
}

/* A dataset already found missing in this directory raises the flag on every
 * use, so every evaluation that needs it fails, but it is reported once. */
static int usable(int *state)
{
  if (*state == MISSING) { failed = 1; return 0; }
  return 1;
}

void pion_ad_abundance(double ABUND[])
{
  int i, element;
  double djunk;

  if (st_abund == UNREAD) {
    FILE *input = open_data("abundance.dat", &st_abund);
    if (input) {
      while (fscanf(input,"%d%lf",&element,&djunk)!=EOF)
        if (element >= 0 && element <= 30) abund[element]=djunk;
      fclose(input);
      st_abund = LOADED;
    }
  }
  usable(&st_abund);
  for (i=1;i<=30;++i) ABUND[i]=abund[i];
}

static void load_oshe(void)
{
  int element;
  double djunk;

  if (st_oshe == UNREAD) {
    FILE *input = open_data("oscillator_he.dat", &st_oshe);
    if (input) {
      while (fscanf(input,"%d%lf",&element,&djunk)!=EOF) {
        if (element >= 0 && element <= 30) { oshe_raw[element]=djunk; oshe_set[element]=1; }
        oshe_last=djunk;
      }
      fclose(input);
      st_oshe = LOADED;
    }
  }
  usable(&st_oshe);
}

void pion_ad_oscillator_he(double oshe[])
{
  int i;

  load_oshe();
  for (i=0;i<=30;++i) if (oshe_set[i]) oshe[i]=cube(10.)*oshe_raw[i];
}

double pion_ad_oscillator_he_lastraw(void)
{
  load_oshe();
  return oshe_last;
}

void pion_ad_temperature(double **Tion)
{
  int i, j, element, electron;
  double djunk;

  if (st_temp == UNREAD) {
    FILE *input = open_data("temperature.dat", &st_temp);
    if (input) {
      while (fscanf(input,"%d%d%lf",&element,&electron,&djunk)!=EOF)
        if (element >= 0 && element <= 30 && electron >= 0 && electron <= 30) {
          temp_val[element][electron]=djunk;
          temp_set[element][electron]=1;
        }
      fclose(input);
      st_temp = LOADED;
    }
  }
  usable(&st_temp);
  /* Tion is pion_dmatrix(1,28,1,28) */
  for (i=1;i<=28;++i) for (j=1;j<=28;++j) if (temp_set[i][j]) Tion[i][j]=temp_val[i][j];
}

const struct VERNER_STRUCT (*pion_ad_verner_full(void))[31]
{
  int element, electron;
  double Eth,Emax,Ezero,s0,ya,P,yw,y0,y1;

  if (st_vfull == UNREAD) {
    FILE *vernerphoto = open_data("verner_photo.dat", &st_vfull);
    if (vernerphoto) {
      while (fscanf(vernerphoto,"%d%d%lf%lf%lf%lf%lf%lf%lf%lf%lf",&element,&electron,&Eth,&Emax,&Ezero,&s0,&ya,&P,&yw,&y0,&y1)!=EOF) {
        if (element < 0 || element > 30 || electron < 0 || electron > 30) continue;
        vfull[element][electron].Eth=Eth;
        vfull[element][electron].Emax=Emax;
        vfull[element][electron].Ezero=Ezero;
        vfull[element][electron].s0=s0;
        vfull[element][electron].ya=ya;
        vfull[element][electron].P=P;
        vfull[element][electron].yw=yw;
        vfull[element][electron].y0=y0;
        vfull[element][electron].y1=y1;
      }
      fclose(vernerphoto);
      st_vfull = LOADED;
    }
  }
  usable(&st_vfull);
  return (const struct VERNER_STRUCT (*)[31]) vfull;
}

static void load_vpart(void)
{
  int element, electron, principal, angular, element_prev=0, k=0;
  double Eth,Ezero,s0,ya,P,yw;

  if (st_vpart != UNREAD) { usable(&st_vpart); return; }
  {
    FILE *vernerpartial = open_data("verner_partial_PIsigmas.dat", &st_vpart);
    if (vernerpartial == NULL) return;
    /* Records for one element are contiguous. k restarts at each new element,
     * as the models' loop did, and npart records how far it got. */
    while (fscanf(vernerpartial,"%d%d%d%d%lf%lf%lf%lf%lf%lf",&element,&electron,&principal,&angular,&Eth,&Ezero,&s0,&ya,&P,&yw)!=EOF) {
      if (element!=element_prev) k=0;
      if (element >= 0 && element <= 28 && k < 125) {
        vpart[element][k].electron=electron;
        vpart[element][k].principal=principal;
        vpart[element][k].angular=angular;
        vpart[element][k].Eth=Eth;
        vpart[element][k].Ezero=Ezero;
        vpart[element][k].s0=s0;
        vpart[element][k].ya=ya;
        vpart[element][k].P=P;
        vpart[element][k].yw=yw;
        if (k+1 > npart[element]) npart[element]=k+1;
      }
      element_prev=element;
      k=k+1;
    }
    fclose(vernerpartial);
    st_vpart = LOADED;
  }
}

const struct VERNER_PARTIAL_STRUCT (*pion_ad_verner_partial(void))[125]
{
  load_vpart();
  return (const struct VERNER_PARTIAL_STRUCT (*)[125]) vpart;
}

const int *pion_ad_verner_partial_count(void)
{
  load_vpart();
  return npart;
}
