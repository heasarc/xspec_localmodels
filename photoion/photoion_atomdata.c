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
#include "photoion_alloc.h"   /* pion_error */
#include "photoion_const.h"   /* HC_KEV_ANGSTROM */

#define cube(X) ((X)*(X)*(X))

/* Longest suffix appended to PHOTOION_DIR anywhere in the package is under
 * 50 characters ("/photoion_dat/verner_partial_PIsigmas.dat" is 41). */
#define PATH_PAD 64

enum { UNREAD = 0, LOADED, MISSING };

static char *datadir = NULL;     /* the directory the cache belongs to */
static int failed = 0;           /* a file was missing since pion_ad_begin */

static int st_abund, st_oshe, st_temp, st_vfull, st_vpart, st_lines, st_highn;

static double abund[31];
static double oshe_raw[31];
static int    oshe_set[31];
static double oshe_last;
static double temp_val[31][31];
static int    temp_set[31][31];
static struct VERNER_STRUCT vfull[31][31];
static struct VERNER_PARTIAL_STRUCT vpart[29][125];
static int    npart[29];
static struct PION_LINE_ROW *line_rows = NULL;
static int    nline_rows = 0;
static struct HYDROGEN_STRUCT hyd[29];
static struct HELIUM_STRUCT   hel[29];
static struct HIGHER_ORDER_STRUCT hn[29][3];   /* [29]: the emission PE-rate block indexes Ni (28) */

static void reset(void)
{
  st_abund = st_oshe = st_temp = st_vfull = st_vpart = st_lines = st_highn = UNREAD;
  free(line_rows);
  line_rows = NULL;
  nline_rows = 0;
  memset(hyd, 0, sizeof hyd);
  memset(hel, 0, sizeof hel);
  memset(hn, 0, sizeof hn);
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

/* line.dat: "%s%d%d%d%s%lf%s%lf%lf" per row, as the models read it. */
static void load_lines(void)
{
  char s[256];
  int element, electron, LINE, cap = 0;
  double WAVE, AMtemp, fMtemp;
  FILE *linedat;

  if (st_lines != UNREAD) { usable(&st_lines); return; }
  linedat = open_data("line.dat", &st_lines);
  if (linedat == NULL) return;
  while (fscanf(linedat,"%255s%d%d%d%255s%lf%255s%lf%lf",s,&element,&electron,&LINE,s,&WAVE,s,&AMtemp,&fMtemp) != EOF) {
    if (nline_rows == cap) {
      struct PION_LINE_ROW *grown;
      cap = cap ? 2*cap : 256;
      grown = realloc(line_rows, cap * sizeof *line_rows);
      if (grown == NULL) { pion_error("photoion_atomdata: out of memory reading line.dat"); break; }
      line_rows = grown;
    }
    line_rows[nline_rows].element = element;
    line_rows[nline_rows].electron = electron;
    line_rows[nline_rows].LINE = LINE;
    line_rows[nline_rows].WAVE = WAVE;
    line_rows[nline_rows].A = AMtemp;
    line_rows[nline_rows].f = fMtemp;
    ++nline_rows;
    if (element < 0 || element > 28) continue;
    if (electron == 1 && LINE >= 0 && LINE < 8) {
      hyd[element].lambda[LINE]=WAVE;
      hyd[element].A[LINE]=AMtemp;
      hyd[element].f[LINE]=fMtemp;
    } else if (electron == 2 && LINE >= 0 && LINE < 11) {
      hel[element].lambda[LINE]=WAVE;
      hel[element].A[LINE]=AMtemp;
      hel[element].f[LINE]=fMtemp;
    }
  }
  fclose(linedat);
  st_lines = LOADED;
}

const struct PION_LINE_ROW *pion_ad_line_rows(int *nrows)
{
  load_lines();
  *nrows = nline_rows;
  return line_rows;
}

const struct HYDROGEN_STRUCT *pion_ad_hydrogen(void)
{
  load_lines();
  return hyd;
}

const struct HELIUM_STRUCT *pion_ad_helium(void)
{
  load_lines();
  return hel;
}

const struct HIGHER_ORDER_STRUCT (*pion_ad_highn(void))[3]
{
  int element, electron, n;
  double Etemp, ftemp;
  FILE *highnfile;

  if (st_highn == UNREAD) {
    highnfile = open_data("highn.dat", &st_highn);
    if (highnfile) {
      while (fscanf(highnfile,"%d%d%d%lf%lf",&element,&electron,&n,&Etemp,&ftemp)!=EOF) {
        if (element < 0 || element > 28 || electron < 0 || electron > 2 || n < 0 || n > 100) continue;
        hn[element][electron].lambda[n]=HC_KEV_ANGSTROM/Etemp*1000.;
        hn[element][electron].f[n]=ftemp;
      }
      fclose(highnfile);
      st_highn = LOADED;
    }
  }
  usable(&st_highn);
  return (const struct HIGHER_ORDER_STRUCT (*)[3]) hn;
}

