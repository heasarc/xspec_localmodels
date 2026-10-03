/* Cache of the per-edge photoionization tables. See the header. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "photoion_xstab.h"

char* FGMSTR(char* name);   /* XSPEC: the value of an xset variable, or "" */

#define CAP_DEFAULT_MB 1000.
#define CAP_MIN_MB     10.

struct entry {
  struct pion_cspline sp;
  struct PION_FAC_PIREC *rec;    /* owner; rec->cache points back here */
  struct entry *prev, *next;     /* LRU list: head is the most recent */
  size_t bytes;
};

static struct entry *head = NULL, *tail = NULL;
static int count = 0;
static size_t used = 0;
static size_t cap = (size_t) (CAP_DEFAULT_MB * 1048576.);
static int warned_invalid = 0, warned_small = 0;

static void unlink_entry(struct entry *e)
{
  if (e->prev) e->prev->next = e->next; else head = e->next;
  if (e->next) e->next->prev = e->prev; else tail = e->prev;
  e->prev = e->next = NULL;
}

static void push_head(struct entry *e)
{
  e->prev = NULL;
  e->next = head;
  if (head) head->prev = e; else tail = e;
  head = e;
}

static void drop(struct entry *e)
{
  unlink_entry(e);
  if (e->rec) e->rec->cache = NULL;
  pion_cspline_free(&e->sp);
  used -= e->bytes;
  --count;
  free(e);
}

/* Evict from the least recently used end until `need` more bytes fit, never
 * evicting `keep`. */
static void make_room(size_t need, const struct entry *keep)
{
  while (tail && tail != keep && used + need > cap) drop(tail);
}

void pion_xstab_begin(void)
{
  const char *s = FGMSTR("PHOTOION_CACHE_MB");
  double mb = CAP_DEFAULT_MB;
  char *end;

  if (s && *s) {
    mb = strtod(s, &end);
    while (*end == ' ' || *end == '\t') ++end;
    if (end == s || *end != '\0' || !(mb == mb)) {
      if (!warned_invalid) printf("PHOTOION: PHOTOION_CACHE_MB=\"%s\" is not a number; using %g\n", s, CAP_DEFAULT_MB);
      warned_invalid = 1;
      mb = CAP_DEFAULT_MB;
    } else if (mb < CAP_MIN_MB) {
      if (!warned_small) printf("PHOTOION: PHOTOION_CACHE_MB=%g is below the minimum; using %g\n", mb, CAP_MIN_MB);
      warned_small = 1;
      mb = CAP_MIN_MB;
    }
  }
  cap = (size_t) (mb * 1048576.);
  make_room(0, NULL);
}

void pion_xstab_flush(void)
{
  while (head) drop(head);
}

struct pion_cspline *pion_xstab_lookup(const struct PION_FAC_PIREC *r)
{
  struct entry *e = (struct entry *) r->cache;
  if (e == NULL) return NULL;
  unlink_entry(e);
  push_head(e);
  return &e->sp;
}

struct pion_cspline *pion_xstab_insert(const struct PION_FAC_PIREC *r, int n)
{
  struct entry *e;
  size_t bytes = sizeof *e + pion_cspline_bytes(n);

  make_room(bytes, NULL);
  e = calloc(1, sizeof *e);
  if (e == NULL) return NULL;
  e->rec = (struct PION_FAC_PIREC *) r;   /* the record is ours to annotate */
  e->rec->cache = e;
  e->bytes = bytes;
  used += bytes;
  ++count;
  push_head(e);
  return &e->sp;
}

void pion_xstab_forget(const struct PION_FAC_PIREC *r)
{
  if (r->cache) drop((struct entry *) r->cache);
}

int pion_xstab_count(void) { return count; }
double pion_xstab_mbytes(void) { return used / 1048576.; }
double pion_xstab_cap_mbytes(void) { return cap / 1048576.; }
