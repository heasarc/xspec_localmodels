/* Cache of the per-edge photoionization tables.
 *
 * Every FAC photoionization record that passes its tau_lim gate needs a
 * 20000-point cross-section table and its spline. Building one is most of
 * an evaluation's cost (50-75% in xiabs, phxi and miabs), and a table depends
 * only on its record, so it is the same in every evaluation, model and
 * block. This keeps them, keyed by record, up to a memory cap:
 *
 *   xset PHOTOION_CACHE_MB <n>    default 1000, minimum 10
 *
 * read at the entry of every model. When a new table would exceed the cap,
 * the least recently used tables are evicted. Eviction only means a later
 * rebuild, so the cap never changes a result. The table in use is never
 * evicted: callers use a table right after getting it, before asking for
 * another. Process lifetime, deliberately not freed (like the parsed data in
 * photoion_atomdata), except on eviction or a PHOTOION_DIR change.
 */
#ifndef PHOTOION_XSTAB_H
#define PHOTOION_XSTAB_H

#include "photoion_cspline.h"
#include "photoion_phys.h"   /* struct PION_FAC_PIREC */

/* Model entry: read PHOTOION_CACHE_MB and evict down to it if it shrank. */
void pion_xstab_begin(void);

/* Drop every table: the records they belong to are about to be freed. */
void pion_xstab_flush(void);

/* The cached table of record r, now the most recently used, or NULL. */
struct pion_cspline *pion_xstab_lookup(const struct PION_FAC_PIREC *r);

/* An empty table slot for record r, made room for under the cap. The caller
 * fills it with pion_cspline_init. NULL only if memory runs out. */
struct pion_cspline *pion_xstab_insert(const struct PION_FAC_PIREC *r, int n);

/* Drop record r's table, e.g. one whose spline failed to build. */
void pion_xstab_forget(const struct PION_FAC_PIREC *r);

/* Current state, for tests and verbose output. */
int pion_xstab_count(void);
double pion_xstab_mbytes(void);
double pion_xstab_cap_mbytes(void);

#endif
