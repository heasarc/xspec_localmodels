# Atomic data used by the photoion models

What each file in this directory contains, where its values come from,
and the known errors and gaps in the data themselves. The models are
described in Kinkhabwala et al. (2003, astro-ph/0304332; "K03" below)
and Kinkhabwala et al. (2002, ApJ 575, 732).

The **Source** entries say how well the origin is established:

- **documented**: a paper says so;
- **verified**: checked number by number against the published table;
- **inferred**: consistent with the data, but not documented anywhere we
  know of.

Corrections and additions made since the original release are listed at the
end.

## Photoionization cross sections

**`verner_photo.dat`**: total ground-state photoionization fits, one row
per ion: `Z N E_th E_max E_0 σ_0
y_a P y_w y_0 y_1`, in eV and Mb.
- Source, verified: Verner, Ferland, Korista & Yakovlev (1996, ApJ 465, 487),
  Table 1. All 185 published rows (every ion of Z = 1–14, 16, 18, 20 and 26)
  are identical to the paper's. That paper covers only the Opacity Project
  elements.
- **Two added rows**, Ni XXVIII (`28 1`) and Ni XXVII (`28 2`), added in
  2026. They come from Verner & Yakovlev (1995, A&AS 109, 125), Table 1
  (1s shell), with the six fit values copied verbatim, plus E_max =
  5×10⁴ eV and y_0 = y_1 = 0. That is exactly how Verner et al. (1996)
  treat H- and He-like ions: they "keep the fits to the cross sections from
  Paper II [Verner & Yakovlev 1995]" (all except He I), and their formula
  reduces to the 1995 one when y_0 = y_1 = 0. Verner's own `phfit2.f`
  returns the same values for Ni.

**`verner_partial_PIsigmas.dat`**: partial photoionization fits for every
subshell of every ion, one row per subshell: `Z N n l E_th E_0 σ_0 y_a P y_w`.
- Source, verified: Verner & Yakovlev (1995, A&AS 109, 125), Table 1. All
  1431 rows are identical to the published table.

**He I** has no table here. Its cross section is the fit of Yan,
Sadeghpour & Dalgarno (1998, ApJ 496, 1044), eq. 14 and Table 4 (verified),
for E ≥ 24.58 eV.

## FAC atomic data (`L_shell/`, `M_shell/`)

Calculated with the Flexible Atomic Code (Gu 2002, 2003) for the elements
C, N, O, Ne, Mg, Al, Si, S, Ar, Ca, Fe and Ni. Source: documented, K03 §4.
File names are `<El><nn>…`, where `nn` is the number of electrons.

| file | contents |
|---|---|
| `<El><nn>a.pi_short` | photoionization: level indices and 2J, threshold [eV], fit parameters p0–p3 for the high-energy tail (FAC `RR_RECORD`), and a 6-point table near threshold: E − E_th [eV], an RR cross section (inferred; the column is not labelled), the PI cross section [10⁻²⁰ cm²], and a fourth column of unknown meaning |
| `<El><nn>a.tr_short` | photoexcitation transitions: levels, 2J, energy [eV], f, A |
| `<El><nn>a.tr_shorter` | a subset of the `tr_short` lines, all with f ≥ 10⁻³ (the selection rule is not documented: e.g. 78 of Fe X's 90 such lines), with extra rate columns |
| `<El><nn>a.rr_short` | radiative-recombination data, same layout as `pi_short` |
| `<El><nn>.dat` | recombination-line and continuum strengths (RR and DR) at 10 temperatures, in groups of 10 rows |
| `trates<El>.dat` | total RR and DR rate coefficients at 10 temperatures, for each L-shell ion |

K03 §4 says recombination includes every level to n ≤ 25 explicitly, n = 45
for interpolation and a hydrogenic sum beyond that, plus two-photon decays
(Drake 1986), which FAC does not calculate.

The meaning of the fourth column of the `pi_short` table is not
documented. In every record sampled, column 3 divided by column 4 is
128.4.

## H- and He-like lines and recombination

**`line.dat`**: wavelengths [Å], A and f for H-like lines 1–6 (Lyα…, then
the continuum edge) and He-like lines 1–9 (f, i, r, …, then the edge), for
C, N, O, F, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe and Ni. Source, documented (K03 §4): FAC, with the important H- and
He-like transitions checked and, where needed, corrected against Verner et
al. (1996), Verner, Verner & Ferland (1996, ADNDT 64, 1) and NIST.
- The continuum-edge wavelengths match the FAC ground-state thresholds in
  `<El>01a`/`<El>02a.pi_short`, wherever both exist: H-like to 2×10⁻⁵,
  He-like to 3×10⁻⁴ (several are given to only 4 figures).
- Known errors and gaps:
  - He-like lines 7–8 (1s6p, 1s7p) are blank (λ, A and f all 0) for Ca, Fe
    and Ni, and line 8 is blank for Ar.

**`highn.dat`**: H- and He-like lines with upper levels n = 6–100:
`Z N n E f`. The energies are documented (K03: photoexcitation to n ≤ 100 for
H- and He-like ions, from FAC). The oscillator strengths are exactly C/n³,
with C = 1.6 for H-like (Bethe & Salpeter) and C = 1000 × the value in
`oscillator_he.dat` for He-like. There are no rows for Ni.

**`oscillator_he.dat`**: the He-like constant C/1000, per element (3.12–3.21
for C to Ni). Source: inferred, the n³ f_n asymptote of the FAC He-like
series; not documented.

**`H_recombination.dat`, `He_recombination.dat`**: for each element, at 11
temperatures, the fractions of recombinations producing each H-like line
(1–5) or He-like line (1–7), the recombination-continuum fraction, and the
total recombination coefficient. Source: documented (K03 §4), the FAC
recombination-cascade calculations described above.

## Abundances, ionization structure, temperatures

**`abundance.dat`**: elemental abundances relative to H. Source, documented
in the model README: Wilms, Allen & McCray (2000, ApJ 542, 914), Table 2,
"ISM" column.

**`xi_ions.dat`**: fractional ionic abundances against the ionization
parameter, for H, He, C, N, O, Ne, Mg, Si, S, Ar, Ca and Fe (no Al or
Ni). Source, documented in the model README: an XSTAR simulation of a
very low-column-density medium irradiated by a Γ = 2 power law. The XSTAR version and its other settings are not recorded.

**`xi.dat`**: an *example* ionization-parameter distribution, not atomic
data. See the main README.

**`temperature.dat`**: an electron temperature [eV] for each ion. Source: not documented, presumably the same XSTAR
run.

## `neutral.tau` and `neutral_original.tau`

**`neutral.tau`**: the opacity per H atom [cm²] of a neutral medium, on a
grid of 100 000 bins from 1 eV to 15 keV, one row per bin, no header:
bin-centre energy [keV], half bin width [keV], opacity. Source, documented:
generated by `photoion_neutral_tau_gen.py` (in the package directory),
whose header describes the method. In short, it is the photoabsorption of
this package's `vneutral` model at its default parameters, plus
scattering by all electrons, n_e σ_T, with n_e = Σ Z × abundance over
`abundance.dat` (1.2038 per H) and σ_T from CODATA 2022. `vneutral` cuts
the high-energy tail of each edge where that edge's optical depth falls
below 10⁻⁵, so its opacity per H grows slightly with column. The table
takes it at an optical depth near 10, which is within about 10⁻³ of the
limit. Regenerate it whenever `vneutral`'s physics or `abundance.dat`
changes.

**`neutral_original.tau`**: the table distributed with the original
models, kept for reference and not read by any model. Same layout. (As
distributed it began with three header lines, `READ SERR 1`,
`@neutral.pco` and `!`, written by XSPEC's `wdata` plot command; they
held no information and were removed.) Its last row is 0. **Source: not
documented; inferred, firmly.** It is the same construction as above,
with three differences:
- He I is the Verner et al. (1996) fit, not Yan et al. (1998), which
  falls 13–23% below Yan et al. at 1–10 keV.
- n_e = 1.2 exactly, and σ_T = 6.6525×10⁻²⁵ cm².
- The edge tails were cut off at a lower column.

The evidence:
- Its energy grid is `vneutral`'s default grid, bin for bin, and its C,
  O, Ne and Fe K edges fall in the same bins as `vneutral`'s.
- Above 14 eV it matches `vneutral` at a low column (optical depth
  0.02–0.2) to 0.2% in the median bin and 1.6% in the worst 0.1%, once the
  Verner He I is substituted. The original 2003 code and today's give
  identical results.
- Below 5 eV, under every ionization threshold, it is 1.2 σ_T to 10⁻⁶.
- It is **not** XSPEC's `tbabs`, current or the 2000 version
  (`TBABSVERSION 1`): those differ by 3–13%, as they include H₂ and dust
  grains. It is close to `phabs` with `abund wilm` and `xsect vern` (0.1%
  at 20–520 eV), which is the same atomic physics.

Compared with it, `neutral.tau` is 1–5% higher at 30 eV–5 keV (mostly He
I) and 0.1–0.9% higher above 7 keV.

## Corrections and additions since the original release

- **2026**: `verner_photo.dat` gains Ni XXVIII and XXVII (above).
- **2026**: `line.dat`'s Ca XX continuum-edge wavelength was a copy of
  S XVI's (3.5483 Å). It is now 2.2667 Å, from the FAC Ca XX threshold
  (5469.8 eV in `Ca01a.pi_short`), the source of the other edges.
- **2026**: the P (Z = 15) rows were removed from `line.dat`. Every one
  was a copy of the S row.
- **2026**: `neutral.tau` regenerated from `vneutral` by
  `photoion_neutral_tau_gen.py`; the original renamed
  `neutral_original.tau`, without its header (above).
- **2026**: `line.dat` repaired: see the photoion PR history (heasarc/xspec_localmodels #2).
