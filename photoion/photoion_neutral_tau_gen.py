#!/usr/bin/env python3
"""Generate photoion_dat/neutral.tau from the vneutral model.

    heainit
    python3 photoion_neutral_tau_gen.py BUILD_DIR [-o OUTPUT]

BUILD_DIR holds a built photoion package (initpackage + hmake), loaded into
XSPEC with "lmod photoion BUILD_DIR". The data are read from the
photoion_dat directory beside this script. OUTPUT defaults to
photoion_dat/neutral.tau there.

neutral.tau is the opacity per H atom [cm^2] of a neutral medium with the
abundances in abundance.dat, on the neutral model's fixed grid of 100 000
bins from 1 eV to 15 keV. One row per bin, no header:

    bin-centre energy [keV]   half bin width [keV]   opacity [cm^2]

The opacity is the sum of two parts:

  * photoabsorption: vneutral at its default parameters (all abundance
    factors 1, sigma_v = 0, so no line absorption), on the same grid
    extended by one bin at the top, both in the model and in XSPEC's
    response, and that extra bin then dropped. Older builds of the models
    never filled the last bin of XSPEC's grid (the original table ends in
    a zero for that reason), so with the extension the last row is right
    to the table's accuracy whichever build generates it;
  * scattering by every electron, bound or free, n_e * sigma_T, where n_e
    is the number of electrons per H atom, sum of Z * abundance over the
    elements in abundance.dat, and sigma_T comes from the CODATA header
    that photoion_const.h includes.

vneutral returns a transmission T = exp(-N_H * sigma), single precision. It
is run at N_H = 0 and at every decade from 1e16 to 1e28 cm^-2, and each bin
takes sigma = -ln(T / T0) / N_H from the largest column whose optical depth
stays below 10. Two effects set that choice:

  * vneutral cuts each edge's high-energy tail where that edge's own optical
    depth falls below tau_lim = 1e-5, so its opacity per H atom grows with
    N_H (at 3 keV: 9.44e-24 cm^2 at tau ~ 0.01, 9.742e-24 at tau ~ 10). The
    table wants the limit with nothing cut off, so the column should be as
    large as possible. At tau ~ 10 the opacity still changes by about 1e-3
    per decade of N_H (at most 4e-3), which is the table's accuracy.
  * XSPEC's single-precision bin edges mix ~1e-7 of each neighbouring bin
    into a bin. Where the neighbour is far more transparent, at an edge,
    that sliver dominates once T falls below ~1e-6 and caps the measured
    optical depth (by 38% at the C I edge at tau ~ 30). At tau <= 10 it is
    at most a few 1e-4. Dividing by the
N_H = 0 run cancels the bin width and anything vneutral adds independently
of N_H. A bin whose optical depth stays below 1e-6 even at N_H = 1e28 is
taken as not absorbing: such a dip is one unit in the last place of a
single-precision number, from the ~1e-7 of a neighbouring bin that XSPEC's
single-precision bin edges mix in.

The script refuses to write a table if the energy grid differs from the
neutral model's, if any opacity is negative, or if vneutral shows any
absorption below 5.5 eV, under the lowest ionization threshold of its
elements (Al I, 5.99 eV).
"""
import argparse
import math
import os
import re
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

# The neutral model's grid (neutral.c: EMIN, EMAX, SPECBINS), in keV.
EMIN, EMAX, NBINS = 0.001, 15.0, 100000

# vneutral's abundance parameters, in the order of photoion_lmodel.dat.
VNEUTRAL_ELEMENTS = 13   # He C N O Ne Mg Al Si S Ar Ca Fe Ni

COLUMNS = [0.0] + [10.0 ** p for p in range(16, 29)]
TAU_LO, TAU_HI = 0.02, 10.0
TAU_FLOOR = 1.e-6   # below this, a single-precision transmission is noise
NO_ABSORPTION_BELOW_KEV = 0.0055


def thomson_cm2():
    """sigma_T [cm^2] from the CODATA header photoion_const.h includes."""
    const = open(os.path.join(HERE, 'photoion_const.h')).read()
    m = re.search(r'#include\s+"(photoion_codata\d+\.h)"', const)
    if not m:
        sys.exit('photoion_const.h: no #include of a photoion_codataYYYY.h')
    header = open(os.path.join(HERE, m.group(1))).read()
    m = re.search(r'#define\s+CODATA_THOMSON_CROSS_SECTION\s+\(([^)]+)\)', header)
    if not m:
        sys.exit('no CODATA_THOMSON_CROSS_SECTION in the CODATA header')
    return float(m.group(1)) * 1.e4   # m^2 -> cm^2


def electrons_per_h():
    """Sum of Z * abundance over abundance.dat."""
    n_e = 0.0
    for line in open(os.path.join(HERE, 'photoion_dat', 'abundance.dat')):
        f = line.split()
        if len(f) >= 2:
            n_e += int(f[0]) * float(f[1])
    return n_e


def run_vneutral(build_dir, workdir):
    """Transmission x bin width for every column; {N_H: [values]}."""
    emax_model = EMIN + (NBINS + 1) * (EMAX - EMIN) / NBINS   # one more bin
    params = ' & '.join(['0'] + ['1'] * VNEUTRAL_ELEMENTS +
                        ['0', '0', repr(EMIN), repr(emax_model), str(NBINS + 1)])
    xcm = ['query yes',
           'xset PHOTOION_DIR %s/' % HERE,
           'lmod photoion %s' % build_dir,
           'dummyrsp %r %r %d lin' % (EMIN, emax_model, NBINS + 1)]
    for i, nh in enumerate(COLUMNS):
        xcm += ['model vneutral*powerlaw & %r & %s & 0 & 1' % (nh, params),
                'set fp [open "col%d.txt" w]' % i,
                'tclout modval',
                'foreach v $xspec_tclout { puts $fp [format "%.9g" $v] }',
                'close $fp']
    xcm.append('exit')
    with open(os.path.join(workdir, 'gen.xcm'), 'w') as f:
        f.write('\n'.join(xcm) + '\n')
    with open(os.path.join(workdir, 'gen.log'), 'w') as log:
        subprocess.run(['xspec', '-', 'gen.xcm'], cwd=workdir, stdout=log,
                       stderr=subprocess.STDOUT, check=False)
    out = {}
    for i, nh in enumerate(COLUMNS):
        path = os.path.join(workdir, 'col%d.txt' % i)
        if not os.path.exists(path):
            sys.exit('XSPEC wrote no output for N_H = %g; see %s'
                     % (nh, os.path.join(workdir, 'gen.log')))
        vals = [float(x) for x in open(path)]
        if len(vals) != NBINS + 1:
            sys.exit('N_H = %g: %d bins, expected %d' % (nh, len(vals), NBINS + 1))
        out[nh] = vals[:NBINS]   # drop the extra bin
    return out


def photoabsorption(cols):
    base = cols[0.0]
    sigma = []
    for k in range(NBINS):
        if base[k] <= 0:
            sys.exit('bin %d: zero flux at N_H = 0' % (k + 1))
        best = None
        for nh in COLUMNS[1:]:
            t = cols[nh][k] / base[k]
            if t <= 0:
                continue
            tau = -math.log(t)
            if TAU_LO < tau < TAU_HI:
                best = (tau, nh)   # COLUMNS ascend: the last one wins
        if best is None:
            t = cols[COLUMNS[-1]][k] / base[k]
            if t > 0 and -math.log(t) < TAU_FLOOR:
                sigma.append(0.0)   # no absorption even at the largest column
                continue
            sys.exit('bin %d: no column gives an optical depth in (%g, %g)'
                     % (k + 1, TAU_LO, TAU_HI))
        sigma.append(best[0] / best[1])
    return sigma


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('build_dir', help='directory of the built photoion package')
    ap.add_argument('-o', '--output',
                    default=os.path.join(HERE, 'photoion_dat', 'neutral.tau'))
    ap.add_argument('--keep', action='store_true',
                    help='keep the XSPEC working directory')
    args = ap.parse_args()

    sigma_t = thomson_cm2()
    n_e = electrons_per_h()
    workdir = tempfile.mkdtemp(prefix='neutral_tau_')
    cols = run_vneutral(os.path.abspath(args.build_dir), workdir)
    sigma = photoabsorption(cols)

    ebin = (EMAX - EMIN) / NBINS   # as neutral.c computes it, in keV
    rows = []
    for k in range(NBINS):
        e = (k + 0.5) * ebin + EMIN
        if sigma[k] < 0:
            sys.exit('bin %d (%g keV): negative opacity' % (k + 1, e))
        if e < NO_ABSORPTION_BELOW_KEV and sigma[k] != 0:
            sys.exit('bin %d (%g keV): vneutral absorbs below every threshold'
                     % (k + 1, e))
        rows.append('%.9g %.9g %.9e' % (e, ebin / 2, sigma[k] + n_e * sigma_t))

    with open(args.output, 'w') as f:
        f.write('\n'.join(rows) + '\n')
    print('wrote %s: %d rows, n_e = %.6f per H, sigma_T = %.10e cm^2'
          % (args.output, len(rows), n_e, sigma_t))
    if args.keep:
        print('XSPEC files kept in', workdir)
    else:
        for name in os.listdir(workdir):
            os.remove(os.path.join(workdir, name))
        os.rmdir(workdir)


if __name__ == '__main__':
    main()
