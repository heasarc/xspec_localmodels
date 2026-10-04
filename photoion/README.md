Photoionization models from Ali Kinkhabwala produced as part of a PhD
thesis at Columbia Astrophysics Laboratory. A detailed description is
available in <a href="http://arxiv.org/abs/astro-ph/0304332">astro-ph/0304332</a> and an
example of their use in the analysis of an XMM-Newton observation of
MCG -6-30-15 is given in <a href="mcg.ps">this unpublished paper</a>.
<p>

The photoion_dat directory holds the atomic data files. Its location is
specified within XSPEC using the command "xset PHOTOION_DIR
directory-name", where directory-name is the directory <i>containing</i>
photoion_dat &mdash; not photoion_dat itself. The models append
"/photoion_dat/" to this value internally, so if photoion_dat lives in
/path/to/photoion, then set
<pre>
  xset PHOTOION_DIR /path/to/photoion
</pre>
and the models will read /path/to/photoion/photoion_dat/abundance.dat,
and so on. Pointing PHOTOION_DIR directly at photoion_dat instead
produces errors of the form "NEUTRAL: Failed to open
.../photoion_dat/photoion_dat/abundance.dat".
<p>

What each data file contains, where its values come from, and the known
errors in the data are described in photoion_dat/README.md.
<p>

PHOTOION_DIR may be any length. (Versions before 2026 kept the
constructed filenames in a fixed 130-character buffer, and a path of
about 100 characters or more aborted XSPEC.)
<p>
The abundance, He oscillator-strength, temperature, Verner
photoionization, line (line.dat), high-n line (highn.dat), FAC
L- and M-shell (pi_short, tr_short), and emission-model recombination
(trates, tr_shorter, .dat, rr_short, H_ and He_recombination.dat),
xi_ions.dat and neutral.tau tables -- every data file the models use --
are read once per XSPEC session, as each is first needed, and kept, not
re-read on every evaluation.
<p>
The photoionization cross-section table built for each absorption edge
(a 20000-point spline, about 0.3 MB) is also kept, so later evaluations
skip rebuilding it, which is most of an evaluation's cost in xiabs, phxi
and miabs. The memory these tables may use is set by
<pre>
  xset PHOTOION_CACHE_MB 1000
</pre>
in megabytes: 1000 by default, 10 at least. It is read at every evaluation,
so it can be changed at any time. When the cap is reached, the least
recently used tables are dropped and rebuilt when next needed. The cap
affects only speed and memory, never the results. A typical xiabs model
uses about 180 MB. Separately, the parsed atomic data take at most about
60 MB. If one of them cannot be opened, the model
prints "PHOTOION: cannot open <i>path</i>" once and returns zeros. Restart
XSPEC, or change PHOTOION_DIR, to retry.
<p>


<br>
<br>
<br>
<h3>NEUTRAL:</h3>
Here's an example of a set of model parameters
for NEUTRAL applied to a power-law spectrum. Abundances are from the
"ISM" column of Table 2 in Wilms, Allen, &amp; McCray 2000, ApJ, 542,
 914. (see below for a description of the parameters.)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  neutral[1]( powerlaw[2] )
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   neutral    N_H      cm^-2     1.0000E+20 frozen
    2    2    1   neutral    sigma_v  km/s        0.00     frozen
    3    3    2   powerlaw   PhoIndex            2.000     +/-   0.000
    4    4    2   powerlaw   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
where the parameters are<br>
1: N_H - Neutral hydrogen column density.<br>
2: sigma_v - Radial velocity width (sigma) of absorbing medium.
If sigma_v=0, line absorption is NOT included.
If sigma_v&gt;0, line absorption IS included.
<p>
<br>
<br>
<h3>VNEUTRAL: </h3>
Here's an example of a set of model parameters
for VNEUTRAL applied to a power-law spectrum. Abundances are from the "ISM" column of Table
2 in Wilms, Allen, &amp; McCray 2000, ApJ, 542, 914. (See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  vneutral[1]( powerlaw[2] )
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   vneutral   N_H      cm^-2     1.0000E+20 frozen
    2    2    1   vneutral   sigma_v  km/s       0.000     frozen
    3    3    1   vneutral   A_He     abund      1.000     frozen
    4    4    1   vneutral   A_C      abund      1.000     frozen
    5    5    1   vneutral   A_N      abund      1.000     frozen
    6    6    1   vneutral   A_O      abund      1.000     frozen
    7    7    1   vneutral   A_Ne     abund      1.000     frozen
    8    8    1   vneutral   A_Mg     abund      1.000     frozen
    9    9    1   vneutral   A_Al     abund      1.000     frozen
   10   10    1   vneutral   A_Si     abund      1.000     frozen
   11   11    1   vneutral   A_S      abund      1.000     frozen
   12   12    1   vneutral   A_Ar     abund      1.000     frozen
   13   13    1   vneutral   A_Ca     abund      1.000     frozen
   14   14    1   vneutral   A_Fe     abund      1.000     frozen
   15   15    1   vneutral   A_Ni     abund      1.000     frozen
   16   16    1   vneutral   redshift            0.000     frozen
   17   17    1   vneutral   v        km/s       0.000     frozen
   18   18    1   vneutral   EMIN     keV       1.0000E-03 frozen
   19   19    1   vneutral   EMAX     keV        15.00     frozen
   20   20    1   vneutral   SPECBINS           1.0000E+05 frozen
   21   21    2   powerlaw   PhoIndex            2.000     +/-   0.000
   22   22    2   powerlaw   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: N_H - Neutral hydrogen column density.<br>
2: sigma_v - Radial velocity width (sigma) of absorbing medium.
If sigma_v=0, line absorption is NOT included.
If sigma_v&gt;0, line absorption IS included.<br>
3: A_He - Overall Helium abundance relative to "solar" default value.<br>
.<br>
.<br>
.<br>
15: A_Ni - Overall Nickel abundance relative to "solar" default value.<br>
16: redshift - Redshift of absorbing medium.<br>
17: v - Radial velocity shift of absorbing medium.<br>
18: EMIN - Minimum energy [keV] for internal grid.<br>
19: EMAX - Maximum energy [keV] for internal grid.<br>
20: SPECBINS - Total number of energy bins (equally-spaced in energy)
for internal grid.<br>
<p>

<br>
<br>
<br>

<h3>SIABS: </h3>
Here's an example of a set of model parameters for SIABS applied to a
power-law spectrum. (See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  siabs[1]( powerlaw[2] )
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   siabs      Z                   8.000     frozen
    2    2    1   siabs      z                   2.000     frozen
    3    3    1   siabs      Nion     cm^-2     1.0000E+17 frozen
    4    4    1   siabs      redshift            0.000     frozen
    5    5    1   siabs      v        km/s       0.000     frozen
    6    6    1   siabs      sigma_v  km/s       100.0     frozen
    7    7    1   siabs      EMIN     keV       1.0000E-03 frozen
    8    8    1   siabs      EMAX     keV        15.00     frozen
    9    9    1   siabs      SPECBINS           1.0000E+05 frozen
   10   10    2   powerlaw   PhoIndex            2.000     +/-   0.000
   11   11    2   powerlaw   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: Z - Atomic number<br>
2: z - Number of electrons<br>
3: Nion - Ion column density [cm^-2]<br>
4: redshift - Redshift of source<br>
5: v - Velocity shift<br>
6: sigma_v - Velocity width (sigma)<br>
7: EMIN - Minimum energy [keV] for internal grid.  <br>
8: EMAX - Maximum energy [keV] for internal grid.  <br>
9: SPECBINS - Total number of energy bins (equally-spaced in energy) for internal grid.<br>
<p>
<br>
<br>
<br>

<h3>XIABS: </h3>
Here's an example of a set of model parameters
for XIABS applied to a power-law spectrum. This models uses a
user-defined distribution in ionization parameter. The file "xi.dat"
must exist in the directory you're running XSPEC in. An example can be
found in photoion_dat/xi.dat. The ionization parameter distribution is 
defined by simply connecting the 
user-defined points in xi.dat with line segments and normalizing. The
"fractional ionic abundances" used were taken from an XSTAR simulation
of an extremely-low-column-density medium irradiated by a Gamma=2 power
law. (See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  xiabs[1]( powerlaw[2] )
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   xiabs      N_H      cm^-2     1.0000E+22 frozen
    2    2    1   xiabs      A_He     abund      1.000     frozen
    3    3    1   xiabs      A_C      abund      1.000     frozen
    4    4    1   xiabs      A_N      abund      1.000     frozen
    5    5    1   xiabs      A_O      abund      1.000     frozen
    6    6    1   xiabs      A_Ne     abund      1.000     frozen
    7    7    1   xiabs      A_Mg     abund      1.000     frozen
    8    8    1   xiabs      A_Al     abund      0.000     frozen
    9    9    1   xiabs      A_Si     abund      1.000     frozen
   10   10    1   xiabs      A_S      abund      1.000     frozen
   11   11    1   xiabs      A_Ar     abund      0.000     frozen
   12   12    1   xiabs      A_Ca     abund      0.000     frozen
   13   13    1   xiabs      A_Fe     abund      1.000     frozen
   14   14    1   xiabs      A_Ni     abund      0.000     frozen
   15   15    1   xiabs      redshift            0.000     frozen
   16   16    1   xiabs      v        km/s       0.000     frozen
   17   17    1   xiabs      sigma_v  km/s       100.0     frozen
   18   18    1   xiabs      EMIN     keV       1.0000E-03 frozen
   19   19    1   xiabs      EMAX     keV        10.00     frozen
   20   20    1   xiabs      SPECBINS           4.0000E+04 frozen
   21   21    1   xiabs      verbose             1.000     frozen
   22   22    2   powerlaw   PhoIndex            2.000     +/-   0.000
   23   23    2   powerlaw   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: N_H - Total hydrogen column density (neutral plus ionized).<br>
2: A_He - Overall Helium abundance relative to "solar" default value.<br>
.<br>
.<br>
.<br>
14: A_Ni - Overall Nickel abundance relative to "solar" default value.<br>
15: redshift - Redshift of source<br>
16: v - Radial velocity shift <br>
17: sigma_v - Radial velocity width (sigma)<br>
18: EMIN - Minimum energy [keV] for internal grid.<br>
19: EMAX - Maximum energy [keV] for internal grid.<br>
20: SPECBINS - Total number of energy bins (equally-spaced in energy) for internal grid.<br>
21: verbose - =1 to output numbers/messages, =0 for no output<br>
<p>
<br>
<br>
<br>

<h3>MIABS: </h3>
Here's an example of a set of model parameters for MIABS applied to a
power-law spectrum. (See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  miabs[1]( powerlaw[2] )
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   miabs      redshift            0.000     frozen
    2    2    1   miabs      v        km/s       0.000     frozen
    3    3    1   miabs      sigma_v  km/s       100.0     frozen
    4    4    1   miabs      EMIN     keV       1.0000E-03 frozen
    5    5    1   miabs      EMAX     keV        10.00     frozen
    6    6    1   miabs      SPECBINS           4.0000E+04 frozen
    7    7    1   miabs      verbose             1.000     frozen
    8    8    1   miabs      LOG?                0.000     frozen
    9    9    1   miabs      COLNORM             1.000     frozen
   10   10    1   miabs      N_e      cm^-2      0.000     frozen
   11   11    1   miabs      H_1      cm^-2      0.000     frozen
   12   12    1   miabs      He_1     cm^-2      0.000     frozen
   13   13    1   miabs      He_2     cm^-2      0.000     frozen
   14   14    1   miabs      C_1      cm^-2      0.000     frozen
   15   15    1   miabs      C_2      cm^-2      0.000     frozen
   16   16    1   miabs      C_3      cm^-2      0.000     frozen
   17   17    1   miabs      C_4      cm^-2      0.000     frozen
   18   18    1   miabs      C_5      cm^-2      0.000     frozen
   19   19    1   miabs      C_6      cm^-2      0.000     frozen
   20   20    1   miabs      N_1      cm^-2      0.000     frozen
   21   21    1   miabs      N_2      cm^-2      0.000     frozen
   22   22    1   miabs      N_3      cm^-2      0.000     frozen
   23   23    1   miabs      N_4      cm^-2      0.000     frozen
   24   24    1   miabs      N_5      cm^-2      0.000     frozen
   25   25    1   miabs      N_6      cm^-2      0.000     frozen
   26   26    1   miabs      N_7      cm^-2      0.000     frozen
   27   27    1   miabs      O_1      cm^-2      0.000     frozen
   28   28    1   miabs      O_2      cm^-2      0.000     frozen
   29   29    1   miabs      O_3      cm^-2      0.000     frozen
   30   30    1   miabs      O_4      cm^-2      0.000     frozen
   31   31    1   miabs      O_5      cm^-2      0.000     frozen
   32   32    1   miabs      O_6      cm^-2      0.000     frozen
   33   33    1   miabs      O_7      cm^-2      0.000     frozen
   34   34    1   miabs      O_8      cm^-2      0.000     frozen
   35   35    1   miabs      Ne_1     cm^-2      0.000     frozen
   36   36    1   miabs      Ne_2     cm^-2      0.000     frozen
   37   37    1   miabs      Ne_3     cm^-2      0.000     frozen
   38   38    1   miabs      Ne_4     cm^-2      0.000     frozen
   39   39    1   miabs      Ne_5     cm^-2      0.000     frozen
   40   40    1   miabs      Ne_6     cm^-2      0.000     frozen
   41   41    1   miabs      Ne_7     cm^-2      0.000     frozen
   42   42    1   miabs      Ne_8     cm^-2      0.000     frozen
   43   43    1   miabs      Ne_9     cm^-2      0.000     frozen
   44   44    1   miabs      Ne_10    cm^-2      0.000     frozen
   45   45    1   miabs      Mg_1     cm^-2      0.000     frozen
   46   46    1   miabs      Mg_2     cm^-2      0.000     frozen
   47   47    1   miabs      Mg_3     cm^-2      0.000     frozen
   48   48    1   miabs      Mg_4     cm^-2      0.000     frozen
   49   49    1   miabs      Mg_5     cm^-2      0.000     frozen
   50   50    1   miabs      Mg_6     cm^-2      0.000     frozen
   51   51    1   miabs      Mg_7     cm^-2      0.000     frozen
   52   52    1   miabs      Mg_8     cm^-2      0.000     frozen
   53   53    1   miabs      Mg_9     cm^-2      0.000     frozen
   54   54    1   miabs      Mg_10    cm^-2      0.000     frozen
   55   55    1   miabs      Mg_11    cm^-2      0.000     frozen
   56   56    1   miabs      Mg_12    cm^-2      0.000     frozen
   57   57    1   miabs      Al_1     cm^-2      0.000     frozen
   58   58    1   miabs      Al_2     cm^-2      0.000     frozen
   59   59    1   miabs      Al_3     cm^-2      0.000     frozen
   60   60    1   miabs      Al_4     cm^-2      0.000     frozen
   61   61    1   miabs      Al_5     cm^-2      0.000     frozen
   62   62    1   miabs      Al_6     cm^-2      0.000     frozen
   63   63    1   miabs      Al_7     cm^-2      0.000     frozen
   64   64    1   miabs      Al_8     cm^-2      0.000     frozen
   65   65    1   miabs      Al_9     cm^-2      0.000     frozen
   66   66    1   miabs      Al_10    cm^-2      0.000     frozen
   67   67    1   miabs      Al_11    cm^-2      0.000     frozen
   68   68    1   miabs      Al_12    cm^-2      0.000     frozen
   69   69    1   miabs      Al_13    cm^-2      0.000     frozen
   70   70    1   miabs      Si_1     cm^-2      0.000     frozen
   71   71    1   miabs      Si_2     cm^-2      0.000     frozen
   72   72    1   miabs      Si_3     cm^-2      0.000     frozen
   73   73    1   miabs      Si_4     cm^-2      0.000     frozen
   74   74    1   miabs      Si_5     cm^-2      0.000     frozen
   75   75    1   miabs      Si_6     cm^-2      0.000     frozen
   76   76    1   miabs      Si_7     cm^-2      0.000     frozen
   77   77    1   miabs      Si_8     cm^-2      0.000     frozen
   78   78    1   miabs      Si_9     cm^-2      0.000     frozen
   79   79    1   miabs      Si_10    cm^-2      0.000     frozen
   80   80    1   miabs      Si_11    cm^-2      0.000     frozen
   81   81    1   miabs      Si_12    cm^-2      0.000     frozen
   82   82    1   miabs      Si_13    cm^-2      0.000     frozen
   83   83    1   miabs      Si_14    cm^-2      0.000     frozen
   84   84    1   miabs      S_1      cm^-2      0.000     frozen
   85   85    1   miabs      S_2      cm^-2      0.000     frozen
   86   86    1   miabs      S_3      cm^-2      0.000     frozen
   87   87    1   miabs      S_4      cm^-2      0.000     frozen
   88   88    1   miabs      S_5      cm^-2      0.000     frozen
   89   89    1   miabs      S_6      cm^-2      0.000     frozen
   90   90    1   miabs      S_7      cm^-2      0.000     frozen
   91   91    1   miabs      S_8      cm^-2      0.000     frozen
   92   92    1   miabs      S_9      cm^-2      0.000     frozen
   93   93    1   miabs      S_10     cm^-2      0.000     frozen
   94   94    1   miabs      S_11     cm^-2      0.000     frozen
   95   95    1   miabs      S_12     cm^-2      0.000     frozen
   96   96    1   miabs      S_13     cm^-2      0.000     frozen
   97   97    1   miabs      S_14     cm^-2      0.000     frozen
   98   98    1   miabs      S_15     cm^-2      0.000     frozen
   99   99    1   miabs      S_16     cm^-2      0.000     frozen
  100  100    1   miabs      Ar_1     cm^-2      0.000     frozen
  101  101    1   miabs      Ar_2     cm^-2      0.000     frozen
  102  102    1   miabs      Ar_3     cm^-2      0.000     frozen
  103  103    1   miabs      Ar_4     cm^-2      0.000     frozen
  104  104    1   miabs      Ar_5     cm^-2      0.000     frozen
  105  105    1   miabs      Ar_6     cm^-2      0.000     frozen
  106  106    1   miabs      Ar_7     cm^-2      0.000     frozen
  107  107    1   miabs      Ar_8     cm^-2      0.000     frozen
  108  108    1   miabs      Ar_9     cm^-2      0.000     frozen
  109  109    1   miabs      Ar_10    cm^-2      0.000     frozen
  110  110    1   miabs      Ar_11    cm^-2      0.000     frozen
  111  111    1   miabs      Ar_12    cm^-2      0.000     frozen
  112  112    1   miabs      Ar_13    cm^-2      0.000     frozen
  113  113    1   miabs      Ar_14    cm^-2      0.000     frozen
  114  114    1   miabs      Ar_15    cm^-2      0.000     frozen
  115  115    1   miabs      Ar_16    cm^-2      0.000     frozen
  116  116    1   miabs      Ar_17    cm^-2      0.000     frozen
  117  117    1   miabs      Ar_18    cm^-2      0.000     frozen
  118  118    1   miabs      Ca_1     cm^-2      0.000     frozen
  119  119    1   miabs      Ca_2     cm^-2      0.000     frozen
  120  120    1   miabs      Ca_3     cm^-2      0.000     frozen
  121  121    1   miabs      Ca_4     cm^-2      0.000     frozen
  122  122    1   miabs      Ca_5     cm^-2      0.000     frozen
  123  123    1   miabs      Ca_6     cm^-2      0.000     frozen
  124  124    1   miabs      Ca_7     cm^-2      0.000     frozen
  125  125    1   miabs      Ca_8     cm^-2      0.000     frozen
  126  126    1   miabs      Ca_9     cm^-2      0.000     frozen
  127  127    1   miabs      Ca_10    cm^-2      0.000     frozen
  128  128    1   miabs      Ca_11    cm^-2      0.000     frozen
  129  129    1   miabs      Ca_12    cm^-2      0.000     frozen
  130  130    1   miabs      Ca_13    cm^-2      0.000     frozen
  131  131    1   miabs      Ca_14    cm^-2      0.000     frozen
  132  132    1   miabs      Ca_15    cm^-2      0.000     frozen
  133  133    1   miabs      Ca_16    cm^-2      0.000     frozen
  134  134    1   miabs      Ca_17    cm^-2      0.000     frozen
  135  135    1   miabs      Ca_18    cm^-2      0.000     frozen
  136  136    1   miabs      Ca_19    cm^-2      0.000     frozen
  137  137    1   miabs      Ca_20    cm^-2      0.000     frozen
  138  138    1   miabs      Fe_1     cm^-2      0.000     frozen
  139  139    1   miabs      Fe_2     cm^-2      0.000     frozen
  140  140    1   miabs      Fe_3     cm^-2      0.000     frozen
  141  141    1   miabs      Fe_4     cm^-2      0.000     frozen
  142  142    1   miabs      Fe_5     cm^-2      0.000     frozen
  143  143    1   miabs      Fe_6     cm^-2      0.000     frozen
  144  144    1   miabs      Fe_7     cm^-2      0.000     frozen
  145  145    1   miabs      Fe_8     cm^-2      0.000     frozen
  146  146    1   miabs      Fe_9     cm^-2      0.000     frozen
  147  147    1   miabs      Fe_10    cm^-2      0.000     frozen
  148  148    1   miabs      Fe_11    cm^-2      0.000     frozen
  149  149    1   miabs      Fe_12    cm^-2      0.000     frozen
  150  150    1   miabs      Fe_13    cm^-2      0.000     frozen
  151  151    1   miabs      Fe_14    cm^-2      0.000     frozen
  152  152    1   miabs      Fe_15    cm^-2      0.000     frozen
  153  153    1   miabs      Fe_16    cm^-2      0.000     frozen
  154  154    1   miabs      Fe_17    cm^-2      0.000     frozen
  155  155    1   miabs      Fe_18    cm^-2      0.000     frozen
  156  156    1   miabs      Fe_19    cm^-2      0.000     frozen
  157  157    1   miabs      Fe_20    cm^-2      0.000     frozen
  158  158    1   miabs      Fe_21    cm^-2      0.000     frozen
  159  159    1   miabs      Fe_22    cm^-2      0.000     frozen
  160  160    1   miabs      Fe_23    cm^-2      0.000     frozen
  161  161    1   miabs      Fe_24    cm^-2      0.000     frozen
  162  162    1   miabs      Fe_25    cm^-2      0.000     frozen
  163  163    1   miabs      Fe_26    cm^-2      0.000     frozen
  164  164    1   miabs      Ni_1     cm^-2      0.000     frozen
  165  165    1   miabs      Ni_2     cm^-2      0.000     frozen
  166  166    1   miabs      Ni_3     cm^-2      0.000     frozen
  167  167    1   miabs      Ni_4     cm^-2      0.000     frozen
  168  168    1   miabs      Ni_5     cm^-2      0.000     frozen
  169  169    1   miabs      Ni_6     cm^-2      0.000     frozen
  170  170    1   miabs      Ni_7     cm^-2      0.000     frozen
  171  171    1   miabs      Ni_8     cm^-2      0.000     frozen
  172  172    1   miabs      Ni_9     cm^-2      0.000     frozen
  173  173    1   miabs      Ni_10    cm^-2      0.000     frozen
  174  174    1   miabs      Ni_11    cm^-2      0.000     frozen
  175  175    1   miabs      Ni_12    cm^-2      0.000     frozen
  176  176    1   miabs      Ni_13    cm^-2      0.000     frozen
  177  177    1   miabs      Ni_14    cm^-2      0.000     frozen
  178  178    1   miabs      Ni_15    cm^-2      0.000     frozen
  179  179    1   miabs      Ni_16    cm^-2      0.000     frozen
  180  180    1   miabs      Ni_17    cm^-2      0.000     frozen
  181  181    1   miabs      Ni_18    cm^-2      0.000     frozen
  182  182    1   miabs      Ni_19    cm^-2      0.000     frozen
  183  183    1   miabs      Ni_20    cm^-2      0.000     frozen
  184  184    1   miabs      Ni_21    cm^-2      0.000     frozen
  185  185    1   miabs      Ni_21    cm^-2      0.000     frozen
  186  186    1   miabs      Ni_23    cm^-2      0.000     frozen
  187  187    1   miabs      Ni_24    cm^-2      0.000     frozen
  188  188    1   miabs      Ni_25    cm^-2      0.000     frozen
  189  189    1   miabs      Ni_26    cm^-2      0.000     frozen
  190  190    1   miabs      Ni_27    cm^-2      0.000     frozen
  191  191    1   miabs      Ni_28    cm^-2      0.000     frozen
  192  192    2   powerlaw   PhoIndex            2.000     +/-   0.000
  193  193    2   powerlaw   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: redshift - Redshift of source<br>
2: v - Radial velocity shift<br>
3: sigma_v - Radial velocity width (sigma)<br>
4: EMIN - Minimum energy [keV] for internal grid.  <br>
5: EMAX - Maximum energy [keV] for internal grid.  <br>
6: SPECBINS - Total number of energy bins (equally-spaced in energy) for internal grid.<br>
7: verbose - =1 to output numbers/messages, =0 for no output<br>
8: LOG? - = 0 to use linear units for column densities, = 1 to use log (base 10) units for column densities<br>
9: COLNORM - Overall column density normalization (default = 1.0). Convenient for
multiplying all the column densities simultaneously by the same factor.<br>
10: N_e - Total electron radial column density (to get Thomson depth).<br>
11: H_1 - Neutral hydrogen radial column density.  Number denotes number of bound
electrons.<br>
12: He_1 - Single-electron helium radial column density.<br>
13: He_2 - Neutral helium radial column density.<br>
14: C_1 - Single-electron C radial column density.<br>
<p>
.<br>
.<br>
.<br>
<p>

<br>
<br>
<br>

<h3>PHSI: </h3>
Here's an example of a set of model parameters for PHSI. 
(See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  phsi[1]
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   phsi       type                1.000     frozen
    2    2    1   phsi       Z                   8.000     frozen
    3    3    1   phsi       z                   2.000     frozen
    4    4    1   phsi       Nion     cm^-2     1.0000E+17 frozen
    5    5    1   phsi       Tion     eV         4.000     frozen
    6    6    1   phsi       redshift            0.000     frozen
    7    7    1   phsi       v_rad    km/s       0.000     frozen
    8    8    1   phsi       v_trans  km/s       0.000     frozen
    9    9    1   phsi       sig_rad  km/s       100.0     frozen
   10   10    1   phsi       sig_tran km/s       100.0     frozen
   11   11    1   phsi       INPUT               0.000     frozen
   12   12    1   phsi       INSHIFT?            0.000     frozen
   13   13    1   phsi       Gamma               2.000     frozen
   14   14    1   phsi       L_EMIN   keV       1.0000E-03 frozen
   15   15    1   phsi       L_EMAX   keV        100.0     frozen
   16   16    1   phsi       L_X      1e30e/s   1.0000E+14 frozen
   17   17    1   phsi       FLUXAVE             1.000     frozen
   18   18    1   phsi       f                  0.1000     frozen
   19   19    1   phsi       D        pc        1.4400E+07 frozen
   20   20    1   phsi       EMIN     keV       1.0000E-03 frozen
   21   21    1   phsi       EMAX     keV        15.00     frozen
   22   22    1   phsi       SPECBINS           1.0000E+05 frozen
   23   23    1   phsi       fileincr           -1.000     frozen
   24   24    1   phsi       verbose             1.000     frozen
   25   25    1   phsi       norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: type - =-1 to give the y-axis as dimensionless total opacity if
x-axis is in wavelength, =0 to give the y-axis as dimensionless total
opacity if x-axis is in energy, = 1 for pure absorption, = 2 for pure 
reemission, = 3 for pure reemission, recombination alone, = 4 for
absorption plus reemission (lower limit), = 5 for absorption plus 
reemission (upper limit), (
The following were designed for cataclysmic variable spectra : see <a href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2003ApJ...586L..77M&amp;db_key=AST&amp;high=3e7a2b9a7e18974" _base_target="_top">Mukai et al. 2003</a>)
= 6 for unobscured intrinsic continuum plus reemission spectrum, 
= 7 for unobscured intrinsic continuum plus reemission spectrum
assuming "infinite" radial velocity width, i.e., lines, but not edges,
are unsaturated at all column densities,
= 8 same as type=7, except without intrinsic continuum<br>
2: Z - Atomic number<br>
3: z - Number of electrons<br>
4: Nion - Ion column density in cm^-2<br>
5: Tion - Electron temperature [eV] for recombination contribution.<br>
6: redshift - Redshift of source<br>
7: v_rad - Radial velocity shift<br>
8: v_trans - Transverse velocity shift <br>
9: sig_rad - Radial velocity width (sigma)<br>
10: sig_tran - Transverse velocity width (sigma)<br>
11: INPUT - For inputting external spectrum (keep at default value of "0").<br>
12: INSHIFT? - For redshifting external spectrum (keep at default value of "0").<br>
13: Gamma - Power-law slope L(E)=AE^(-Gamma).<br>
14: L_EMIN - Low-energy limit [eV] to power law.<br>
15: L_EMAX - High-energy limit [eV] to power law.<br>
16: L_X - Total rest-frame luminosity (from L_EMIN [eV] to L_EMAX [eV]) in 
10^30 ergs/s. For non-zero redshift, cosmological correction is applied.<br>
17: FLUXAVE - This is the average flux of the intrinsic continuum (default = 1).  For highly
variable sources like Sy1 galaxies, this allows the user to determine the
"average" flux level to determine the proper level of reemission. <br>
18: f - Covering factor: f=Omega/4*Pi<br>
19: D - Distance to source in parsec.
If D is set to "0.", the distance is computed from the redshift with the Hubble law, using the cosmology set by XSPEC's <tt>cosmo</tt> command (H0 and lambda0; q0 is not used). The universe is taken to be flat, Omega_m = 1 - lambda0, so lambda0 must be less than 1. With XSPEC's default (H0 = 70, lambda0 = 0.73) that is Omega_m = 0.27. Versions before 2026 used a fixed H0 = 71 km/s/Mpc instead.<br>
20: EMIN - Minimum energy [keV] for internal grid.  This grid has nothing to do with the input luminosity spectrum.
For type&gt;=2 (calculation of reemission spectrum), make sure that energy 
range includes all regions with significant photoelectric absorption.<br>
</font>
21: EMAX - Maximum energy [keV] for internal grid.  This grid has nothing to do with the input luminosity spectrum.<br>
For type&gt;=2 (calculation of reemission spectrum), make sure that energy 
range includes all regions with significant photoelectric absorption (
EMAX &gt;= 15.0 keV should be sufficient).<br>
22: SPECBINS - Total number of energy bins (equally-spaced in energy) for internal grid.<br>
23: fileincr - &lt; 0 for no output files, &gt;= 0 for output files
are produced.  E.g., fileincr=22 would produce four files with output columns as follows:
E_spectrum_22.qdp (Observed E [keV], half-bin width [keV], and
spectrum [photons/cm^2/s/keV]), 
l_spectrum_22.qdp (Observed lambda [Angstrom], half-bin width [A], and
spectrum [photons/cm^2/s/A]), 
E_output_22.qdp (Observed E [eV], tau, L(E)/(4*Pi*D^2)
[photons/cm^2/s/eV], type1 spectrum [ph/cm^2/s/eV], type2 spectrum ["],
type3 spectrum ["], type4 spectrum ["], type5 spectrum ["], type6
spectrum ["], photoexcitation spectrum ["], RR spectrum ["], DR
spectrum ["]) <br>
l_output_22.qdp (Observed lambda [Angstrom], tau, L(lambda)/(4*Pi*D^2)
[photons/cm^2/s/A], type1 spectrum [ph/cm^2/s/A], type2 spectrum ["],
type3 spectrum ["], type4 spectrum ["], type5 spectrum ["], type6
spectrum ["], photoexcitation spectrum ["], RR spectrum ["], DR
spectrum ["]) <br>
24: verbose - =1 for output numbers/messages, =0 for no output numbers/messages<br>
<p>
<br>
<br>
<br>

<h3>PHXI: </h3>
Here's an example of a set of model parameters
for PHXI. This models uses
a user-defined distribution in ionization parameter. The file "xi.dat"
must exist in the directory you're running XSPEC in. See
photoion_dat/xi.dat for an example of this file. The
ionization parameter distribution is defined by simply connecting the
user-defined points in xi.dat with line segments and normalizing. The
"fractional ionic abundances" used were taken from an XSTAR simulation
of an extremely-low-column-density medium irradiated by a Gamma=2 power
law. (See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  phxi[1]
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   phxi       type                1.000     frozen
    2    2    1   phxi       N_H      cm^-2     1.0000E+22 frozen
    3    3    1   phxi       A_He     abund      1.000     frozen
    4    4    1   phxi       A_C      abund      1.000     frozen
    5    5    1   phxi       A_N      abund      1.000     frozen
    6    6    1   phxi       A_O      abund      1.000     frozen
    7    7    1   phxi       A_Ne     abund      1.000     frozen
    8    8    1   phxi       A_Mg     abund      1.000     frozen
    9    9    1   phxi       A_Al     abund      0.000     frozen
   10   10    1   phxi       A_Si     abund      1.000     frozen
   11   11    1   phxi       A_S      abund      1.000     frozen
   12   12    1   phxi       A_Ar     abund      0.000     frozen
   13   13    1   phxi       A_Ca     abund      0.000     frozen
   14   14    1   phxi       A_Fe     abund      1.000     frozen
   15   15    1   phxi       A_Ni     abund      0.000     frozen
   16   16    1   phxi       redshift            0.000     frozen
   17   17    1   phxi       v_rad    km/s       0.000     frozen
   18   18    1   phxi       v_trans  km/s       0.000     frozen
   19   19    1   phxi       sig_rad  km/s       100.0     frozen
   20   20    1   phxi       sig_tran km/s       100.0     frozen
   21   21    1   phxi       INPUT               0.000     frozen
   22   22    1   phxi       INSHIFT?            0.000     frozen
   23   23    1   phxi       Gamma               2.000     frozen
   24   24    1   phxi       L_EMIN   keV       1.0000E-03 frozen
   25   25    1   phxi       L_EMAX   keV        100.0     frozen
   26   26    1   phxi       L_X      1e30e/s   1.0000E+14 frozen
   27   27    1   phxi       FLUXAVE             1.000     frozen
   28   28    1   phxi       f                  0.1000     frozen
   29   29    1   phxi       D        pc        1.4400E+07 frozen
   30   30    1   phxi       EMIN     keV       1.0000E-03 frozen
   31   31    1   phxi       EMAX     keV        15.00     frozen
   32   32    1   phxi       SPECBINS           1.0000E+05 frozen
   33   33    1   phxi       fileincr           -1.000     frozen
   34   34    1   phxi       verbose             1.000     frozen
   35   35    1   phxi       norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: type - =-1 to give the y-axis as dimensionless total opacity if
x-axis is in wavelength, =0 to give the y-axis as dimensionless total
opacity if x-axis is in energy, = 1 for pure absorption, = 2 for pure 
reemission, = 3 for pure reemission, recombination alone, = 4 for
absorption plus reemission (lower limit), = 5 for absorption plus 
reemission (upper limit), (
The following were designed for cataclysmic variable spectra : see <a href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2003ApJ...586L..77M&amp;db_key=AST&amp;high=3e7a2b9a7e18974" _base_target="_top">Mukai et al. 2003</a>)
= 6 for unobscured intrinsic continuum plus reemission spectrum, 
= 7 for unobscured intrinsic continuum plus reemission spectrum
assuming "infinite" radial velocity width, i.e., lines, but not edges,
are unsaturated at all column densities,
= 8 same as type=7, except without intrinsic continuum<br>
2: N_H - Total hydrogen column density (neutral plus ionized).<br>
3: A_He - Overall helium abundance relative to "solar".<br>
.<br>
.<br>
.<br>
16: redshift - Redshift of source<br>
17: v_rad - Radial velocity shift <br>
18: v_trans - Transverse velocity shift <br>
19: sig_rad - Radial velocity width (sigma)<br>
20: sig_tran - Transverse velocity width (sigma)<br>
21: INPUT - For inputting external spectrum (keep at default value of "0").<br>
22: INSHIFT? - For redshifting external spectrum (keep at default value of "0").<br>
23: Gamma - Power-law slope L(E)=AE^(-Gamma).<br>
24: L_EMIN - Low-energy limit [eV] to power law.<br>
25: L_EMAX - High-energy limit [eV] to power law.<br>
26: L_X - Total rest-frame luminosity (from L_EMIN [eV] to L_EMAX [eV]) in 
10^30 ergs/s. For non-zero redshift, cosmological correction is applied.<br>
27: FLUXAVE - This is the average flux of the intrinsic continuum (default = 1).  For highly
variable sources like Sy1 galaxies, this allows the user to determine the
"average" flux level to determine the proper level of reemission.  <br>
28: f - Covering factor: f=Omega/4*Pi<br>
29: D - Distance to source in parsec<br>
If D is set to "0.", the distance is computed from the redshift with the Hubble law, using the cosmology set by XSPEC's <tt>cosmo</tt> command (H0 and lambda0; q0 is not used). The universe is taken to be flat, Omega_m = 1 - lambda0, so lambda0 must be less than 1. With XSPEC's default (H0 = 70, lambda0 = 0.73) that is Omega_m = 0.27. Versions before 2026 used a fixed H0 = 71 km/s/Mpc instead.<br>
30: EMIN - Minimum energy [keV] for internal grid.  This grid has nothing to do with the input luminosity spectrum.
For type&gt;=2 (calculation of reemission spectrum), make sure that energy 
range includes all regions with significant photoelectric absorption.<br>
31: EMAX - Maximum energy [keV] for internal grid.  This grid has nothing to do with the input luminosity spectrum.
For type&gt;=2 (calculation of reemission spectrum), make sure that energy 
range includes all regions with significant photoelectric absorption.<br>
EMAX &gt;= 15.0 keV should be sufficient.<br>
32: SPECBINS - Total number of energy bins (equally-spaced in energy) for internal grid.<br>
33: fileincr - &lt; 0 for no output files, &gt;= 0 for output files
are produced.  E.g., fileincr=22 would produce four files with output columns as follows:
E_spectrum_22.qdp (Observed E [keV], half-bin width [keV], and
spectrum [photons/cm^2/s/keV]), 
l_spectrum_22.qdp (Observed lambda [Angstrom], half-bin width [A], and
spectrum [photons/cm^2/s/A]), 
E_output_22.qdp (Observed E [eV], tau, L(E)/(4*Pi*D^2)
[photons/cm^2/s/eV], type1 spectrum [ph/cm^2/s/eV], type2 spectrum ["],
type3 spectrum ["], type4 spectrum ["], type5 spectrum ["], type6
spectrum ["], photoexcitation spectrum ["], RR spectrum ["], DR
spectrum ["]) <br>
l_output_22.qdp (Observed lambda [Angstrom], tau, L(lambda)/(4*Pi*D^2)
[photons/cm^2/s/A], type1 spectrum [ph/cm^2/s/A], type2 spectrum ["],
type3 spectrum ["], type4 spectrum ["], type5 spectrum ["], type6
spectrum ["], photoexcitation spectrum ["], RR spectrum ["], DR
spectrum ["]) <br>
34: verbose - =1 for output numbers/messages, =0 for no output numbers/messages<br>
35: norm - XSPEC internal 'normalization' parameter.  Leave this at
'1.000.'<br>
<p>
<br>
<br>
<br>

<h3>PHOTOION</h3>: 
Here's an example of a set of model parameters for PHOTOION.
(See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  photoion[1]
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   photoion   type                1.000     frozen
    2    2    1   photoion   redshift            0.000     frozen
    3    3    1   photoion   v_rad    km/s       0.000     frozen
    4    4    1   photoion   v_trans  km/s       0.000     frozen
    5    5    1   photoion   sig_rad  km/s       100.0     frozen
    6    6    1   photoion   sig_tran km/s       100.0     frozen
    7    7    1   photoion   INPUT               0.000     frozen
    8    8    1   photoion   INSHIFT?            0.000     frozen
    9    9    1   photoion   Gamma               2.000     frozen
   10   10    1   photoion   L_EMIN   keV       1.0000E-03 frozen
   11   11    1   photoion   L_EMAX   keV        100.0     frozen
   12   12    1   photoion   L_X      1e30e/s   1.0000E+14 frozen
   13   13    1   photoion   FLUXAVE             1.000     frozen
   14   14    1   photoion   f                  0.1000     frozen
   15   15    1   photoion   D        pc        1.4400E+07 frozen
   16   16    1   photoion   EMIN     keV       1.0000E-03 frozen
   17   17    1   photoion   EMAX     keV        15.00     frozen
   18   18    1   photoion   SPECBINS           1.0000E+05 frozen
   19   19    1   photoion   fileincr           -1.000     frozen
   20   20    1   photoion   verbose             1.000     frozen
   21   21    1   photoion   COLNORM             1.000     frozen
   22   22    1   photoion   N_e      cm^-2      0.000     frozen
   23   23    1   photoion   H_1      cm^-2      0.000     frozen
   24   24    1   photoion   He_1     cm^-2      0.000     frozen
   25   25    1   photoion   He_2     cm^-2      0.000     frozen
   26   26    1   photoion   C_1      cm^-2      0.000     frozen
   27   27    1   photoion   C_2      cm^-2      0.000     frozen
   28   28    1   photoion   C_3      cm^-2      0.000     frozen
   29   29    1   photoion   C_4      cm^-2      0.000     frozen
   30   30    1   photoion   C_5      cm^-2      0.000     frozen
   31   31    1   photoion   C_6      cm^-2      0.000     frozen
   32   32    1   photoion   N_1      cm^-2      0.000     frozen
   33   33    1   photoion   N_2      cm^-2      0.000     frozen
   34   34    1   photoion   N_3      cm^-2      0.000     frozen
   35   35    1   photoion   N_4      cm^-2      0.000     frozen
   36   36    1   photoion   N_5      cm^-2      0.000     frozen
   37   37    1   photoion   N_6      cm^-2      0.000     frozen
   38   38    1   photoion   N_7      cm^-2      0.000     frozen
   39   39    1   photoion   O_1      cm^-2      0.000     frozen
   40   40    1   photoion   O_2      cm^-2      0.000     frozen
   41   41    1   photoion   O_3      cm^-2      0.000     frozen
   42   42    1   photoion   O_4      cm^-2      0.000     frozen
   43   43    1   photoion   O_5      cm^-2      0.000     frozen
   44   44    1   photoion   O_6      cm^-2      0.000     frozen
   45   45    1   photoion   O_7      cm^-2      0.000     frozen
   46   46    1   photoion   O_8      cm^-2      0.000     frozen
   47   47    1   photoion   Ne_1     cm^-2      0.000     frozen
   48   48    1   photoion   Ne_2     cm^-2      0.000     frozen
   49   49    1   photoion   Ne_3     cm^-2      0.000     frozen
   50   50    1   photoion   Ne_4     cm^-2      0.000     frozen
   51   51    1   photoion   Ne_5     cm^-2      0.000     frozen
   52   52    1   photoion   Ne_6     cm^-2      0.000     frozen
   53   53    1   photoion   Ne_7     cm^-2      0.000     frozen
   54   54    1   photoion   Ne_8     cm^-2      0.000     frozen
   55   55    1   photoion   Ne_9     cm^-2      0.000     frozen
   56   56    1   photoion   Ne_10    cm^-2      0.000     frozen
   57   57    1   photoion   Mg_1     cm^-2      0.000     frozen
   58   58    1   photoion   Mg_2     cm^-2      0.000     frozen
   59   59    1   photoion   Mg_3     cm^-2      0.000     frozen
   60   60    1   photoion   Mg_4     cm^-2      0.000     frozen
   61   61    1   photoion   Mg_5     cm^-2      0.000     frozen
   62   62    1   photoion   Mg_6     cm^-2      0.000     frozen
   63   63    1   photoion   Mg_7     cm^-2      0.000     frozen
   64   64    1   photoion   Mg_8     cm^-2      0.000     frozen
   65   65    1   photoion   Mg_9     cm^-2      0.000     frozen
   66   66    1   photoion   Mg_10    cm^-2      0.000     frozen
   67   67    1   photoion   Mg_11    cm^-2      0.000     frozen
   68   68    1   photoion   Mg_12    cm^-2      0.000     frozen
   69   69    1   photoion   Al_1     cm^-2      0.000     frozen
   70   70    1   photoion   Al_2     cm^-2      0.000     frozen
   71   71    1   photoion   Al_3     cm^-2      0.000     frozen
   72   72    1   photoion   Al_4     cm^-2      0.000     frozen
   73   73    1   photoion   Al_5     cm^-2      0.000     frozen
   74   74    1   photoion   Al_6     cm^-2      0.000     frozen
   75   75    1   photoion   Al_7     cm^-2      0.000     frozen
   76   76    1   photoion   Al_8     cm^-2      0.000     frozen
   77   77    1   photoion   Al_9     cm^-2      0.000     frozen
   78   78    1   photoion   Al_10    cm^-2      0.000     frozen
   79   79    1   photoion   Al_11    cm^-2      0.000     frozen
   80   80    1   photoion   Al_12    cm^-2      0.000     frozen
   81   81    1   photoion   Al_13    cm^-2      0.000     frozen
   82   82    1   photoion   Si_1     cm^-2      0.000     frozen
   83   83    1   photoion   Si_2     cm^-2      0.000     frozen
   84   84    1   photoion   Si_3     cm^-2      0.000     frozen
   85   85    1   photoion   Si_4     cm^-2      0.000     frozen
   86   86    1   photoion   Si_5     cm^-2      0.000     frozen
   87   87    1   photoion   Si_6     cm^-2      0.000     frozen
   88   88    1   photoion   Si_7     cm^-2      0.000     frozen
   89   89    1   photoion   Si_8     cm^-2      0.000     frozen
   90   90    1   photoion   Si_9     cm^-2      0.000     frozen
   91   91    1   photoion   Si_10    cm^-2      0.000     frozen
   92   92    1   photoion   Si_11    cm^-2      0.000     frozen
   93   93    1   photoion   Si_12    cm^-2      0.000     frozen
   94   94    1   photoion   Si_13    cm^-2      0.000     frozen
   95   95    1   photoion   Si_14    cm^-2      0.000     frozen
   96   96    1   photoion   S_1      cm^-2      0.000     frozen
   97   97    1   photoion   S_2      cm^-2      0.000     frozen
   98   98    1   photoion   S_3      cm^-2      0.000     frozen
   99   99    1   photoion   S_4      cm^-2      0.000     frozen
  100  100    1   photoion   S_5      cm^-2      0.000     frozen
  101  101    1   photoion   S_6      cm^-2      0.000     frozen
  102  102    1   photoion   S_7      cm^-2      0.000     frozen
  103  103    1   photoion   S_8      cm^-2      0.000     frozen
  104  104    1   photoion   S_9      cm^-2      0.000     frozen
  105  105    1   photoion   S_10     cm^-2      0.000     frozen
  106  106    1   photoion   S_11     cm^-2      0.000     frozen
  107  107    1   photoion   S_12     cm^-2      0.000     frozen
  108  108    1   photoion   S_13     cm^-2      0.000     frozen
  109  109    1   photoion   S_14     cm^-2      0.000     frozen
  110  110    1   photoion   S_15     cm^-2      0.000     frozen
  111  111    1   photoion   S_16     cm^-2      0.000     frozen
  112  112    1   photoion   Ar_1     cm^-2      0.000     frozen
  113  113    1   photoion   Ar_2     cm^-2      0.000     frozen
  114  114    1   photoion   Ar_3     cm^-2      0.000     frozen
  115  115    1   photoion   Ar_4     cm^-2      0.000     frozen
  116  116    1   photoion   Ar_5     cm^-2      0.000     frozen
  117  117    1   photoion   Ar_6     cm^-2      0.000     frozen
  118  118    1   photoion   Ar_7     cm^-2      0.000     frozen
  119  119    1   photoion   Ar_8     cm^-2      0.000     frozen
  120  120    1   photoion   Ar_9     cm^-2      0.000     frozen
  121  121    1   photoion   Ar_10    cm^-2      0.000     frozen
  122  122    1   photoion   Ar_11    cm^-2      0.000     frozen
  123  123    1   photoion   Ar_12    cm^-2      0.000     frozen
  124  124    1   photoion   Ar_13    cm^-2      0.000     frozen
  125  125    1   photoion   Ar_14    cm^-2      0.000     frozen
  126  126    1   photoion   Ar_15    cm^-2      0.000     frozen
  127  127    1   photoion   Ar_16    cm^-2      0.000     frozen
  128  128    1   photoion   Ar_17    cm^-2      0.000     frozen
  129  129    1   photoion   Ar_18    cm^-2      0.000     frozen
  130  130    1   photoion   Ca_1     cm^-2      0.000     frozen
  131  131    1   photoion   Ca_2     cm^-2      0.000     frozen
  132  132    1   photoion   Ca_3     cm^-2      0.000     frozen
  133  133    1   photoion   Ca_4     cm^-2      0.000     frozen
  134  134    1   photoion   Ca_5     cm^-2      0.000     frozen
  135  135    1   photoion   Ca_6     cm^-2      0.000     frozen
  136  136    1   photoion   Ca_7     cm^-2      0.000     frozen
  137  137    1   photoion   Ca_8     cm^-2      0.000     frozen
  138  138    1   photoion   Ca_9     cm^-2      0.000     frozen
  139  139    1   photoion   Ca_10    cm^-2      0.000     frozen
  140  140    1   photoion   Ca_11    cm^-2      0.000     frozen
  141  141    1   photoion   Ca_12    cm^-2      0.000     frozen
  142  142    1   photoion   Ca_13    cm^-2      0.000     frozen
  143  143    1   photoion   Ca_14    cm^-2      0.000     frozen
  144  144    1   photoion   Ca_15    cm^-2      0.000     frozen
  145  145    1   photoion   Ca_16    cm^-2      0.000     frozen
  146  146    1   photoion   Ca_17    cm^-2      0.000     frozen
  147  147    1   photoion   Ca_18    cm^-2      0.000     frozen
  148  148    1   photoion   Ca_19    cm^-2      0.000     frozen
  149  149    1   photoion   Ca_20    cm^-2      0.000     frozen
  150  150    1   photoion   Fe_1     cm^-2      0.000     frozen
  151  151    1   photoion   Fe_2     cm^-2      0.000     frozen
  152  152    1   photoion   Fe_3     cm^-2      0.000     frozen
  153  153    1   photoion   Fe_4     cm^-2      0.000     frozen
  154  154    1   photoion   Fe_5     cm^-2      0.000     frozen
  155  155    1   photoion   Fe_6     cm^-2      0.000     frozen
  156  156    1   photoion   Fe_7     cm^-2      0.000     frozen
  157  157    1   photoion   Fe_8     cm^-2      0.000     frozen
  158  158    1   photoion   Fe_9     cm^-2      0.000     frozen
  159  159    1   photoion   Fe_10    cm^-2      0.000     frozen
  160  160    1   photoion   Fe_11    cm^-2      0.000     frozen
  161  161    1   photoion   Fe_12    cm^-2      0.000     frozen
  162  162    1   photoion   Fe_13    cm^-2      0.000     frozen
  163  163    1   photoion   Fe_14    cm^-2      0.000     frozen
  164  164    1   photoion   Fe_15    cm^-2      0.000     frozen
  165  165    1   photoion   Fe_16    cm^-2      0.000     frozen
  166  166    1   photoion   Fe_17    cm^-2      0.000     frozen
  167  167    1   photoion   Fe_18    cm^-2      0.000     frozen
  168  168    1   photoion   Fe_19    cm^-2      0.000     frozen
  169  169    1   photoion   Fe_20    cm^-2      0.000     frozen
  170  170    1   photoion   Fe_21    cm^-2      0.000     frozen
  171  171    1   photoion   Fe_22    cm^-2      0.000     frozen
  172  172    1   photoion   Fe_23    cm^-2      0.000     frozen
  173  173    1   photoion   Fe_24    cm^-2      0.000     frozen
  174  174    1   photoion   Fe_25    cm^-2      0.000     frozen
  175  175    1   photoion   Fe_26    cm^-2      0.000     frozen
  176  176    1   photoion   Ni_1     cm^-2      0.000     frozen
  177  177    1   photoion   Ni_2     cm^-2      0.000     frozen
  178  178    1   photoion   Ni_3     cm^-2      0.000     frozen
  179  179    1   photoion   Ni_4     cm^-2      0.000     frozen
  180  180    1   photoion   Ni_5     cm^-2      0.000     frozen
  181  181    1   photoion   Ni_6     cm^-2      0.000     frozen
  182  182    1   photoion   Ni_7     cm^-2      0.000     frozen
  183  183    1   photoion   Ni_8     cm^-2      0.000     frozen
  184  184    1   photoion   Ni_9     cm^-2      0.000     frozen
  185  185    1   photoion   Ni_10    cm^-2      0.000     frozen
  186  186    1   photoion   Ni_11    cm^-2      0.000     frozen
  187  187    1   photoion   Ni_12    cm^-2      0.000     frozen
  188  188    1   photoion   Ni_13    cm^-2      0.000     frozen
  189  189    1   photoion   Ni_14    cm^-2      0.000     frozen
  190  190    1   photoion   Ni_15    cm^-2      0.000     frozen
  191  191    1   photoion   Ni_16    cm^-2      0.000     frozen
  192  192    1   photoion   Ni_17    cm^-2      0.000     frozen
  193  193    1   photoion   Ni_18    cm^-2      0.000     frozen
  194  194    1   photoion   Ni_19    cm^-2      0.000     frozen
  195  195    1   photoion   Ni_20    cm^-2      0.000     frozen
  196  196    1   photoion   Ni_21    cm^-2      0.000     frozen
  197  197    1   photoion   Ni_21    cm^-2      0.000     frozen
  198  198    1   photoion   Ni_23    cm^-2      0.000     frozen
  199  199    1   photoion   Ni_24    cm^-2      0.000     frozen
  200  200    1   photoion   Ni_25    cm^-2      0.000     frozen
  201  201    1   photoion   Ni_26    cm^-2      0.000     frozen
  202  202    1   photoion   Ni_27    cm^-2      0.000     frozen
  203  203    1   photoion   Ni_28    cm^-2      0.000     frozen
  204  204    1   photoion   C_1_T    eV         4.000     frozen
  205  205    1   photoion   C_2_T    eV         2.500     frozen
  206  206    1   photoion   N_1_T    eV         4.000     frozen
  207  207    1   photoion   N_2_T    eV         3.000     frozen
  208  208    1   photoion   O_1_T    eV         10.00     frozen
  209  209    1   photoion   O_2_T    eV         4.000     frozen
  210  210    1   photoion   Ne_1_T   eV         10.00     frozen
  211  211    1   photoion   Ne_2_T   eV         10.00     frozen
  212  212    1   photoion   Ne_3_T   eV         10.00     frozen
  213  213    1   photoion   Ne_4_T   eV         10.00     frozen
  214  214    1   photoion   Ne_5_T   eV         10.00     frozen
  215  215    1   photoion   Ne_6_T   eV         10.00     frozen
  216  216    1   photoion   Ne_7_T   eV         10.00     frozen
  217  217    1   photoion   Ne_8_T   eV         10.00     frozen
  218  218    1   photoion   Ne_9_T   eV         10.00     frozen
  219  219    1   photoion   Ne_10_T  eV         10.00     frozen
  220  220    1   photoion   Mg_1_T   eV         10.00     frozen
  221  221    1   photoion   Mg_2_T   eV         10.00     frozen
  222  222    1   photoion   Mg_3_T   eV         10.00     frozen
  223  223    1   photoion   Mg_4_T   eV         10.00     frozen
  224  224    1   photoion   Mg_5_T   eV         10.00     frozen
  225  225    1   photoion   Mg_6_T   eV         10.00     frozen
  226  226    1   photoion   Mg_7_T   eV         10.00     frozen
  227  227    1   photoion   Mg_8_T   eV         10.00     frozen
  228  228    1   photoion   Mg_9_T   eV         10.00     frozen
  229  229    1   photoion   Mg_10_T  eV         10.00     frozen
  230  230    1   photoion   Al_1_T   eV         10.00     frozen
  231  231    1   photoion   Al_2_T   eV         10.00     frozen
  232  232    1   photoion   Al_3_T   eV         10.00     frozen
  233  233    1   photoion   Al_4_T   eV         10.00     frozen
  234  234    1   photoion   Al_5_T   eV         10.00     frozen
  235  235    1   photoion   Al_6_T   eV         10.00     frozen
  236  236    1   photoion   Al_7_T   eV         10.00     frozen
  237  237    1   photoion   Al_8_T   eV         10.00     frozen
  238  238    1   photoion   Al_9_T   eV         10.00     frozen
  239  239    1   photoion   Al_10_T  eV         10.00     frozen
  240  240    1   photoion   Si_1_T   eV         10.00     frozen
  241  241    1   photoion   Si_2_T   eV         10.00     frozen
  242  242    1   photoion   Si_3_T   eV         10.00     frozen
  243  243    1   photoion   Si_4_T   eV         10.00     frozen
  244  244    1   photoion   Si_5_T   eV         10.00     frozen
  245  245    1   photoion   Si_6_T   eV         10.00     frozen
  246  246    1   photoion   Si_7_T   eV         10.00     frozen
  247  247    1   photoion   Si_8_T   eV         10.00     frozen
  248  248    1   photoion   Si_9_T   eV         10.00     frozen
  249  249    1   photoion   Si_10_T  eV         10.00     frozen
  250  250    1   photoion   S_1_T    eV         10.00     frozen
  251  251    1   photoion   S_2_T    eV         10.00     frozen
  252  252    1   photoion   S_3_T    eV         10.00     frozen
  253  253    1   photoion   S_4_T    eV         10.00     frozen
  254  254    1   photoion   S_5_T    eV         10.00     frozen
  255  255    1   photoion   S_6_T    eV         10.00     frozen
  256  256    1   photoion   S_7_T    eV         10.00     frozen
  257  257    1   photoion   S_8_T    eV         10.00     frozen
  258  258    1   photoion   S_9_T    eV         10.00     frozen
  259  259    1   photoion   S_10_T   eV         10.00     frozen
  260  260    1   photoion   Ar_1_T   eV         10.00     frozen
  261  261    1   photoion   Ar_2_T   eV         10.00     frozen
  262  262    1   photoion   Ar_3_T   eV         10.00     frozen
  263  263    1   photoion   Ar_4_T   eV         10.00     frozen
  264  264    1   photoion   Ar_5_T   eV         10.00     frozen
  265  265    1   photoion   Ar_6_T   eV         10.00     frozen
  266  266    1   photoion   Ar_7_T   eV         10.00     frozen
  267  267    1   photoion   Ar_8_T   eV         10.00     frozen
  268  268    1   photoion   Ar_9_T   eV         10.00     frozen
  269  269    1   photoion   Ar_10_T  eV         10.00     frozen
  270  270    1   photoion   Ca_1_T   eV         10.00     frozen
  271  271    1   photoion   Ca_2_T   eV         10.00     frozen
  272  272    1   photoion   Ca_3_T   eV         10.00     frozen
  273  273    1   photoion   Ca_4_T   eV         10.00     frozen
  274  274    1   photoion   Ca_5_T   eV         10.00     frozen
  275  275    1   photoion   Ca_6_T   eV         10.00     frozen
  276  276    1   photoion   Ca_7_T   eV         10.00     frozen
  277  277    1   photoion   Ca_8_T   eV         10.00     frozen
  278  278    1   photoion   Ca_9_T   eV         10.00     frozen
  279  279    1   photoion   Ca_10_T  eV         10.00     frozen
  280  280    1   photoion   Fe_1_T   eV         10.00     frozen
  281  281    1   photoion   Fe_2_T   eV         10.00     frozen
  282  282    1   photoion   Fe_3_T   eV         10.00     frozen
  283  283    1   photoion   Fe_4_T   eV         10.00     frozen
  284  284    1   photoion   Fe_5_T   eV         10.00     frozen
  285  285    1   photoion   Fe_6_T   eV         10.00     frozen
  286  286    1   photoion   Fe_7_T   eV         10.00     frozen
  287  287    1   photoion   Fe_8_T   eV         10.00     frozen
  288  288    1   photoion   Fe_9_T   eV         10.00     frozen
  289  289    1   photoion   Fe_10_T  eV         10.00     frozen
  290  290    1   photoion   Ni_1_T   eV         10.00     frozen
  291  291    1   photoion   Ni_2_T   eV         10.00     frozen
  292  292    1   photoion   Ni_3_T   eV         10.00     frozen
  293  293    1   photoion   Ni_4_T   eV         10.00     frozen
  294  294    1   photoion   Ni_5_T   eV         10.00     frozen
  295  295    1   photoion   Ni_6_T   eV         10.00     frozen
  296  296    1   photoion   Ni_7_T   eV         10.00     frozen
  297  297    1   photoion   Ni_8_T   eV         10.00     frozen
  298  298    1   photoion   Ni_9_T   eV         10.00     frozen
  299  299    1   photoion   Ni_10_T  eV         10.00     frozen
  300  300    1   photoion   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: type - =-1 to give the y-axis as dimensionless total opacity if
x-axis is in wavelength, =0 to give the y-axis as dimensionless total
opacity if x-axis is in energy, = 1 for pure absorption, = 2 for pure 
reemission, = 3 for pure reemission, recombination alone, = 4 for
absorption plus reemission (lower limit), = 5 for absorption plus 
reemission (upper limit), (
The following were designed for cataclysmic variable spectra : see <a href="http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=2003ApJ...586L..77M&amp;db_key=AST&amp;high=3e7a2b9a7e18974" _base_target="_top">Mukai et al. 2003</a>)
= 6 for unobscured intrinsic continuum plus reemission spectrum, 
= 7 for unobscured intrinsic continuum plus reemission spectrum
assuming "infinite" radial velocity width, i.e., lines, but not edges,
are unsaturated at all column densities,
= 8 same as type=7, except without intrinsic continuum<br>
2: redshift - Redshift of source<br>
3: v_rad - Radial velocity shift <br>
4: v_trans - Transverse velocity shift <br>
5: sig_rad - Radial velocity width (sigma)<br>
6: sig_tran - Transverse velocity width (sigma)<br>
7: INPUT - For inputting external spectrum (keep at default value of "0").<br>
8: INSHIFT? - For redshifting external spectrum (keep at default value of "0").<br>
9: Gamma - Power-law slope L(E)=AE^(-Gamma).<br>
10: L_EMIN - Low-energy limit [eV] to power law.<br>
11: L_EMAX - High-energy limit [eV] to power law.<br>
12: L_X - Total rest-frame luminosity (from L_EMIN [eV] to L_EMAX [eV]) in 
10^30 ergs/s. For non-zero redshift, cosmological correction is applied.<br>
13: FLUXAVE - This is the average flux of the intrinsic continuum (default = 1).  For highly
variable sources like Sy1 galaxies, this allows the user to determine the
"average" flux level to determine the proper level of reemission.  <br>
14: f - Covering factor: f=Omega/4*Pi<br>
15: D - Distance to source in parsec. 
If D is set to "0.", the distance is computed from the redshift with the Hubble law, using the cosmology set by XSPEC's <tt>cosmo</tt> command (H0 and lambda0; q0 is not used). The universe is taken to be flat, Omega_m = 1 - lambda0, so lambda0 must be less than 1. With XSPEC's default (H0 = 70, lambda0 = 0.73) that is Omega_m = 0.27. Versions before 2026 used a fixed H0 = 71 km/s/Mpc instead.<br>
16: EMIN - Minimum energy [keV] for internal grid.  This grid has nothing to do with the input luminosity spectrum.
For type&gt;=2 (calculation of reemission spectrum), make sure that energy 
range includes all regions with significant photoelectric absorption.<br>
17: EMAX - Maximum energy [keV] for internal grid.  This grid has nothing to do with the input luminosity spectrum.
For type&gt;=2 (calculation of reemission spectrum), make sure that energy 
range includes all regions with significant photoelectric absorption.<br>
EMAX &gt;= 15.0 keV should be sufficient.<br>
18: SPECBINS - Total number of energy bins (equally-spaced in energy) for internal grid.<br>
19: fileincr - &lt; 0 for no output files, &gt;= 0 for output files
are produced.  E.g., fileincr=22 would produce four files with output columns as follows:
E_spectrum_22.qdp (Observed E [keV], half-bin width [keV], and
spectrum [photons/cm^2/s/keV]), 
l_spectrum_22.qdp (Observed lambda [Angstrom], half-bin width [A], and
spectrum [photons/cm^2/s/A]), 
E_output_22.qdp (Observed E [eV], tau, L(E)/(4*Pi*D^2)
[photons/cm^2/s/eV], type1 spectrum [ph/cm^2/s/eV], type2 spectrum ["],
type3 spectrum ["], type4 spectrum ["], type5 spectrum ["], type6
spectrum ["], photoexcitation spectrum ["], RR spectrum ["], DR
spectrum ["]) <br>
l_output_22.qdp (Observed lambda [Angstrom], tau, L(lambda)/(4*Pi*D^2)
[photons/cm^2/s/A], type1 spectrum [ph/cm^2/s/A], type2 spectrum ["],
type3 spectrum ["], type4 spectrum ["], type5 spectrum ["], type6
spectrum ["], photoexcitation spectrum ["], RR spectrum ["], DR
spectrum ["]) <br>
20: verbose - =1 for output numbers/messages, =0 for no output numbers/messages<br>
21: COLNORM - Overall column density normalization (default = 1.0). Convenient for
multiplying all the column densities simultaneously by the same factor.
<br>
22: N_e - Total electron radial column density (to get Thomson depth).<br>
23: H_1 - Neutral hydrogen radial column density.  Number denotes number of bound electrons.<br>
24: He_1 - Single-electron helium radial column density.  <br>
25: He_2 - Neutral helium radial column density.  <br>
26: C_1 - H-like C radial column density.  <br>
27: C_2 - He-like C radial column density.  <br>
.<br>
.<br>
.<br>
204: C_1_T - Electron temperature for recombinations forming H-like C.  <br>
204: C_2_T - Electron temperature for recombinations forming He-like C.  <br>
.<br>
.<br>
.<br>
300: norm - XSPEC internal 'normalization' parameter.  Leave this at '1.000.'<br></font>
<p>
<br>
<br>
<br>

<h3>ADDEXT: </h3>
Here's an example of a set of model parameters
for ADDEXT. ADDEXT allows 
the user to input an external spectrum located in the file
"addext.qdp".
(See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  addext[1]
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   addext     E_or_l              0.000     frozen
    2    2    1   addext     redshift            0.000     frozen
    3    3    1   addext     v        km/s       0.000     frozen
    4    4    1   addext     norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: E_or_l - =0 implies external file "addext.qdp" is in energy units [keV], =1 implies external file "addext.qdp" is in wavelength units [Angstrom].<br>
2: redshift - for redshifting external spectrum.<br>
3: v - Velocity shift.<br>
<p>
<br>
<br>
<br>

<h3>MULEXT: </h3>
 Here's an example of a set of model parameters
for MULEXT. MULEXT allows
the user to multiply any spectrum by an external "opacity" spectrum
located in file "mulext.qdp".
(See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  mulext[1]( powerlaw[2] )
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   mulext     E_or_l              0.000     frozen
    2    2    1   mulext     redshift            0.000     frozen
    3    3    1   mulext     v        km/s       0.000     frozen
    4    4    2   powerlaw   PhoIndex            2.000     +/-   0.000
    5    5    2   powerlaw   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: E_or_l - =0 implies external file "mulext.qdp" is in energy units
 [keV], =1 implies external file "mulext.qdp" is in wavelength units [Angstrom].<br>
2: redshift - For redshifting external spectrum.<br>
3: v - Velocity shift.<br>

<p>
<br>
<br>
<br>

<h3>TAUEXT: </h3>
 Here's an example of a set of model parameters
for TAUEXT. TAUEXT allows
the user to multiply any spectrum by an external "optical depth"
spectrum located in file "tauext.qdp".
(See below for a description of the parameters)
<pre>
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
  Model:  tauext[1]( powerlaw[2] )
  Model Fit Model Component  Parameter  Unit     Value
  par   par comp
    1    1    1   tauext     E_or_l              0.000     frozen
    2    2    1   tauext     redshift            0.000     frozen
    3    3    1   tauext     v        km/s       0.000     frozen
    4    4    1   tauext     tau_norm            1.000     frozen
    5    5    2   powerlaw   PhoIndex            2.000     +/-   0.000
    6    6    2   powerlaw   norm                1.000     +/-   0.000
  ---------------------------------------------------------------------------
  ---------------------------------------------------------------------------
</pre>
1: E_or_l - =0 implies external file "tauext.qdp" is in energy units
 [keV], =1 implies external file "tauext.qdp" is in wavelength units [Angstrom].<br>
2: redshift - For redshifting external spectrum.<br>
3: v - Velocity shift.<br>
4: tau_norm - Factor to multiply "tauext.qdp" values by.<br>

<p>
