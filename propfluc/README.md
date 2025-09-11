# Propagating fluctuations in the mass accretion rate of a precessing flow as a power spectral model for black hole binaries: propfluc

This model and description are due to <a href="https://www2.physics.ox.ac.uk/contacts/people/ingrama">Adam Ingram</a>.

## Description

<p>
The power spectrum of black hole binaries (BHBs) in the rise to outburst
generally consists of a quasi-periodic oscillation (QPO) and additional
band limited noise. This can be phenomenologically modelled using a number
of broad Lorentzians for the band limited noise and a narrow Lorentzian
for each harmonic of the QPO. Here, we assume a physical origin for this
variability in a truncated disc / hot inner flow geometry. Mass accretion
rate fluctuations are generated everywhere in the flow (primarily) at the
local viscous frequency and propagate inwards towards the black hole. We
also assume that the entire flow is precessing due to frame-dragging
(Lense-Thirring precession), which gives rise to a QPO. The key assumption
is the surface density profile as this sets the precession (and therefore QPO)
frequency but also, by mass conservation, sets the viscous frequency as a
function of radius. Further details of the model, propfluc, are included in
<a href="http://arxiv.org/abs/1108.0789">Ingram & Done (2011)</a>,
particularly in the appendix. This version of the model is calculated
analytically using the formalism of <a href=" http://adsabs.harvard.edu/abs/2013MNRAS.434.1476I">Ingram & van der Klis (2013)</a> and
is consequently a <b>lot</b> faster than the original model.

<p>
We assume that the surface density is given by a smoothly broken power
law, consistent with the results of general relativistic magneto
hydrodynamic (GRMHD) simulations. The break occurs at the bending wave
radius, rbw, with the power law dependence on radius parametrised by
zeta for r&gt;&gt;rbw and by lambda for r&lt;&lt;rbw. The sharpness of
the break is given by another parameter, kappa. The QPO frequency is
calculated self-consistently and the power spectrum of the QPO is
represented by Lorentzians centred at f_QPO, 2f_QPO, 3f_QPO and
1/2f_QPO to represent the 1st, 2nd, 3rd and sub harmonics
respectively.

<p>
There is some extra functionality to this version of the
model. Parameters 15 and 16 are respectively `hard' and `soft' band
emissivity indices (gamma_h and gamma_s). Parameter 20 (mode)
specifies whether the output is the hard band power spectrum (mode=1),
soft band power spectrum (mode=2) or the time lag between the two
energy bands (mode=3). If gamma_h &gt; gamma_s, the model predicts
hard lags. Note, if mode=1, gamma_s should be fixed and if mode=2,
gamma_h should be fixed (to anything, the fit is insensitive to the
parameter). Parameter 21 (conv) specifies whether the QPO signal is
added to the broad band noise signal (conv=0) or multiplied with it
(conv=1). The latter is more physically motivated (<a href="
http://adsabs.harvard.edu/abs/2013MNRAS.434.1476I">Ingram & van der
Klis 2013</a>). Handy tips for working with the model are provided in
<a href="http://adsabs.harvard.edu/abs/2014MNRAS.440.2882R">Rapisarda,
Ingram & van der Klis (2014)</a>.

<p>

As this is a model for the power spectrum (or lag spectrum if mode=3)
rather than the spectral energy distribution, the process for loading
the data into XSPEC is slightly different. First of all, a power
spectrum can easily be created from a light curve using <a
href="https://heasarc.gsfc.nasa.gov/docs/software/xronos/examples/powspec.html">powspec</a>
from the <a
href="https://heasarc.gsfc.nasa.gov/docs/software/xronos/xronos.html">XRONOS</a>
package (for example). The power spectrum will then typically be
written in the form<br>
f, df, P, dP<br>
where f is the frequency, P the power and df and dP denote the corresponding
error. XSPEC, however, expects a .pha file in the form<br>
Emin, Emax, F(Emax-Emin), dF(Emax-Emin)<br>
where Emin and Emax are the lower and upper bands of each energy bin and F
is the flux. It is therefore necessary to create a data file with the inputs<br>
f-df, f+df, 2Pdf, 2dPdf.<br>
This data file can then be converted into a .pha file using <a href="https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/ftflx2xsp.html">ftflx2xsp</a> which
will also generate a diagonal response (.rsp) file. The data can then be
read into XSPEC in the usual way, using the command data 1:1 'filename'.pha.
The command ip euf will then present the data and model in terms of
fP plotted against f (even though, by default the axes will be labelled in
the units of EF vs E). More information about this procedure can be found
in appendix A in <a href="http://arxiv.org/abs/1108.0789">Ingram &
Done (2011)</a>. A similar process can be used for converting a lag
spectrum into a format which can be read in by XSPEC.

## Parameters in propfluc

<p>
<table>
<tr><td>1.</td><td>Sigma0</td><td>normalisation of surface density profile.</td></tr>
<tr><td>2.</td><td>rbw</td><td>bending wave radius.</td></tr>
<tr><td>3.</td><td>kappa</td><td>surface density profile parameter.</td></tr>
<tr><td>4.</td><td>lambda</td><td>surface density profile parameter.</td></tr>
<tr><td>5.</td><td>zeta</td><td>surface density profile parameter.</td></tr>
<tr><td>6.</td><td>Fvar</td><td>fractional variability per decade in radius.</td></tr>
<tr><td>7.</td><td>ro</td><td>truncation radius.</td></tr>
<tr><td>8.</td><td>ri</td><td>inner radius of the flow.</td></tr>
<tr><td>9.</td><td>Q</td><td>quality factor, Q=centroid/FWHM, of the QPO harmonics.</td></tr>
<tr><td>10.</td><td>Qsub</td><td>quality factor of the subharmonic (can be different from Q).</td></tr>
<tr><td>11.</td><td>n_qpo</td><td>normalisation of fundamental.</td></tr>
<tr><td>12.</td><td>n_2qpo</td><td>normalisation of 2nd harmonic.</td></tr>
<tr><td>13.</td><td>n_3qpo</td><td>normalisation of 3rd harmonic.</td></tr>
<tr><td>14.</td><td>n_05qpo</td><td>normalisation of sub-harmonic.</td></tr>
<tr><td>15.</td><td>gamma_h</td><td>hard band emissivity index.</td></tr>
<tr><td>16.</td><td>gamma_s</td><td>soft band emissivity indenx.</td></tr>
<tr><td>17.</td><td>M</td><td>BH mass in Solar masses.</td></tr>
<tr><td>18.</td><td>a</td><td>dimensionless spin parameter.</td></tr>
<tr><td>19.</td><td>Ndec</td><td>number of rings per decade in radius.</td></tr>
<tr><td>20.</td><td>mode</td><td>1 = hard band PSD, 2 = soft band PSD, 3 = lag spectrum.</td></tr>
<tr><td>21.</td><td>conv</td><td>0 = add QPO, 1 = multiply QPO.</td></tr>
<tr><td>22.</td><td>Norm</td><td>included by XSPEC - fix to unity.</td></tr>
</table>
The final four parameters are all just settings, not free parameters.
<p>


## Installing the model

This model requires heasoft 6.26 or newer, with Xspec patch 12.10.1i or newer.


