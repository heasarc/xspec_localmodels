# Shakura–Sunyaev Spectral Energy Distribution: SSsed

## Description

The SSsed model was developed by Kubota et al. (<a href="https://academic.oup.com/mnras/article/528/2/1668/7513775">2024 MNRAS 528, 1668</a>) to describe spectral energy distribution (SED) of black hole binaries (BHBs). This model is a revised version of the agnsed model (Kubota &amp; Done 2018, MNRAS, 489, 524), but it was tuned to BHBs, and especially to their intermediate spectra where there are clearly two Compton components as well as (truncated) disc. The key concept is that the flow is radially stratified such that the accretion power is emitted as (colour-corrected) black body radiation at $r > r_{cor}$, while it is emitted as inverse-Comptonization by both the thermal and non-thermal corona at $r < r_{cor}$ The seed photons are assumed to be emitted from the underlying passive disc at $r < r_{cor}$. All the emission is constrained by the standard disc emissivity by Shakura &amp; Sunyaev (1973), with $\dot M$ constant with radius. If the components of the outer disc are visible, both parameters $r_{cor}$ and $r_{in}$ are determined independently. In cases where the outer disc is not visible, such as in the bright hard state, caution is needed in interpreting the obtained value of $r_{cor}$.


## Parameters in SSsed

Spectral parameters of the SSsed model are summarised in Table 1. The model has some switching parameters:

* If parameter 6 is negative, the model gives the inner hot Comptonisation component.
* If parameter 7 is negative, the model gives the Comptonisation component in the passive-disc corona region.
* If parameter 9 is negative, the model gives the outer disc.
* If parameter 12 is −1, the code will use the self gravity radius as calculate from Laor &amp; Netzer (1989, MNRAS, 238, 897).
* Colour correction is included when parameter 14 is set to 1, while it is not included when this parameter is set to 0. For BHB spectra, this parameter should be fixed at 1.

<table>
<caption>Table 1</caption>
<tr> <td>1.</td>   <td>mass</td>            <td>black hole mass in solar masses</td></tr>
<tr> <td>2.</td>   <td>dist</td>            <td>comoving (proper) distance in kpc</td></tr>
<tr> <td>3.</td>   <td>logmdot</td>         <td>$log (\dot m)$ where $\dot m$ = $\dot M / \dot M_{Edd}$ 
                                                and where $\dot M_{Edd}c^2 = L_{Edd}$</td></tr>
<tr> <td>4.</td>   <td>Rin</td>             <td>$r_{in}$, inner most radius of the accretion flow in $r_g$</td></tr>
<tr> <td>5.</td>   <td>cosi</td>            <td>$cos(i)$, inclination angle of the disc</td></tr>
<tr> <td>6.</td>   <td>kTe_th</td>          <td>$kT_{e,th}$, electron temperature for thermal corona in keV. 
                                                If this parameter is negative, the model gives 
                                                the inner hot Comptonisation component.</td></tr>
<tr> <td>7.</td>   <td>kTe_nt</td>          <td>$kT_{e,nth}$, apparent electron temperature for non-thermal 
                                                corona in keV which is recommended to be fixed at 300 keV 
                                                to mimic non-thermal electron distribution. If this 
                                                parameter is negative, the model gives the Comptonisation component 
                                                in the passive-disc corona region.</td></tr>
<tr> <td>8.</td>   <td>Gamma_th</td>        <td>$Γ_{th}$, photon index of inner hot corona. 
                                                If this parameter is negative, then only the inner 
                                                Compton component is used.</td></tr>
<tr> <td>9.</td>   <td>Gamma_nt</td>        <td>$Γ_{nth}$, photon index of disc-corona. 
                                                If this parameter is negative, 
                                                the model gives the outer disc.</td></tr>
<tr> <td>10.</td>  <td>frac_th</td>         <td>$f_{th}$, fraction of the hot Comptonising component to 
                                                the total Comptonisation</td></tr>
<tr> <td>11.</td>  <td>Rcor</td>            <td>$r_{cor}$, outer radius of the disc-corona region in $r_g$</td></tr>
<tr> <td>12.</td>  <td>logrout</td>         <td>$log (r_{out})$, outer radius of accretion disc in $r_g$.
                                                 If this parameter is −1, the code will use the 
                                                self gravity radius as calculated from Laor &amp; Netzer (1989)</td></tr>
<tr> <td>13.</td>  <td>redshift</td>        <td> must be fixed</td></tr>
<tr> <td>14.</td>  <td>color_cor</td>       <td>$f_{col}$, switching parameter for colour correction 
                                                (0: no colour correction, 
                                                 1: colour correction factor is calculated by the same way 
                                                as optxagnf (Done et al. 2012 MNRAS, 420, 1848)</td></tr>
<tr> <td>15.</td>  <td>norm</td>            <td>must be fixed at 1</td></tr>
</table>



## Installing the model

Follow the guidance in <a href="https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSappendixLocal.html">Appendix C in the Xspec manual</a> for installing local models.

For example: 

<code>$ xspec
XSPEC12>initpackage sssed /path/to/xspec_localmodels/SSsed/lmodel_sssed.dat /path/to/localmodels/xspec_localmodels/SSsed/
XSPEC12>lmod sssed /path/to/xspec_localmodels/SSsed/
</code>


