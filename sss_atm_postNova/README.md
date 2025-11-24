# sss_atm_PostNova: spectra of hot white dwarf in LTE and hydrostatic approximations for Classical Novae in their Super-Soft Stage

## Description

The models represent soft X-ray emission (0.07-1.4 keV) of the atmospheres of hot white dwarf, which are main emission component of the Classical Novae in their Super-soft stage. 

Super-soft sources are close binary systems with accreting white dwarfs. 
Accretion rates are so high that (quasi-)steady state thermonuclear burning on the white dwarf surfaces curried out. 
As a result the surfaces of such white dwarfs could be very hot, with effective temperatures up to 800-1000 kK. 

Classical Novae may undergo Super-Soft stage following the outburst. In this case, part of the matter, which has not been ejected from the system, 
may undergo stable thermonuclear burning on the surface of the white dwarf, producing soft X-ray emission. 
In the first approximation, the X-ray emission can be approximated by hot white dwarf atmosphere spectra.

Presented here are spectra of hot white dwarf model atmospheres were computed assuming plane-parallel and LTE approximations:
* The atmospheres were considered in hydrostatic equilibrium. Formally, the radiation pressure force exceed the gravity in the
  outer layers of the model atmospheres.  We postulated that in that layers the gas pressure equals to 10% the total pressure. 
* We considered 15 most abundant chemical elements: H, He, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe and Ni. 
* About 20 000 spectral lines of the considered elements were taken into account.
  The line parameters were taken from CHIANTI database (Dere et al. 1997, A&AS, 125, 149).

The chemical composition was adjusted for the Classical Nova AT 2018bej in the LMC in its Super-soft stage:
* The number fractions of H/He, $A({\rm H}), A({\rm He}) = 0.46, 0.54$.
* The number fraction of Carbon in solar units is the parameter of the model `chemC` = $A({\rm C}) = 0.0, 0.1, 0.2, 0.3, 0.4, 0.5$.
* The number fraction of Nitrogen in solar units is connected to Carbon:
```math
  A({\rm N}) = 0.5 + (0.5 - {\rm chemC}) \cdot \frac{A(C)_{\odot}}{A(N)_{\odot}} = 0.5 + (0.5 - {\rm chemC})\cdot 3.236.
```
  Thus some part of (0.5 − chemC) of Carbon is transformed into nitrogen during the CNO-cycle. 
  The abundance of carbon `chemC` changes from 0.0 (all C→N) to 0.5 (LMC composition).
* The abundance of the heavier elements is equal to 0.5 * solar.

The details are described in the paper Tavleev et al. 2024, A&A, 689, A335. Please, refer to this paper if you use this model.

The spectra are computed for 21 effective temperatures from 500 to 1000 kK with a step of 25 kK.


## Parameters

<p>
Fitting model parameters are:
</p>

| Parameter | Description |
| --------- | ----------- |
| chemC | the number fraction of Carbon in solar units | 
| Teff  | the effective temperature in kK              | 
| dlg   | the relative surface gravity                 |
| norm  | the normalization norm in units (1 km / 10 kpc)^2 |

The third parameter is the relative surface gravity dlg. This parameter is connected to the surface gravity g 
and the critical gravity for the given effective temperature g_Edd, dlg = g - log g_Edd, 
where log g_Edd = 4.818 + 4*log(Teff / 100 kK). The models with dlg = 0.1, 0.2, 0.4, 0.6, 1.0, and 1.4 were computed.


## Table model file
  
The model is available as an additive tabular model and stored in one file:

[sss_atm_postNova.fits](https://heasarc.gsfc.nasa.gov/FTP/software/xspec/localmodels/sss_atm_postNova/sss_atm_postNova.fits)

## Using the model

The model can be used by
 
<code>XSPEC12>model atable{sss_atm_postNova.fits}
</code>

## Contact

Contact the email listed in the FITS primary header for further information.