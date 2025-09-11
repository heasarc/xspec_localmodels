# sss_atm: spectra of hot white dwarf in LTE and hydrostatic approximations for Super-Soft X-ray Sources

## Description

The sss_atm model represents soft X-ray emission (0.07-1.4 keV) of the atmospheres of hot white dwarf, which are main emission component of the super-soft X-ray sources. Super-soft sources are close binary systems with accreting white dwarfs. Accretion rates are so high that (quasi-)steady state thermonuclear burning on the white dwarf surfaces occurred. As a result the surfaces of such white dwarfs could be very hot, with effective temperatures up to 800-1000 kK. In the first approximation, their emission can be approximated with hot white dwarf atmosphere spectra.

Presented here are spectra of hot white dwarf model atmospheres which were computed assuming plane-parallel and LTE approximations. 
* The atmospheres were considered in hydrostatic equilibrium. Formally, the radiation pressure force exceeded the gravity in the outer layers of the model atmospheres.  We postulated that in those layers the gas pressure equals to 10% of the total pressure. 
* We considered 15 most abundant chemical elements, assuming solar H/He mix and solar (or proportional to solar) abundances of other elements, C, N, O, Ne, Na, Mg, Al, Si, S, Ar, Ca, Fe and Ni. Three heavy element abundances were considered, namely the solar one (A=1), one half of the solar (A=0.5, correspond to LMC abundance) and one tenth of the solar (A=0.1, correspond to SMC abundance). About 20,000 spectral lines of the considered elements were taken into account.
* The line parameters were taken from CHIANTI database (Dere et al. 1997, A&AS, 125, 149).

The details are described in the paper Suleimanov et al. 2024, A&A (arXiv: 2403.13557). Please, refer to this paper if you use this model.

The model grid for each chemical composition is computed for 37 effective temperatures from 100 to 1000 kK with a step of 25 kK.


## Parameters

<p>
Fitting model parameters are:
</p>

| Parameter | Description |
| --------- | ----------- |
| Teff | the effective temperature in kK | 
| dlg  | the relative surface gravity    |
| norm | the normalization norm in units (1 km / 10 kpc)^2 |

The second parameter is the relative surface gravity dlg. This parameter is connected to the surface gravity g and the critical gravity for the given effective temperature g_Edd, dlg = log g - log g_Edd, where log g_Edd = 4.818 + 4*log (Teff /100 kK). The models with dlg = 0.1, 0.2, 0.4, 0.6, 1.0, 1.4, 1.8 and 2.2 were computed.


## Table model files
  
The table model is available as an additive tabular model and stored in three files depending on the chemical composition:

| Table file name  | heavy element abundances  |
| ---------------  | ------------------------  |
| [sss_atm_sol.fits](https://heasarc.gsfc.nasa.gov/FTP/software/xspec/localmodels/sss_atm/sss_atm_sol.fits) | for solar one (A=1)       |
| [sss_atm_LMC.fits](https://heasarc.gsfc.nasa.gov/FTP/software/xspec/localmodels/sss_atm/sss_atm_LMC.fits) | for LMC abundance (A=0.5) |
| [sss_atm_SMC.fits](https://heasarc.gsfc.nasa.gov/FTP/software/xspec/localmodels/sss_atm/sss_atm_SMC.fits) | for SMC abundance (A=0.1) |


## Using the model

These files can be used by e.g.
 
<code>XSPEC12>model atable{sss_atm_sol.fits}
</code>

## Contact

Contact the email listed in the FITS primary header for further information.

