## ipolar30kk: emission from the accretion column of a magnetized white dwarf in intermediate polars (extension of ipolar models)

## Description

The models in the ipolar30kk.fits additive table model file represent the hard X-ray emission of the post-shock structures (PSR), or accretion columns, of magnetized white dwarfs in intermediate polars. In this version of the models we use mass-radius relation for the white dwarfs with the thick hydrogen envelopes with the surface temperature 30kK (Fontaine et al. 2001). The new models provide white dwarf masses that are approximately 0.04 solar masses higher in comparison with using "ipolar" models.

Accretion discs in this kind of objects are disrupted by the magnetic field of the compact object at the magnetospheric radius. Accreting matter then follows the magnetic field lines and forms a shock above the white dwarf surface. Heated matter cools down radiating the hard X-ray emission and settles onto the white dwarf surface.

The models of PSRs were computed in pseudo-dipole geometry ignoring the cyclotron cooling. The shape of a PSR spectrum depends on the free-fall velocity at the white dwarf surface, which in turn depends on two parameters: the white dwarf mass and the magnetospheric radius.

Only free-free emission from the 15 most abundant chemical elements is taken into account assuming solar element abundances. No emission lines or photo-recombination radiation are taken into account. As a result, the models should only be used to fit the hard X-ray emission (E > 15 keV) of intermediate polars.

The details are described in the papers [Suleimanov et al. 2016, A&A, 591, A35](https://ui.adsabs.harvard.edu/abs/2016A%26A...591A..35S/abstract),
[Suleimanov et al. 2025, A&A (arXiv:2506.03711)](https://arxiv.org/abs/2506.03711). Please, refer to these papers if you use this model.

The model grid is computed for 56 white dwarf masses from 0.3 to 1.4 solar masses with the step 0.02 solar masses and for 40 values of the relative magnetospheric radius (R_m/R_wd) which are uniformly distributed over the range of 1.5 to 60.

## Parameters

| Parameter | Description |
| --------- | ----------- |
| fall_heigh | the relative magnetospheric radius (R_m / R_wd) |
| M | the white dwarf mass in solar masses |
| norm | the normalization f (R_wd / d)^2, where d is the distance in cm, and f is the fraction of the white dwarf surface occupied by the PSR footprint. |

## Table Model File

The table model file can be downloaded as [ipolar30kk.fits](https://heasarc.gsfc.nasa.gov/FTP/software/xspec/localmodels/ipolar30kk/ipolar30kk.fits)

## Using the model

The file can be used by eg

<code>XSPEC12> model atable{ipolar30kk.fits}.</code>

