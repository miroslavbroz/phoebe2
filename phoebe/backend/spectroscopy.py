#!/usr/bin/env python3

"""
spectroscopy.py

Spectroscopic module for Phoebe.

Reference: Nemravová et al. (2016, A&A 594, A55).
Reference: Brož (2017, ApJS 230, 19).

"""

import sys
import numpy as np
from astropy import units as u

from phoebe.backend import pyterpolmini

fluxes = None

def spe_integrate(b, system, wavelengths=None, info={}, k=None):
    """
    Compute monochromatic flux F_nu.

    Note: All wavelengths are computed at once, but returned sequentially 0, 1, 2, ...

    """
    global fluxes

    j = info['original_index']
    if k > 0:
        return {'flux': fluxes[j]} 

    meshes = system.meshes
    components = info['component']
    dataset = info['dataset']

    visibilities = meshes.get_column_flat('visibilities', components)

    if np.all(visibilities==0):
        return {'flux': np.nan}

    # Note: intensity should be per-wavelength!
    abs_intensities = meshes.get_column_flat('abs_intensities:{}'.format(dataset), components)
    mus = meshes.get_column_flat('mus', components)
    areas = meshes.get_column_flat('areas_si', components)
    rvs = (meshes.get_column_flat("rvs:{}".format(dataset), components)*u.solRad/u.d).to(u.m/u.s).value
    teffs = meshes.get_column_flat('teffs', components)
    loggs = meshes.get_column_flat('loggs', components)

    Lum = abs_intensities*areas*mus*visibilities		# J s^-1 m^-1

    sg = pyterpolmini.SyntheticGrid(flux_type='relative', debug=False)

    z = 1.0					# 1
    step = 0.01					# Ang
    angstroms = wavelengths*1.0e10		# Ang
    fluxes = np.zeros(len(wavelengths))		# 1

    for i in range(len(abs_intensities)):
        if Lum[i] == 0.0:
            continue

        props = {'teff': teffs[i], 'logg': loggs[i], 'z': z}

        c = sg.get_synthetic_spectrum(props, angstroms, order=2, step=step, padding=20.0)

        wave_ = pyterpolmini.doppler_shift(c.wave, rvs[i]*1.0e-3)
        intens_ = pyterpolmini.interpolate_spectrum(wave_, c.intens, angstroms)

        fluxes += Lum[i]*intens_

    Lumtot = np.sum(Lum)
    fluxes /= Lumtot

    return {'flux': fluxes[0]}


spe = spe_integrate


