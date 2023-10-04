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
from astropy import constants as c

from phoebe.backend import pyterpolmini

fluxes = None

def planck(T, lambda_):
    """
    Planck function, i.e., black-body intensity in J s^-1 sr^-1 m^-2 m^-1 units.

    """
    h = c.h.value
    c_ = c.c.value
    k_B = c.k_B.value
    return 2.0*h*c_**2/lambda_**5 / (np.exp(h*c_/(lambda_*k_B*T))-1.0)


def spe_simple(b, system, wavelengths=None, info={}, k=None):
    """
    Compute monochromatic flux F_nu.
    A simple model of uniform disk(s).

    b           .. Bundle object
    system      .. System object
    wavelengths .. wavelengths [m]
    info        .. dictionary w. 'original_index' to get wavelengths
    k           .. index to run computation

    Note: Applicable to detached non-eclipsing binaries.

    """
    global fluxes

    j = info['original_index']
    if k > 0:
        return {'flux': fluxes[j]}

    components = info['component']
    dataset = info['dataset']

    step = 0.01					# Ang
    angstroms = wavelengths*1.0e10		# Ang
    Lumtot = np.zeros(len(wavelengths))		# W
    fluxes = np.zeros(len(wavelengths))		# 1

    sg = pyterpolmini.SyntheticGrid(flux_type='relative', debug=False)

    for i, body in enumerate(system.bodies):

        rv = (system.vzi[i]*u.solRad/u.day).to('km/s').value	# km/s
        area = np.pi*(body.requiv*u.solRad.to('m'))**2		# m^2
        teff = body.teff					# K
        mass = body.masses[body.ind_self]			# M_S
        tmp = c.G*mass*u.solMass/(body.requiv*u.solRad)**2	# si
        logg = np.log10(tmp.cgs.value)				# cgs
        z = 10.0**body.abun					# 1

        props = {'teff': teff, 'logg': logg, 'z': z}

        s = sg.get_synthetic_spectrum(props, angstroms, order=2, step=step, padding=20.0)

        wave_ = pyterpolmini.doppler_shift(s.wave, rv)
        intens_ = pyterpolmini.interpolate_spectrum(wave_, s.intens, angstroms)

        for l in range(len(wavelengths)):

            Lum_lambda = area * planck(wavelengths[l], body.teff)

            fluxes[l] += Lum_lambda*intens_[l]
            Lumtot[l] += Lum_lambda

    fluxes /= Lumtot

    return {'flux': fluxes[0]}


def spe_integrate(b, system, wavelengths=None, info={}, k=None):
    """
    Compute monochromatic flux F_nu.
    A complex model w. integration over meshes.

    Note: See spe_simple().

    Note: All wavelengths are computed at once, but returned sequentially 0, 1, 2, ...

    Note: Applicable to contact or eclipsing binaries.

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

    # Note: metallicity should be per-triangle!
    z = 1.0					# 1
    step = 0.01					# Ang
    angstroms = wavelengths*1.0e10		# Ang
    fluxes = np.zeros(len(wavelengths))		# 1

    for i in range(len(abs_intensities)):
        if Lum[i] == 0.0:
            continue

        props = {'teff': teffs[i], 'logg': loggs[i], 'z': z}

        s = sg.get_synthetic_spectrum(props, angstroms, order=2, step=step, padding=20.0)

        wave_ = pyterpolmini.doppler_shift(s.wave, rvs[i]*1.0e-3)
        intens_ = pyterpolmini.interpolate_spectrum(wave_, s.intens, angstroms)

        fluxes += Lum[i]*intens_

    Lumtot = np.sum(Lum)
    fluxes /= Lumtot

    return {'flux': fluxes[0]}


spe = spe_integrate


