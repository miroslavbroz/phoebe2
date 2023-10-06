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

sg = None
sg2 = None
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
    Compute relative monochromatic flux F_nu.
    A simple model of uniform disk(s).

    b           .. Bundle object
    system      .. System object
    wavelengths .. wavelengths [m]
    info        .. dictionary w. 'original_index' to get wavelengths
    k           .. index to run computation

    Note: All wavelengths are computed at once, but returned sequentially 0, 1, 2, ...

    Note: Applicable to detached non-eclipsing binaries.

    """
    global sg
    global fluxes

    if sg is None:
        sg = pyterpolmini.SyntheticGrid(flux_type='relative', debug=False)

    j = info['original_index']
    if k > 0:
        return {'flux': fluxes[j]}

    components = info['component']
    dataset = info['dataset']

    step = 0.01					# Ang
    angstroms = wavelengths*1.0e10		# Ang
    Lumtot = np.zeros(len(wavelengths))		# W
    fluxes = np.zeros(len(wavelengths))		# 1

    for i, body in enumerate(system.bodies):

        rv = -(system.vzi[i]*u.solRad/u.day).to('km/s').value		# km/s
        area = np.pi*(body.requiv*u.solRad.to('m'))**2			# m^2
        teff = body.teff						# K
        mass = body.masses[body.ind_self]				# M_S
        tmp = c.G*mass*u.solMass/(body.requiv*u.solRad)**2		# si
        logg = np.log10(tmp.cgs.value)					# cgs
        omega = body.freq_rot/u.day.to('s')				# rad/s
        sini = body.polar_direction_xyz[2]				# 1
        vrot = (omega*body.requiv*u.solRad.to('m')*sini)*1.0e-3		# km/s
        z = 10.0**body.abun						# 1

#        print("rv = ", rv)
#        print("teff = ", teff)
#        print("logg = ", logg)
#        print("omega = ", omega)
#        print("sini = ", sini)
#        print("vrot = ", vrot)
#        print("z = ", z)

        props = {'teff': teff, 'logg': logg, 'z': z}

        s = sg.get_synthetic_spectrum(props, angstroms, order=2, step=step, padding=20.0)

        wave_ = pyterpolmini.doppler_shift(s.wave, rv)
        intens_ = pyterpolmini.rotational_broadening(wave_, s.intens, vrot)
        intens__ = pyterpolmini.interpolate_spectrum(wave_, intens_, angstroms)

        for l in range(len(wavelengths)):

            Lum_lambda = area * planck(wavelengths[l], body.teff)

            fluxes[l] += Lum_lambda*intens__[l]
            Lumtot[l] += Lum_lambda

    fluxes /= Lumtot

    return {'flux': fluxes[0]}


def spe_integrate(b, system, wavelengths=None, info={}, k=None):
    """
    Compute relative monochromatic flux F_nu.
    A complex model w. integration over meshes.

    Note: See spe_simple().

    Note: Applicable to contact or eclipsing binaries.

    """
    global sg
    global fluxes

    if sg is None:
        sg = pyterpolmini.SyntheticGrid(flux_type='relative', debug=False)

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
    zs = 10.0**meshes.get_column_flat('abuns', components)

    Lum = abs_intensities*areas*mus*visibilities	# J s^-1 m^-1

    step = 0.01						# Ang
    angstroms = wavelengths*1.0e10			# Ang
    fluxes = np.zeros(len(wavelengths))			# 1

    for i in range(len(Lum)):
        if Lum[i] == 0.0:
            continue

        props = {'teff': teffs[i], 'logg': loggs[i], 'z': zs[i]}

        s = sg.get_synthetic_spectrum(props, angstroms, order=2, step=step, padding=20.0)

        rv = rvs[i]*1.0e-3							# km/s
        wave_ = pyterpolmini.doppler_shift(s.wave, rv)				# Ang
        intens_ = pyterpolmini.interpolate_spectrum(wave_, s.intens, angstroms)	# 1

        fluxes += Lum[i]*intens_

    Lumtot = np.sum(Lum)
    fluxes /= Lumtot

    return {'flux': fluxes[0]}


def sed_integrate(b, system, wavelengths=None, bandwidths=None, info={}, k=None):
    """
    Compute absolute monochromatic flux F_nu.

    Note: See spe_integrate().

    """
    global sg2
    global fluxes

    if sg2 is None:
        sg2 = pyterpolmini.SyntheticGrid(flux_type='absolute', debug=False)

    j = info['original_index']
    if k > 0:
        return {'flux': fluxes[j]}

    meshes = system.meshes
    components = info['component']
    dataset = info['dataset']

    visibilities = meshes.get_column_flat('visibilities', components)

    if np.all(visibilities==0):
        return {'flux': np.nan}

    mus = meshes.get_column_flat('mus', components)
    areas = meshes.get_column_flat('areas_si', components)
    rvs = (meshes.get_column_flat("rvs:{}".format(dataset), components)*u.solRad/u.d).to(u.m/u.s).value
    teffs = meshes.get_column_flat('teffs', components)
    loggs = meshes.get_column_flat('loggs', components)
    zs = 10.0**meshes.get_column_flat('abuns', components)

    d = system.distance				# m
    Lum = areas*mus*visibilities		# m^2
    Lum /= d**2					# 1

    # Note: a factor 1/pi is needed to obtain the solar values:
    # F_lambda ~ 2.e9 W m^-2 m^-1 (at Earth, 550 nm; Verbunt 2008)
    # F = \int F_lambda dlambda = 1363 W m^-2 (Kopp & Lean 2011)
    Lum /= np.pi				# 1

    step = 0.1					# Ang
    angstroms = wavelengths*1.0e10		# Ang
    fluxes = np.zeros(len(wavelengths))		# 1

    for i in range(len(Lum)):
        if Lum[i] == 0.0:
            continue

        props = {'teff': teffs[i], 'logg': loggs[i], 'z': zs[i]}

        s = sg2.get_synthetic_spectrum(props, angstroms, order=2, step=step, padding=20.0)

        rv = rvs[i]*1.0e-3							# km/s
        wave_ = pyterpolmini.doppler_shift(s.wave, rv)				# Ang
        intens_ = pyterpolmini.interpolate_spectrum(wave_, s.intens, angstroms)	# erg s^-1 cm^-2 Ang^-1 
        intens_ *= 1.0e7							# W m^-2 m^-1

        fluxes += Lum[i]*intens_

    return {'flux': fluxes[0]}


spe = spe_integrate
sed = sed_integrate


