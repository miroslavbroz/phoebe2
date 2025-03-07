#!/usr/bin/env python3

"""
spectroscopy.py

Spectroscopic module for Phoebe.

Reference: Nemravová et al. (2016, A&A 594, A55).
Reference: Brož (2017, ApJS 230, 19).

"""

import numpy as np

from phoebe import u, c
from phoebe import conf
from phoebe.backend import pyterpolmu

sg = None
sg2 = None
fluxes = None

def planck(lambda_, T=None):
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
        sg = pyterpolmu.SyntheticGrid(gridlist='gridlist')

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
        cosi = body.polar_direction_uvw[2]				# 1
        sini = np.sqrt(1.0-cosi**2)					# 1
        vrot = (omega*body.requiv*u.solRad.to('m')*sini)*1.0e-3		# km/s
        z = 10.0**body.abun						# 1

        props = [teff, logg, z]

        s = sg.get_synthetic_spectrum(props, angstroms, step=step, padding=20.0)

        wave_ = pyterpolmu.doppler_shift(s.wave, rv)				# Ang
        intens_ = pyterpolmu.rotational_broadening(wave_, s.intens, vrot)	# 1
        intens__ = pyterpolmu.interpolate_spectrum(wave_, intens_, angstroms)	# 1

        Lum = planck(wavelengths, T=teff)
        fluxes += Lum*area*intens__
        Lumtot += Lum*area

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
        sg = pyterpolmu.SyntheticGrid(gridlist='gridlist')

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
    # Note: intensity includes limb darkening (as passband)
    abs_intensities = meshes.get_column_flat('abs_intensities:{}'.format(dataset), components)
    mus = meshes.get_column_flat('mus', components)
    areas = meshes.get_column_flat('areas_si', components)
    rvs = (meshes.get_column_flat("rvs:{}".format(dataset), components)*u.solRad/u.d).to(u.m/u.s).value
    drvs = (meshes.get_column_flat("drvs", components)*u.solRad/u.d).to(u.m/u.s).value
    teffs = meshes.get_column_flat('teffs', components)
    loggs = meshes.get_column_flat('loggs', components)
    zs = 10.0**meshes.get_column_flat('abuns', components)

    fwhm = info['fwhm'] if info['use_instrumental'] else 0.0
    drvs *= 1.0 if info['use_rotational'] else 0.0

    Lum = abs_intensities*areas*mus*visibilities	# J s^-1 m^-1

    step = 0.01						# Ang
    angstroms = wavelengths*1.0e10			# Ang
    fluxes = np.zeros(len(wavelengths))			# 1

    for i in range(len(Lum)):
        if Lum[i] == 0.0:
            continue

        props = [teffs[i], loggs[i], zs[i]]

        s = sg.get_synthetic_spectrum(props, angstroms, step=step, padding=20.0)

        rv = rvs[i]*1.0e-3							# km/s
        drv = drvs[i]*1.0e-3							# km/s
        wave_ = pyterpolmu.doppler_shift(s.wave, rv)				# Ang
        intens_ = pyterpolmu.instrumental_broadening(wave_, s.intens, fwhm)	# 1
        intens__ = pyterpolmu.rotational_broadening(wave_, intens_, drv)	# 1
        intens___ = pyterpolmu.interpolate_spectrum(wave_, intens__, angstroms)	# 1

        fluxes += Lum[i]*intens___

        if conf.devel:
            f = open("spectroscopy.tmp", "a")
            np.savetxt(f, np.c_[len(angstroms)*[i], angstroms, intens___])
            f.write("\n")
            f.close()

    Lumtot = np.sum(Lum)
    fluxes /= Lumtot

    return {'flux': fluxes[0]}

########################################################################

def sed_simple(b, system, wavelengths=None, info={}, k=None):
    """
    Compute absolute monochromatic flux F_nu.

    Note: See spe_simple().

    """
    global sg2
    global fluxes

    if sg2 is None:
        sg2 = pyterpolmu.SyntheticGrid(gridlist='gridlist_ABS')

    j = info['original_index']
    if k > 0:
        return {'flux': fluxes[j]}

    components = info['component']
    dataset = info['dataset']

    d = system.distance				# m
    Lum = 1.0					# 1
    Lum /= d**2					# m^-2
    Lum /= np.pi				# m^-2

    step = 0.1					# Ang
    angstroms = wavelengths*1.0e10		# Ang
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

        props = [teff, logg, z]

        s = sg2.get_synthetic_spectrum(props, angstroms, step=step, padding=20.0)

        wave_ = pyterpolmu.doppler_shift(s.wave, rv)
        intens_ = pyterpolmu.rotational_broadening(wave_, s.intens, vrot)
        intens__ = pyterpolmu.interpolate_spectrum(wave_, intens_, angstroms)		# erg s^-1 cm^-2 Ang^-1
        intens__ *= 1.0e7								# W m^-2 m^-1

        fluxes += Lum*area*intens__

    return {'flux': fluxes[0]}


def sed_integrate(b, system, wavelengths=None, bandwidths=None, info={}, k=None):
    """
    Compute absolute monochromatic flux F_nu.

    Note: See spe_integrate().

    """
    global sg2
    global fluxes

    if sg2 is None:
        sg2 = pyterpolmu.SyntheticGrid(gridlist='gridlist_ABS')

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
    lds = meshes.get_column_flat("lds:{}".format(dataset), components)
    rvs = (meshes.get_column_flat("rvs:{}".format(dataset), components)*u.solRad/u.d).to(u.m/u.s).value
    drvs = (meshes.get_column_flat("drvs", components)*u.solRad/u.d).to(u.m/u.s).value
    teffs = meshes.get_column_flat('teffs', components)
    loggs = meshes.get_column_flat('loggs', components)
    zs = 10.0**meshes.get_column_flat('abuns', components)

    fwhm = info['fwhm'] if info['use_instrumental'] else 0.0
    drvs *= 1.0 if info['use_rotational'] else 0.0

    d = system.distance				# m
    Lum = lds*areas*mus*visibilities		# m^2
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

        props = [teffs[i], loggs[i], zs[i]]

        s = sg2.get_synthetic_spectrum(props, angstroms, step=step, padding=20.0)

        rv = rvs[i]*1.0e-3							# km/s
        drv = drvs[i]*1.0e-3							# km/s
        wave_ = pyterpolmu.doppler_shift(s.wave, rv)				# Ang
        intens_ = pyterpolmu.instrumental_broadening(wave_, s.intens, fwhm)	# erg s^-1 cm^-2 Ang^-1
        intens__ = pyterpolmu.rotational_broadening(wave_, intens_, drv)	# erg s^-1 cm^-2 Ang^-1
        intens___ = pyterpolmu.interpolate_spectrum(wave_, intens__, angstroms)	# erg s^-1 cm^-2 Ang^-1
        intens___ *= 1.0e7							# W m^-2 m^-1

        fluxes += Lum[i]*intens___

    return {'flux': fluxes[0]}


spe = spe_integrate
sed = sed_integrate


