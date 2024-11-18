#!/usr/bin/env python3

"""
pyterpolmu.py
A C-accelerated version of Pyterpolmini.

Using Ndpolator from https://github.com/aprsa/ndpolator.

Reference: Nemravová et al. (2016, A&A 594, A55).
Reference: Brož (2017, ApJS 230, 19).

"""

import os
import numpy as np
import ndpolator

from astropy.constants import c
from scipy.signal import fftconvolve

ZERO_TOLERANCE = 1.0e-6

def doppler_shift(wave, rv):

    return wave * (1.0 + rv*1.0e3/c.value)

def interpolate_spectrum(wave, intens, wave_):
    """
    Fast linear interpolation

    Note: In previous version, cubic splines were used:

    tck = splrep(wave, intens, k=3)
    intens_ = splev(wave_, tck) 

    """

    intens_ = np.interp(wave_, wave, intens)

    return intens_

def instrumental_broadening(wave, intens, width=0.25, type='fwhm'):
    """
    Instrumental broadening; a convolution with a normal distribution.

    :param wave: wavelengths
    :param intens: intensities
    :param width: width in A
    :param type: either 'fwhh', or 'sigma'
    :return intens: the broadened spectrum

    """
    if width < ZERO_TOLERANCE:
        return intens

    if type == 'fwhm':
        sigma = width/2.3548
    elif type == 'sigma':
        sigma = width
    else:
        raise ValueError(("Unrecognised type='{}'").format(type))

    # Make sure the wavelengths are equidistant
    delta = np.diff(wave).min()
    Delta = wave.ptp()
    n = int(Delta/delta) + 1
    wave_ = np.linspace(wave[0], wave[-1], n)

    intens_ = interpolate_spectrum(wave, intens, wave_)

    # Construct the kernel!
    delta = wave_[1]-wave_[0]
    n_kernel = int(2*4*sigma/delta)

    if n_kernel == 0:
        raise ValueError(("Spectrum resolution too low for instrumental broadening (delta={}, width={}").format(delta, width))

    if n_kernel > n:
        raise ValueError(("Spectrum range too low for instrumental broadening"))

    wave_k = np.arange(n_kernel)*delta
    wave_k -= wave_k[-1]/2.0
    kernel = np.exp(-(wave_k)**2/(2.0*sigma**2))
    kernel /= sum(kernel)

    # Convolve the flux!
    intens_conv = fftconvolve(1.0-intens_, kernel, mode='same')

    if n_kernel%2 == 1:
        offset = 0.0
    else:
        offset = dwave/2.0

    intens = np.interp(wave+offset, wave_, 1.0-intens_conv, left=1, right=1)
    return intens

def rotational_broadening(wave, intens, vrot, epsilon=0.6):
    """
    Rotational broadening.

    :param wave: wavelengths
    :param intens: intensities
    :param vrot: projected rotational velocity in km/s
    :param epsilon: coefficient of linear limb-darkening
    :return intens: the rotated spectrum

    """
    if vrot < ZERO_TOLERANCE:
        return intens

    # Make sure the RVs are equidistant
    wave_log = np.log(wave)
    rv = np.linspace(wave_log[0], wave_log[-1], len(wave))
    step = rv[1] - rv[0]

    intens_rv = interpolate_spectrum(wave_log, intens, rv)

    vrot *= 1.0e3/c.value

    # Construct the kernel!
    n = int(np.ceil(2.0*vrot/step))
    rv_ker = np.arange(n)*step
    rv_ker = rv_ker - rv_ker[-1]/2.0
    y = 1.0 - (rv_ker/vrot)**2

    kernel = (2.0*(1.0-epsilon)*np.sqrt(y) + np.pi*epsilon/2.0*y) / (np.pi*vrot*(1.0-epsilon/3.0))
    kernel /= kernel.sum()

    # Convolve the flux!
    intens_conv = fftconvolve(1 - intens_rv, kernel, mode='same')

    if n % 2 == 1:
        rv = np.arange(len(intens_conv))*step + rv[0]
    else:
        rv = np.arange(len(intens_conv))*step + rv[0] - step/2.0

    wave_conv = np.exp(rv)

    intens = interpolate_spectrum(wave_conv, 1.0-intens_conv, wave)
    return intens


class Spectrum():
    """A synthetic spectrum."""

    def __init__(self, wave=None, intens=None):
        """
        Init a spectrum.

        """
        self.wave = wave
        self.intens = intens

    def load_spectrum(self, f=None):
        """
        Loads a spectrum.

        :param f: filename

        """
        if f is not None:
            self.filename = f

        print(("Loading file: " + str(self.filename)))

        # check if a binary file exists and -if true- load it
        # otherwise, load ascii (i.e., SLOW) and save as binary
        binary_file = self.filename + '.npy'

        if os.path.isfile(binary_file):
            data = np.load(binary_file, mmap_mode='r')
        else:
            data = np.loadtxt(self.filename, usecols=[0, 1], unpack=True)
            np.save(binary_file, data)

        self.wave = data[0]
        self.intens = data[1]

    def write_spectrum(self, f='spectrum.dat', fmt='%12.6f %12.8e'):
        """
        Writes a spectrum.

        :param f: filename
        :param fmt: format

        """
        header = "wave intens"
        np.savetxt(f, np.column_stack([self.wave, self.intens]), fmt=fmt, header=header)

    def truncate_spectrum(self, wmin=None, wmax=None):
        """
        Truncates a spectrum.

        :param wmin: minimum wavelength
        :param wmax: maximum wavelength
        """
        if wmin==None or wmax==None:
            return

        n = len(self.wave)
        w1 = self.wave[0]
        w2 = self.wave[-1]

        i = int((wmin-w1)/(w2-w1)*n + 0.0)
        j = int((wmax-w1)/(w2-w1)*n + 0.5)

        self.wave = self.wave[i:j]
        self.intens = self.intens[i:j]


class SyntheticGrid():
    """A grid of synthetic spectra."""

    def __init__(self, gridlist='gridlist', wmin=None, wmax=None):
        """
        Setup the grid.

        :param gridlist: list of files, their teff, logg, z, mu, ...

        """

        # read all filenames 
        self.files = np.loadtxt(gridlist, usecols=[0], dtype='str')

        # read wave's (cf. allocation below)
        s = Spectrum()
        s.load_spectrum(self.files[0])
        s.truncate_spectrum(wmin, wmax)
        self.wave = s.wave

        # read teffs, loggs -> N x M grid
        # Note: One needs to have 2 values of z!
        a = np.loadtxt(gridlist, usecols=[1, 2, 3], unpack=True)
        b = []
        for i in range(len(a)):
            tmp = np.unique(a[i])
            if len(tmp) < 2:
                tmp = np.array((tmp[0], tmp[0]+ZERO_TOLERANCE))
            b.append(tmp)
        b = np.array(b, dtype=object)
        c = np.arange(0,len(a[0]))

        grid = np.empty((len(b[0]), len(b[1]), len(b[2]), len(self.wave))) * np.nan

        # assign intens's <- cf. "voids"!
        for i in range(len(b[0])):
            for j in range(len(b[1])):
                for k in range(len(b[2])):
                    i_ = np.where(a[0] == b[0][i])[0]
                    j_ = np.where(a[1][i_] == b[1][j])[0]
                    k_ = np.where(a[2][i_][j_] == b[2][k])[0]
                    l_ = c[i_][j_][k_]
                    if len(l_) == 0:
                        continue
                    filename = self.files[l_][0]

                    s = Spectrum()
                    s.load_spectrum(filename)
                    s.truncate_spectrum(wmin, wmax)

                    grid[i, j, k, :] = s.intens

        # ndpolator instance
        self.ndp = ndpolator.Ndpolator(basic_axes=(b[0], b[1], b[2]))
        self.ndp.register(table='main', grid=grid, associated_axes=None)

    def get_synthetic_spectrum(self, props, wave, step=0.01, padding=20.0):
        """
        Computes a synthetic (interpolated) spectrum.

        Cf. SyntheticSpectrum().

        :param props: an array of values (teff, logg, z), at which to interpolate
        :param wave: wavelengths in A
        :param step: step in A
        :param padding: padding in A

        """
        props = np.array([props])
        wave = np.array(wave)

        interps = self.ndp.ndpolate(table='main', query_pts=props, extrapolation_method='linear')

        s = Spectrum()
        s.wave = self.wave
        s.intens = interps['interps'][0]

        return s

def main():
    teff = 6000.0
    logg = 4.25
    z = 1.0
    wmin = 6500.0
    wmax = 6600.0

    sg = SyntheticGrid(wmin=wmin, wmax=wmax)

    s = sg.get_synthetic_spectrum([teff, logg], [wmin, wmax], padding=0.0)

    s.write_spectrum('pyterpolmu.out')

if __name__ == "__main__":
    main()


