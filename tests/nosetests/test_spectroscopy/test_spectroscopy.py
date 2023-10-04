#!/usr/bin/env python3

import os
import sys
import numpy as np
import phoebe
from astropy import units

dir_ = os.path.dirname(os.path.realpath(__file__))

logger = phoebe.logger(clevel='INFO')
#logger = phoebe.logger(clevel='DEBUG')

b = phoebe.default_binary()

times, wavelengths, fluxes, sigmas = np.loadtxt(os.path.join(dir_, "Spe.dat"), usecols=[0, 1, 2, 3], unpack=True)

b.add_dataset('spe', times=times, wavelengths=wavelengths, fluxes=fluxes, sigmas=sigmas, passband='Johnson:R')

print("b['times@spe01@spe@dataset'] = ", b['times@spe01@spe@dataset'])
print("b['wavelengths@spe01@spe@dataset'] = ", b['wavelengths@spe01@spe@dataset'])
print("b['fluxes@spe01@spe@dataset'] = ", b['fluxes@spe01@spe@dataset'])

b.set_value('distance', context='system', value=100*units.pc)
b.set_value('ntriangles@primary', context='compute', value=100)
b.set_value('ntriangles@secondary', context='compute', value=100)

# Note: Using a tiny test grid.
from phoebe.backend import spectroscopy
from phoebe.backend import pyterpolmini

pyterpolmini.grid_directory = os.path.join(dir_, "grids")
pyterpolmini.grid_dict = dict(
    identification=[
        'test',
        ],
    directories=[
        ['TEST'],
        ],
    families=[
        ['test'],
        ],
    columns=[
        ['filename', 'teff', 'logg', 'z'],
        ],
    )

spectroscopy.sg = pyterpolmini.SyntheticGrid(mode='test', flux_type='relative', debug=False)

b.run_compute()

f = open('twigs.txt', 'w')
for twig in b.twigs:
  f.write("%s\n" % (twig))
f.close()

print("b['spe@spe01@phoebe01@latest@spe@model'] = ", b['spe@spe01@phoebe01@latest@spe@model'])
print("b['times@spe01@phoebe01@latest@spe@model'] = ", b['times@spe01@phoebe01@latest@spe@model'])
print("b['wavelengths@spe01@phoebe01@latest@spe@model'] = ", b['wavelengths@spe01@phoebe01@latest@spe@model'])
print("b['fluxes@spe01@phoebe01@latest@spe@model'] = ", b['fluxes@spe01@phoebe01@latest@spe@model'])

times = b['times@spe01@phoebe01@latest@spe@model'].value
wavelengths = b['wavelengths@spe01@phoebe01@latest@spe@model'].value
fluxes = b['fluxes@spe01@phoebe01@latest@spe@model'].value

np.savetxt('test_spectroscopy.out', np.c_[times, wavelengths, fluxes], header='times wavelenghts fluxes')

#b.plot(show=True)
#b.plot(x='wavelengths', marker='.', linestyle='none', show=True)


