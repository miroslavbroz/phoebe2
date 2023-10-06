#!/usr/bin/env python3

import os
import sys
import numpy as np
import phoebe
from astropy import units

dir_ = os.path.dirname(os.path.realpath(__file__))

logger = phoebe.logger(clevel='INFO')
#logger = phoebe.logger(clevel='DEBUG')

b = phoebe.default_star()

times, wavelengths, fluxes, sigmas = np.loadtxt(os.path.join(dir_, "Sed.dat"), usecols=[0, 1, 2, 3], unpack=True)

b.add_dataset('sed', times=times, wavelengths=wavelengths, fluxes=fluxes, sigmas=sigmas)

print("b['times@sed01@sed@dataset'] = ", b['times@sed01@sed@dataset'])
print("b['wavelengths@sed01@sed@dataset'] = ", b['wavelengths@sed01@sed@dataset'])
print("b['fluxes@sed01@sed@dataset'] = ", b['fluxes@sed01@sed@dataset'])

b.set_value('distance', context='system', value=1*units.au)
b.set_value('teff@starA', context='component', value=5770.)
b.set_value('ntriangles@starA', context='compute', value=100)

# Note: Using a tiny test grid.
from phoebe.backend import spectroscopy
from phoebe.backend import pyterpolmini

pyterpolmini.grid_directory_ABS = os.path.join(dir_, "grids_ABS")
pyterpolmini.grid_dict_ABS = dict(
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

spectroscopy.sg2 = pyterpolmini.SyntheticGrid(mode='test', flux_type='absolute', debug=False)

b.run_compute()

f = open('twigs.txt', 'w')
for twig in b.twigs:
  f.write("%s\n" % (twig))
f.close()

print("b['sed@sed01@phoebe01@latest@sed@model'] = ", b['sed@sed01@phoebe01@latest@sed@model'])
print("b['times@sed01@phoebe01@latest@sed@model'] = ", b['times@sed01@phoebe01@latest@sed@model'])
print("b['wavelengths@sed01@phoebe01@latest@sed@model'] = ", b['wavelengths@sed01@phoebe01@latest@sed@model'])
print("b['fluxes@sed01@phoebe01@latest@sed@model'] = ", b['fluxes@sed01@phoebe01@latest@sed@model'])

times = b['times@sed01@phoebe01@latest@sed@model'].value
wavelengths = b['wavelengths@sed01@phoebe01@latest@sed@model'].value
fluxes = b['fluxes@sed01@phoebe01@latest@sed@model'].value

np.savetxt('test_spectroscopy2.out', np.c_[times, wavelengths, fluxes], header='times wavelenghts fluxes')

#b.plot(show=True)
#b.plot(x='wavelengths', marker='.', linestyle='none', show=True)


