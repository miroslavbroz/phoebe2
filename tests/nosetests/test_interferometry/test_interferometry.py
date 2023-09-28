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

times, u, v, wavelengths, vises, sigmas = np.loadtxt(os.path.join(dir_, "Vis.dat"), usecols=[0, 1, 2, 3, 5, 6], unpack=True)

b.add_dataset('vis', times=times, u=u, v=v, wavelengths=wavelengths, vises=vises, sigmas=sigmas)

b.set_value('distance', context = 'system', value=100*units.pc)

b.run_compute()

f = open('twigs.txt', 'w')
for twig in b.twigs:
  f.write("%s\n" % (twig))
f.close()

# Note: Model (synthetic) baselines and wavelengths must be copied ex-post!
b.set_value('u@latest@model', value=b.get_value('u@dataset'), ignore_readonly=True)
b.set_value('v@latest@model', value=b.get_value('v@dataset'), ignore_readonly=True)
b.set_value('wavelengths@latest@model', value=b.get_value('wavelengths@dataset'), ignore_readonly=True)

print("b['vis@vis01@phoebe01@latest@vis@model'] = ", b['vis@vis01@phoebe01@latest@vis@model'])
print("b['times@vis01@phoebe01@latest@vis@model'] = ", b['times@vis01@phoebe01@latest@vis@model'])
print("b['u@vis01@phoebe01@latest@vis@model'] = ", b['u@vis01@phoebe01@latest@vis@model'])
print("b['v@vis01@phoebe01@latest@vis@model'] = ", b['v@vis01@phoebe01@latest@vis@model'])
print("b['wavelengths@vis01@phoebe01@latest@vis@model'] = ", b['wavelengths@vis01@phoebe01@latest@vis@model'])
print("b['vises@vis01@phoebe01@latest@vis@model'] = ", b['vises@vis01@phoebe01@latest@vis@model'])

times = b['times@vis01@phoebe01@latest@vis@model'].value
u = b['u@vis01@phoebe01@latest@vis@model'].value
v = b['v@vis01@phoebe01@latest@vis@model'].value
wavelengths = b['wavelengths@vis01@phoebe01@latest@vis@model'].value
vises = b['vises@vis01@phoebe01@latest@vis@model'].value

np.savetxt('test_interferometry.out', np.c_[times, u, v, wavelengths, vises], header='times u v wavelenghts vises')

b.plot(show=True)
b.plot(x='u', show=True)
b.plot(x='v', show=True)


