#!/usr/bin/env python3

f = open("Spe.dat", "w")
f.write("# time wavelength flux sigma\n")

t1 = 0.0            # d
t2 = 0.25           # d
dt = 0.25           # d
wave1 = 6500.0e-10  # m
wave2 = 6600.0e-10  # m
dwave = 1.0e-10     # m
flux = 1.0          # 1
sigma = 0.01        # 1

t = t1
while t < t2+0.5*dt:

    wave = wave1
    while wave < wave2+0.5*dwave:
       f.write("%.8f  %.8e  %.8f  %.8f\n" % (t, wave, flux, sigma))
       wave += dwave

    t += dt

