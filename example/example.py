#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy import table, units as u
import cjam

# read in the positions from the example and add appropriate units
pos = table.QTable.read("xy.dat", format="ascii", names=["x","y"])
pos["x"].unit = u.arcsec
pos["y"].unit = u.arcsec

# read in the tracer density MGE from the example and add appropriate units
tracer_mge = table.QTable.read("mge.dat", format="ascii",
    names=["n","i","s","q"])
tracer_mge["i"].unit = 1/u.pc**2
tracer_mge["s"].unit = u.arcsec

# set anisotropy for each component of tracer MGE
beta = np.ones(len(tracer_mge))*0.01

# set rotation for each component of tracer MGE
kappa = np.array([0, 0.2, 0.5, 1.1, 0.6, 0, 0, 0])

# set distance
distance = 4.6*u.kpc

# set inclination
inclination = 0.87*u.rad

# set potential MGE, assumes mass-follows-tracers with mass/tracer ratio of 2.6
potential_mge = tracer_mge.copy()
potential_mge["i"] *= 2.6*u.Msun

# set BH mass and radius of BH gaussian
BH_mass = 100*u.Msun
BH_radius = 1*u.arcsec


# first run the model without rotation to test second moment calculations
# should run quite fast, will have 0 for all first moments
# second moments will be the same as the example as they do not depend on kappa
print("\nNo Rotation\n")
moments = cjam.axisymmetric(pos["x"], pos["y"], tracer_mge, potential_mge,
    distance, beta=beta, kappa=0, incl=inclination, mbh=BH_mass, rbh=BH_radius)

# set dispersions from first and second moments
moments["sigma_x"] = np.sqrt(moments["v2xx"]-moments["vx"]**2)
moments["sigma_y"] = np.sqrt(moments["v2yy"]-moments["vy"]**2)
moments["sigma_z"] = np.sqrt(moments["v2zz"]-moments["vz"]**2)


# now run with rotation to test first moment calculations
# will run a bit slower, should match the example!
print("\nWith Rotation\n")
moments_rot = cjam.axisymmetric(pos["x"], pos["y"], tracer_mge, potential_mge,
    distance, beta=beta, kappa=kappa, incl=inclination, mbh=BH_mass,
    rbh=BH_radius)

# set dispersions from first and second moments
moments_rot["sigma_x"] = np.sqrt(moments_rot["v2xx"]-moments_rot["vx"]**2)
moments_rot["sigma_y"] = np.sqrt(moments_rot["v2yy"]-moments_rot["vy"]**2)
moments_rot["sigma_z"] = np.sqrt(moments_rot["v2zz"]-moments_rot["vz"]**2)

# write moments out to file
table.Table(moments_rot).write("new_moments.dat", format="ascii.ecsv")

# read in the example moments file
example = table.QTable.read("moments.dat", format="ascii",
    names=["vx","vy","vz","v2xx","v2yy","v2zz","v2xy","v2xz","v2yz"])
example["vx"].unit = u.km/u.s
example["vy"].unit = u.km/u.s
example["vz"].unit = u.km/u.s
example["v2xx"].unit = (u.km/u.s)**2
example["v2yy"].unit = (u.km/u.s)**2
example["v2zz"].unit = (u.km/u.s)**2
example["v2xy"].unit = (u.km/u.s)**2
example["v2xz"].unit = (u.km/u.s)**2
example["v2yz"].unit = (u.km/u.s)**2

# convert proper motion components into mas/yr using distance
kms2masyr = (u.rad/distance).to("(mas s)/(km yr)")
example.replace_column("vx", example["vx"]*kms2masyr)
example.replace_column("vy", example["vy"]*kms2masyr)
example.replace_column("v2xx", example["v2xx"]*kms2masyr**2)
example.replace_column("v2yy", example["v2yy"]*kms2masyr**2)
example.replace_column("v2xy", example["v2xy"]*kms2masyr**2)
example.replace_column("v2xz", example["v2xz"]*kms2masyr)
example.replace_column("v2yz", example["v2yz"]*kms2masyr)


plt.figure(figsize=(7,6))
plt.subplots_adjust(0.1, 0.1, 0.98, 0.98, 0.5, 0.5)

plt.subplot(3,3,1)
plt.scatter(example["vx"], moments_rot["vx"], s=9)
plt.xlabel("example vx")
plt.ylabel("calculated vx")

plt.subplot(3,3,2)
plt.scatter(example["vy"], moments_rot["vy"], s=9)
plt.xlabel("example vy")
plt.ylabel("calculated vy")

plt.subplot(3,3,3)
plt.scatter(example["vz"], moments_rot["vz"], s=9)
plt.xlabel("example vz")
plt.ylabel("calculated vz")

plt.subplot(3,3,4)
plt.scatter(example["v2xx"], moments_rot["v2xx"], s=9)
plt.xlabel("example v2xx")
plt.ylabel("calculated v2xx")

plt.subplot(3,3,5)
plt.scatter(example["v2yy"], moments_rot["v2yy"], s=9)
plt.xlabel("example v2yy")
plt.ylabel("calculated v2yy")

plt.subplot(3,3,6)
plt.scatter(example["v2zz"], moments_rot["v2zz"], s=9)
plt.xlabel("example v2zz")
plt.ylabel("calculated v2zz")

plt.subplot(3,3,7)
plt.scatter(example["v2xy"], moments_rot["v2xy"], s=9)
plt.xlabel("example v2xy")
plt.ylabel("calculated v2xy")

plt.subplot(3,3,8)
plt.scatter(example["v2xz"], moments_rot["v2xz"], s=9)
plt.xlabel("example v2xz")
plt.ylabel("calculated v2xz")

plt.subplot(3,3,9)
plt.scatter(example["v2yz"], moments_rot["v2yz"], s=9)
plt.xlabel("example v2yz")
plt.ylabel("calculated v2yz")

plt.show()
