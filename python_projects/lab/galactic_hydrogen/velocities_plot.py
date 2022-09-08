# -*- coding: utf-8 -*-
"""
Created on Fri Nov 20 11:14:11 2020

plotting the velocitites obtained from the galactic hydrogen spectra against 
calculated distances 

@author: Ksenija
"""
import numpy as np
import matplotlib.pyplot as plt 

R_0 = 7.94 #kpc
unc_R0 = 0.42 #systematic 
kpc_to_m = 3.08567758 * 10 ** 19 #meters
VELOCITY_OF_SUN = 235 #km/s 
#(took from http://adsabs.harvard.edu/full/1996AstL...22..455K)
UNC_IN_SUN_VELOCITY = 20 #km/s systematic!
beam_width = 2.143 #deg systematic as well 

# http://adsabs.harvard.edu/full/1987AJ.....93.1090C 7.8 +/- 0.7 kpc
# R0 = 7.94 ± 0.42 kpc https://iopscience.iop.org/article/10.1086/380188/meta

FILE_NAME = 'velocities_data.txt'

INPUT_FILE = open(FILE_NAME, 'r')

longitudes = np.array([])
velocities = np.array([])
velocity_uncertainties = np.array([])

for line in INPUT_FILE:
    if line[0] != '%':

        split_up = line.split(',')
        longitudes = np.append(longitudes, float(split_up[0]))
        velocities = np.append(velocities, float(split_up[1]))
        velocity_uncertainties = np.append(velocity_uncertainties, 
                                           float(split_up[2]))

INPUT_FILE.close()

distances = R_0 * np.sin(np.deg2rad(longitudes))
real_velocities = velocities + VELOCITY_OF_SUN * np.sin(np.deg2rad(longitudes))


plt.errorbar(distances, real_velocities, velocity_uncertainties, 
             linestyle='None')
plt.scatter(distances, real_velocities)
plt.title('Velocity of Hydrogen clouds against their discance from the galactic center') 
plt.xlabel('Radial distance from the galactic center (kpc)')
plt.ylabel('orbital speed (km/s)')
plt.savefig('rotation_curve.png', dpi=300, bbox_inches='tight')
plt.show()
