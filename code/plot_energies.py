#!/usr/bin/env python

#
# PYTHON SCRIPT FOR PLOTTING THE ENERGY LEVELS, USING MATPLOTLIB
# Made by Simon Nilsson 2016-11-10
#


import matplotlib.pylab as plt
import numpy as np

filename = "energies.dat"


data = np.loadtxt(filename)

plt.figure(figsize=(8,6))
plt.plot(data[:,0], data[:,1], '-')

# Set labels
plt.xlabel('Time / [ps]')
plt.ylabel('Energy / [eV]')
plt.title('$E(t)$')


# Axis limits
#plt.xlim([0, 100])

# Tick fontsize
plt.xticks(fontisize=12)
plt.yticks(fontsize=12)


# Save and display the plot
plt.savefig('energies.pdf')
plt.show()
