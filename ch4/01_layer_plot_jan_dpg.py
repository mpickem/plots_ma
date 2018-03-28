#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

plt.rc_context({'axes.edgecolor':'white', 'xtick.color':'white', 'ytick.color':'white', 'axes.labelcolor':'white', 'axes.linewidth':'2.4'})
f = plt.figure(figsize=(9,4))
data = np.genfromtxt('path_2layer.hk', usecols=(1,2,3,4,5,6), dtype=np.float64, unpack=True)

ax1 = f.add_subplot(111)
ax1.set_xlim(0,120)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_ylabel(r'$E\;[\mathrm{eV}]$', fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
plt.plot(np.arange(121), np.zeros(121), c='white', lw='3')
plt.plot(np.arange(121), data[0], 'r', lw='6', label='top layer')
plt.plot(np.arange(121), data[1], 'r', lw='6')
plt.plot(np.arange(121), data[2], 'r', lw='6')
plt.plot(np.arange(121), data[3], 'lightblue', lw='6', label='bottom layer')
plt.plot(np.arange(121), data[4], 'lightblue', lw='6')
plt.plot(np.arange(121), data[5], 'lightblue', lw='6')
# plt.plot(np.arange(121), np.zeros(121), 'black' , lw='1')
plt.axvline(x=40, color='white', lw='3')
plt.axvline(x=80, color='white', lw='3')

legend = plt.legend(loc = 'best', frameon=False, fontsize=16)
plt.setp(legend.get_texts(), color='w')
plt.savefig('ham_2layer_JAN.png', format='png', bbox_inches='tight', dpi=300, rasterized=True, transparent=True)
plt.show()
