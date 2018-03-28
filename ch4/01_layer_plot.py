#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

f = plt.figure(figsize=(9,4))
data = np.genfromtxt('path_2layer.hk', usecols=(1,2,3,4,5,6), dtype=np.float64, unpack=True)

ax1 = f.add_subplot(111)
ax1.set_xlim(0,120)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_ylabel(r'$E\;[\mathrm{eV}]$', fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
plt.plot(np.arange(121), data[0], 'r', lw='2', label='top layer')
plt.plot(np.arange(121), data[1], 'r', lw='2')
plt.plot(np.arange(121), data[2], 'r', lw='2')
plt.plot(np.arange(121), data[3], 'b', lw='2', label='bottom layer')
plt.plot(np.arange(121), data[4], 'b', lw='2')
plt.plot(np.arange(121), data[5], 'b', lw='2')
# plt.plot(np.arange(121), np.zeros(121), 'black' , lw='1')
plt.axvline(x=40, color='black', lw='2')
plt.axvline(x=80, color='black', lw='2')
plt.plot(np.arange(121), np.zeros(121), c='black', lw='1')

plt.legend(loc = 'best', frameon=False, fontsize=16)
plt.savefig('ham_2layer.eps', format='eps', bbox_inches='tight', dpi=600, rasterized=True)
plt.show()
