#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

f = plt.figure(figsize=(9,4))
data = np.genfromtxt('bulk_extracted_path.hk', usecols=(0,1,2), dtype=np.float64, unpack=True)

ax1 = f.add_subplot(111)
ax1.set_xlim(0,160)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_ylabel(r'$E\;[\mathrm{eV}]$', fontsize=14)
plt.xticks([0,40,80,120,160], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$', r'$R$'], fontsize=14)
plt.plot(np.arange(161), data[0], 'red', lw='2', label=r'$d_{xy}$')
plt.plot(np.arange(161), data[1], 'darkmagenta', lw='2', label=r'$d_{xz}$')
plt.plot(np.arange(161), data[2], 'cornflowerblue', lw='2', label=r'$d_{yz}$')
# plt.plot(np.arange(121), np.zeros(121), 'black' , lw='1')
plt.axvline(x=40, color='black', lw='2')
plt.axvline(x=80, color='black', lw='2')
plt.axvline(x=120, color='black', lw='2')
plt.axhline(y=0, c='black', lw='1')

plt.legend(loc = 'best', frameon=False, fontsize=16)
plt.savefig('ham_bulk.eps', format='eps', bbox_inches='tight', dpi=600, rasterized=True)
plt.show()
