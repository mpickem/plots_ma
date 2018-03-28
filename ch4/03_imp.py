#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

o111 = np.genfromtxt('occ_111.dat', comments='#', usecols=(6,))
o121 = np.genfromtxt('occ_121.dat', comments='#', usecols=(6,))
o211 = np.genfromtxt('occ_211.dat', comments='#', usecols=(6,))
o221 = np.genfromtxt('occ_221.dat', comments='#', usecols=(6,))


f = plt.figure(figsize=(8,3))
ax1 = f.add_subplot(111)

# ax1.set_title(r'$\mathrm{top\;layer}$', fontsize=16)
ax1.plot(np.linspace(2.8,3.8,11), 2*o111[:11]+4*o121[:11], 'r', lw='2', marker='^', ms='7.0', label=r'$\mathrm{top\;layer}$')
ax1.plot(np.linspace(3.8,2.8,11), 2*o111[11:]+4*o121[11:], 'r', lw='2', marker='v', ms='7.0')
ax1.plot(np.linspace(2.8,3.8,11), 2*o211[:11]+4*o221[:11], 'b', lw='2', marker='^', ms='7.0', label=r'$\mathrm{bottom\;layer}$')
ax1.plot(np.linspace(3.8,2.8,11), 2*o211[11:]+4*o221[11:], 'b', lw='2', marker='v', ms='7.0')
ax1.plot(np.linspace(3.8,2.8,11), np.ones(11), c='black', lw='1')
ax1.set_xlim(2.8,3.8)
ax1.set_xlabel(r"$U'\;[eV]$")
# ax1.set_ylim(0,0.05)
ax1.set_ylabel(r'$n_\mathrm{imp}$', fontsize=14, labelpad=12)
plt.axvspan(3.0, 3.3, facecolor='black', alpha=0.2, zorder=0)
plt.legend(loc='best', frameon=False, fontsize=14)

plt.tight_layout()
plt.savefig('layer-charge-transfer.svg', bbox_inches='tight', format='svg', dpi=600, rasterized=True)
plt.show()
