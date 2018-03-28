#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

f111 = np.genfromtxt('Gtau_111.dat', comments='#', usecols=(6,))
f121 = np.genfromtxt('Gtau_121.dat', comments='#', usecols=(6,))
f211 = np.genfromtxt('Gtau_211.dat', comments='#', usecols=(6,))
f221 = np.genfromtxt('Gtau_221.dat', comments='#', usecols=(6,))

plt.rc_context({'axes.edgecolor':'white', 'xtick.color':'white', 'ytick.color':'white', 'axes.labelcolor':'white', 'axes.linewidth':'2.4', 'xtick.labelsize':16, 'ytick.labelsize':16, 'axes.labelpad':10})

f = plt.figure(figsize=(12,9))

ax1 = f.add_subplot(321)
ax1.set_title(r'$\mathrm{top\;layer}$', fontsize=18, color='white')
ax1.plot(np.linspace(2.8,3.8,11), f111[:11], 'r', lw='6', marker='^', ms='9.0', label=r'$d_{xy}$')
ax1.plot(np.linspace(3.8,2.8,11), f111[11:], 'r', lw='6', marker='v', ms='9.0')
ax1.plot(np.linspace(2.8,3.8,11), f121[:11], 'lightblue', lw='6', marker='^', ms='9.0', label=r'$d_{xz}, d_{yz}$')
ax1.plot(np.linspace(3.8,2.8,11), f121[11:], 'lightblue', lw='6', marker='v', ms='9.0')
ax1.set_xlim(2.8,3.8)
ax1.set_ylim(0,0.05)
ax1.set_ylabel(r'$G(\tau=\beta / 2)$', fontsize=16, labelpad=12)
plt.axvspan(3.0, 3.3, facecolor='white', alpha=0.5)
legend = plt.legend(loc='upper right', frameon=False, fontsize=14)
plt.setp(legend.get_texts(), color='w')

ax1 = f.add_subplot(322)
ax1.set_title(r'$\mathrm{bottom\;layer}$', fontsize=18, color='white')
ax1.plot(np.linspace(2.8,3.8,11), f211[:11], 'r', lw='6', marker='^', ms='9.0', label=r'$d_{xy}$')
ax1.plot(np.linspace(3.8,2.8,11), f211[11:], 'r', lw='6', marker='v', ms='9.0')
ax1.plot(np.linspace(2.8,3.8,11), f221[:11], 'lightblue', lw='6', marker='^', ms='9.0', label=r'$d_{xz}, d_{yz}$')
ax1.plot(np.linspace(3.8,2.8,11), f221[11:], 'lightblue', lw='6', marker='v', ms='9.0')
ax1.set_xlim(2.8,3.8)
ax1.set_ylim(0,0.05)
plt.axvspan(3.0, 3.3, facecolor='white', alpha=0.5)
# ax1.set_ylabel(r'$G(\tau=\beta / 2)$')
legend = plt.legend(loc='upper right', frameon=False)
plt.setp(legend.get_texts(), color='w')


o111 = np.genfromtxt('occ_111.dat', comments='#', usecols=(6,))
o121 = np.genfromtxt('occ_121.dat', comments='#', usecols=(6,))
o211 = np.genfromtxt('occ_211.dat', comments='#', usecols=(6,))
o221 = np.genfromtxt('occ_221.dat', comments='#', usecols=(6,))

ax1 = f.add_subplot(323)
ax1.plot(np.linspace(2.8,3.8,11), o111[:11], 'r', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), o111[11:], 'r', lw='6', marker='v', ms='9.0')
ax1.plot(np.linspace(2.8,3.8,11), o121[:11], 'lightblue', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), o121[11:], 'lightblue', lw='6', marker='v', ms='9.0')
ax1.set_xlim(2.8,3.8)
ax1.set_ylim(0,0.5)
ax1.set_ylabel(r'$\left\langle \hat{c}_{i,\uparrow}^{\dagger}\hat{c}_{i,\uparrow} \right\rangle$', fontsize=16, labelpad=12)
plt.axvspan(3.0, 3.3, facecolor='white', alpha=0.5)

ax1 = f.add_subplot(324)
ax1.plot(np.linspace(2.8,3.8,11), o211[:11], 'r', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), o211[11:], 'r', lw='6', marker='v', ms='9.0')
ax1.plot(np.linspace(2.8,3.8,11), o221[:11], 'lightblue', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), o221[11:], 'lightblue', lw='6', marker='v', ms='9.0')
ax1.set_ylim(0,0.5)
ax1.set_xlim(2.8,3.8)
plt.axvspan(3.0, 3.3, facecolor='white', alpha=0.5)

d111 = np.genfromtxt('docc_111.dat', comments='#', usecols=(6,))
d121 = np.genfromtxt('docc_121.dat', comments='#', usecols=(6,))
d211 = np.genfromtxt('docc_211.dat', comments='#', usecols=(6,))
d221 = np.genfromtxt('docc_221.dat', comments='#', usecols=(6,))

ax1 = f.add_subplot(325)
ax1.plot(np.linspace(2.8,3.8,11), d111[:11], 'r', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), d111[11:], 'r', lw='6', marker='v', ms='9.0')
ax1.plot(np.linspace(2.8,3.8,11), d121[:11], 'lightblue', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), d121[11:], 'lightblue', lw='6', marker='v', ms='9.0')
ax1.set_xlim(2.8,3.8)
ax1.set_ylim(0,0.008)
plt.yticks([0,0.002,0.004,0.006,0.008])
ax1.set_xlabel(r"$U'\;[eV]$", fontsize=16)
ax1.set_ylabel(r'$\left\langle \hat{c}_{i,\uparrow}^{\dagger}\hat{c}_{i,\uparrow}\hat{c}_{i,\downarrow}^{\dagger}\hat{c}_{i,\downarrow} \right\rangle$', fontsize=16)
plt.axvspan(3.0, 3.3, facecolor='white', alpha=0.5)

ax1 = f.add_subplot(326)
ax1.plot(np.linspace(2.8,3.8,11), d211[:11], 'r', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), d211[11:], 'r', lw='6', marker='v', ms='9.0')
ax1.plot(np.linspace(2.8,3.8,11), d221[:11], 'lightblue', lw='6', marker='^', ms='9.0')
ax1.plot(np.linspace(3.8,2.8,11), d221[11:], 'lightblue', lw='6', marker='v', ms='9.0')
ax1.set_xlim(2.8,3.8)
ax1.set_ylim(0,0.008)
plt.yticks([0,0.002,0.004,0.006,0.008])
ax1.set_xlabel(r"$U'\;[eV]$", fontsize=16)
plt.axvspan(3.0, 3.3, facecolor='white', alpha=0.5, zorder=0)


plt.tight_layout()
plt.savefig('layer-b38-hysteresis-JAN-DPG.png', format='png', dpi=300, Rasterized=True, transparent=True)
plt.show()
