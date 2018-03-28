#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys

# make matplot math look like latex
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

# vq60 = h5py.File('../Output-u43-full-threeleg-gamma-newgiw/adga-20171221-174036.423-output.hdf5','r')
vq60 = h5py.File('../Output-u43-wo-gamma-newgiw/adga-20180125-211135.125-output.hdf5','r')

beta = 38.0
fmats = np.arange(np.pi/beta, 100, 2*np.pi/beta)

g = plt.figure(figsize=(11,6))
ax1 = g.add_subplot(231)
ax1.set_title(r'$\Gamma = (0,0,0)$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{Im}[\Sigma(\mathbf{k},i\nu)]$', fontsize=14)
# ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.text(2,-0.20,'$Z=0.38, \gamma=0.009$', color='red', fontsize=14)
ax1.text(2,-0.12,'$Z=0.40, \gamma=0.058$', color='black', fontsize=14)
ax1.set_ylim(-1.2,0.0)
# ax1.set_yticks([-1,-0.8,-0.6,-0.4,-0.2,0])
plt.plot(fmats[:60], vq60['input/siw'][0,2000:2060].imag, lw='4', ls='--', c='black')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][0,0,0,0,0,60:].imag, c='red', label='$d_{xy}$', lw='2')
plt.plot(fmats[:60], vq60['input/siw'][1,2000:2060].imag, lw='4', ls='--', c='gray')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][1,1,0,0,0,60:].imag, c='blue', label='$d_{xz}, d_{yz}$', lw='2')
plt.legend(loc='best', frameon=False)


ax1 = g.add_subplot(232)
ax1.set_title(r'$X = (0,\pi,0)$', fontsize=14)
# ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-1.2,0.0)
ax1.text(2,-0.12,'$Z=0.36, \gamma=0.068$', color='red', fontsize=14)
# ax1.set_yticks([-1,-0.8,-0.6,-0.4,-0.2,0])
plt.plot(fmats[:60], vq60['input/siw'][0,2000:2060].imag, lw='4', ls='--', c='black')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][0,0,0,40,0,60:].imag, c='red', label='$d_{xy}$', lw='2')
plt.plot(fmats[:60], vq60['input/siw'][1,2000:2060].imag, lw='4', ls='--', c='gray')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][1,1,0,40,0,60:].imag, c='blue', label='$d_{xz}$', lw='2')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][2,2,0,40,0,60:].imag, c='green', label='$d_{yz}$', lw='2')
plt.legend(loc='lower right', frameon=False)

ax1 = g.add_subplot(233)
ax1.set_title(r'$M = (\pi,\pi,0)$', fontsize=14)
ax1.text(2,-0.12,'$Z=0.35, \gamma=0.190$', color='red', fontsize=14)
# ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-1.2,0.0)
# ax1.set_yticks([-1,-0.8,-0.6,-0.4,-0.2,0])
plt.plot(fmats[:60], vq60['input/siw'][0,2000:2060].imag, lw='4', ls='--', c='black')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][0,0,40,40,0,60:].imag, c='red', label='$d_{xy}$', lw='2')
plt.plot(fmats[:60], vq60['input/siw'][1,2000:2060].imag, lw='4', ls='--', c='gray')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][1,1,40,40,0,60:].imag, c='blue', label='$d_{xz}, d_{yz}$', lw='2')
plt.legend(loc='lower right', frameon=False)

ax1 = g.add_subplot(234)
# ax1.set_title(r'$\Gamma = (0,0,0)$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{Re}[\Sigma(\mathbf{k},i\nu)]$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(0.5,2.5)
# ax1.set_yticks([-1,-0.8,-0.6,-0.4,-0.2,0])
plt.plot(fmats[:60], vq60['input/siw'][0,2000:2060].real, lw='4', ls='--', c='black')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][0,0,0,0,0,60:].real, c='red', label='$d_{xy}$', lw='2')
plt.plot(fmats[:60], vq60['input/siw'][1,2000:2060].real, lw='4', ls='--', c='gray')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][1,1,0,0,0,60:].real, c='blue', label='$d_{xz}, d_{yz}$', lw='2')
plt.legend(loc='best', frameon=False)


ax1 = g.add_subplot(235)
# ax1.set_title(r'$X = (0,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(0.5,2.5)
# ax1.set_yticks([-1,-0.8,-0.6,-0.4,-0.2,0])
plt.plot(fmats[:60], vq60['input/siw'][0,2000:2060].real, lw='4', ls='--', c='black')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][0,0,0,40,0,60:].real, c='red', label='$d_{xy}$', lw='2')
plt.plot(fmats[:60], vq60['input/siw'][1,2000:2060].real, lw='4', ls='--', c='gray')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][1,1,0,40,0,60:].real, c='blue', label='$d_{xz}$', lw='2')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][2,2,0,40,0,60:].real, c='green', label='$d_{yz}$', lw='2')
plt.legend(loc='lower right', frameon=False)

ax1 = g.add_subplot(236)
# ax1.set_title(r'$M = (\pi,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(0.5,2.5)
# ax1.set_yticks([-1,-0.8,-0.6,-0.4,-0.2,0])
plt.plot(fmats[:60], vq60['input/siw'][0,2000:2060].real, lw='4', ls='--', c='black')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][0,0,40,40,0,60:].real, c='red', label='$d_{xy}$', lw='2')
plt.plot(fmats[:60], vq60['input/siw'][1,2000:2060].real, lw='4', ls='--', c='gray')
plt.plot(fmats[:60], vq60['selfenergy/nonloc/dga'][1,1,40,40,0,60:].real, c='blue', label='$d_{xz}, d_{yz}$', lw='2')
plt.legend(loc='lower right', frameon=False)

plt.tight_layout()
plt.savefig('layer-2x3-ineq0.eps', bbox_inches='tight', format='eps', dpi=600, rasterized=True)

plt.show()
