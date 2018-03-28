#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import h5py
import sys

# make matplot math look like latex
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

asy = h5py.File('../Output-u43-full-threeleg-gamma-newgiw/adga-20171221-174036.423-output.hdf5','r')

frequencies = 60
band = 0

beta = 38.0
fmats = np.linspace(np.pi/beta, 401*np.pi/beta, 201)
siwdga = np.zeros((6,6,80,80,1,2*frequencies), dtype=np.complex128)
asy['selfenergy/nonloc/dga'].read_direct(siwdga)
asy.close()


f = plt.figure(figsize=(11,6))
ax1 = f.add_subplot(231)
ax1.set_title(r'$\Gamma = (0,0,0)$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{Im}\Sigma(\mathbf{k},i\nu)$', fontsize=14)
# ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.06,0.02)
ax1.set_yticks([-0.06,-0.04,-0.02,0,0.02])
# plt.plot(fmats[:60], siwdga[0,1,0,0,0,60:].imag, label='$\mathrm{band: }1-2$', lw='2')
# plt.plot(fmats[:60], siwdga[0,2,0,0,0,60:].imag, label='$\mathrm{band: }1-3$', lw='2')
plt.plot(fmats[:60], siwdga[0,3,0,0,0,60:].imag, label='$d_{xy}(t) : d_{xy}(b)$', lw='2')
# plt.plot(fmats[:60], siwdga[0,4,0,0,0,60:].imag, label='$\mathrm{band: }1-5$', lw='2')
# plt.plot(fmats[:60], siwdga[0,5,0,0,0,60:].imag, label='$\mathrm{band: }1-6$', lw='2')

# plt.plot(fmats[:60], siwdga[1,2,0,0,0,60:].imag, label='$\mathrm{band: }2-3$', lw='2')
# plt.plot(fmats[:60], siwdga[1,3,0,0,0,60:].imag, label='$\mathrm{band: }2-4$', lw='2')
plt.plot(fmats[:60], siwdga[1,4,0,0,0,60:].imag, label='$d_{xz}(t) : d_{xz}(b)$', lw='2')
# plt.plot(fmats[:60], siwdga[1,5,0,0,0,60:].imag, label='$\mathrm{band: }2-6$', lw='2')

# plt.plot(fmats[:60], siwdga[2,3,0,0,0,60:].imag, label='$\mathrm{band: }3-4$', lw='2')
# plt.plot(fmats[:60], siwdga[2,4,0,0,0,60:].imag, label='$\mathrm{band: }3-5$', lw='2')
plt.plot(fmats[:60], siwdga[2,5,0,0,0,60:].imag, label='$d_{yz}(t) : d_{yz}(b)$', lw='2')

# plt.plot(fmats[:60], siwdga[3,4,0,0,0,60:].imag, label='$\mathrm{band: }4-5$', lw='2')
# plt.plot(fmats[:60], siwdga[3,5,0,0,0,60:].imag, label='$\mathrm{band: }4-6$', lw='2')

# plt.plot(fmats[:60], siwdga[4,5,0,0,0,60:].imag, label='$\mathrm{band: }5-6$', lw='2')
plt.legend(loc='lower right', frameon=False)

ax1 = f.add_subplot(232)
ax1.set_title(r'$X = (0,\pi,0)$', fontsize=14)
# ax1.set_ylabel(r'$\mathrm{Im}\Sigma(\mathbf{k},i\nu)$', fontsize=14)
# ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.06,0.02)
ax1.set_yticks([-0.06,-0.04,-0.02,0,0.02])
plt.plot(fmats[:60], siwdga[0,3,0,40,0,60:].imag, lw='2')
plt.plot(fmats[:60], siwdga[1,4,0,40,0,60:].imag, lw='2')
plt.plot(fmats[:60], siwdga[2,5,0,40,0,60:].imag, lw='2')
plt.legend(loc='best', frameon=False)

ax1 = f.add_subplot(233)
ax1.set_title(r'$M = (\pi,\pi,0)$', fontsize=14)
# ax1.set_ylabel(r'$\mathrm{Im}\Sigma(\mathbf{k},i\nu)$', fontsize=14)
# ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.06,0.02)
ax1.set_yticks([-0.06, -0.04,-0.02,0,0.02])
plt.plot(fmats[:60], siwdga[0,3,40,40,0,60:].imag, lw='2')
plt.plot(fmats[:60], siwdga[1,4,40,40,0,60:].imag, lw='2')
plt.plot(fmats[:60], siwdga[2,5,40,40,0,60:].imag, lw='2')
plt.legend(loc='best', frameon=False)

ax1 = f.add_subplot(234)
ax1.set_ylabel(r'$\mathrm{Re}\Sigma(\mathbf{k},i\nu)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.05,0.01)
ax1.set_yticks([-0.05, -0.03,-0.01,0.01])
# plt.plot(fmats[:60], siwdga[0,1,0,0,0,60:].imag, label='$\mathrm{band: }1-2$', lw='2')
# plt.plot(fmats[:60], siwdga[0,2,0,0,0,60:].imag, label='$\mathrm{band: }1-3$', lw='2')
plt.plot(fmats[:60], siwdga[0,3,0,0,0,60:].real, lw='2')
# plt.plot(fmats[:60], siwdga[0,4,0,0,0,60:].imag, label='$\mathrm{band: }1-5$', lw='2')
# plt.plot(fmats[:60], siwdga[0,5,0,0,0,60:].imag, label='$\mathrm{band: }1-6$', lw='2')

# plt.plot(fmats[:60], siwdga[1,2,0,0,0,60:].imag, label='$\mathrm{band: }2-3$', lw='2')
# plt.plot(fmats[:60], siwdga[1,3,0,0,0,60:].imag, label='$\mathrm{band: }2-4$', lw='2')
plt.plot(fmats[:60], siwdga[1,4,0,0,0,60:].real, lw='2')
# plt.plot(fmats[:60], siwdga[1,5,0,0,0,60:].imag, label='$\mathrm{band: }2-6$', lw='2')

# plt.plot(fmats[:60], siwdga[2,3,0,0,0,60:].imag, label='$\mathrm{band: }3-4$', lw='2')
# plt.plot(fmats[:60], siwdga[2,4,0,0,0,60:].imag, label='$\mathrm{band: }3-5$', lw='2')
plt.plot(fmats[:60], siwdga[2,5,0,0,0,60:].real, lw='2')

# plt.plot(fmats[:60], siwdga[3,4,0,0,0,60:].imag, label='$\mathrm{band: }4-5$', lw='2')
# plt.plot(fmats[:60], siwdga[3,5,0,0,0,60:].imag, label='$\mathrm{band: }4-6$', lw='2')

# plt.plot(fmats[:60], siwdga[4,5,0,0,0,60:].imag, label='$\mathrm{band: }5-6$', lw='2')
plt.legend(loc='lower right', frameon=False)

ax1 = f.add_subplot(235)
# ax1.set_ylabel(r'$\mathrm{Im}\Sigma(\mathbf{k},i\nu)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.05,0.01)
ax1.set_yticks([-0.05, -0.03,-0.01,0.01])
plt.plot(fmats[:60], siwdga[0,3,0,40,0,60:].real, lw='2')
plt.plot(fmats[:60], siwdga[1,4,0,40,0,60:].real, lw='2')
plt.plot(fmats[:60], siwdga[2,5,0,40,0,60:].real, lw='2')
plt.legend(loc='best', frameon=False)

ax1 = f.add_subplot(236)
# ax1.set_ylabel(r'$\mathrm{Im}\Sigma(\mathbf{k},i\nu)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.05,0.01)
ax1.set_yticks([-0.05, -0.03,-0.01,0.01])
plt.plot(fmats[:60], siwdga[0,3,40,40,0,60:].real, lw='2')
plt.plot(fmats[:60], siwdga[1,4,40,40,0,60:].real, lw='2')
plt.plot(fmats[:60], siwdga[2,5,40,40,0,60:].real, lw='2')
plt.legend(loc='best', frameon=False)


# print(fmats)

# gamma    = np.zeros((6,6,81,81),dtype = np.float64)
# Z        = np.zeros((6,6,81,81),dtype = np.float64)

# for iband in xrange(6):
#   for jband in xrange(6):
#     for ikx in xrange(80):
#         for iky in xrange(80):
#             for ikz in xrange(1):
#                 # this can be made more elegantly but
#                 # since we do the loop anyways ...
#                 polynomial = poly.polyfit(fmats[:2], siwdga[iband,jband,ikx,iky,ikz,frequencies:frequencies+2].imag, deg=1)
#                 gamma[iband,jband,ikx,iky] = -polynomial[0]
#                 Z[iband,jband,ikx,iky] = (1-polynomial[1])**(-1)

# sys.exit()

plt.tight_layout()
plt.savefig('layer-offdiag.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()
