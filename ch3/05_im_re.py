#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import numpy.polynomial.polynomial as poly
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import h5py
import sys

# make matplot math look like latex
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

asy = h5py.File('output-srvo3-beta38-n240/adga-20171209-114252.312-output.hdf5','r')
beta = 38.0
fmats = np.arange(np.pi/beta,100,2*np.pi/beta)

nfreq = 128

## Polynomial fit to extract Z and gamma -- we use a polynomial third order for 5 frequencies
polynomial = poly.polyfit(fmats[:2], asy['input/siw'][0,2000:2002].imag, deg=1)
print('DMFT - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:2], asy['selfenergy/nonloc/dga'][0,0,0,0,0,240:242].imag, deg=1)
print('DGA-G-xy-xz-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:2], asy['selfenergy/nonloc/dga'][0,0,0,10,0,240:242].imag, deg=1)
print('DGA-X-xy-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:2], asy['selfenergy/nonloc/dga'][1,1,0,10,0,240:242].imag, deg=1)
print('DGA-X-xz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:2], asy['selfenergy/nonloc/dga'][0,0,10,10,0,240:242].imag, deg=1)
print('DGA-M-xy - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:2], asy['selfenergy/nonloc/dga'][1,1,10,10,0,240:242].imag, deg=1)
print('DGA-M-xz-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:2], asy['selfenergy/nonloc/dga'][0,0,10,10,10,240:242].imag, deg=1)
print('DGA-R-xy-xz-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))

if False: # plot the approximation if one wants to
    f = plt.figure()
    ax1 = f.add_subplot(111)
    plt.plot(fmats[:16], siw_new_g_i[:16])
    plt.plot(np.arange(55)/10, poly.polyval(np.arange(55)/10, polynomial, tensor=False))
    plt.show()
    sys.exit()

f = plt.figure(figsize=(12,6))
ax1 = f.add_subplot(241)
ax1.set_title(r'$\Gamma = (0,0,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{Im}[\Sigma(i\nu)]$', fontsize=14)
ax1.set_xlim((0,20))
ax1.set_ylim((-0.9,0))
plt.yticks([-0.9,-0.7,-0.5,-0.3,-0.1])
ax1.text(5,-0.16,'$Z=0.48, \gamma=0.019$', color='blue', fontsize=14)
ax1.text(5,-0.1,'$Z=0.48, \gamma=0.015$', color='gray', fontsize=14)
p1 = plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].imag, label=r'$\mathrm{DMFT}$', lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,0,0,240:240+nfreq].imag, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='', c='blue')

plt.legend(loc='lower right', frameon=False)

ax1 = f.add_subplot(242)
ax1.set_title(r'$X = (0,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.text(5,-0.16,'$Z=0.48, \gamma=0.019$', color='red', fontsize=14)
ax1.text(5,-0.1,'$Z=0.47, \gamma=0.018$', color='blue', fontsize=14)
ax1.set_xlim((0,20))
ax1.set_ylim((-0.9,0))
plt.yticks([-0.9,-0.7,-0.5,-0.3,-0.1])
p1 = plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].imag, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,10,0,240:240+nfreq].imag, label='$d_{xy},d_{yz}$', lw='2', marker='',c='blue')
p3 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,0,10,0,240:240+nfreq].imag, label='$d_{xz}$', lw='2', marker='',c='red')
plt.legend(loc='lower right', frameon=False)

ax1 = f.add_subplot(243)
ax1.set_title(r'$M = (\pi,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.text(5,-0.16,'$Z=0.46, \gamma=0.006$', color='red', fontsize=14)
ax1.text(5,-0.10,'$Z=0.47, \gamma=0.018$', color='blue', fontsize=14)
ax1.set_xlim((0,20))
ax1.set_ylim((-0.9,0))
plt.yticks([-0.9,-0.7,-0.5,-0.3,-0.1])
p1 = plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].imag, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,10,10,0,240:240+nfreq].imag, label='$d_{xz},d_{yz}$', lw='2', marker='',c='blue')
p3 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,10,10,0,240:240+nfreq].imag, label='$d_{xy}$', lw='2', marker='',c='red')
plt.legend(loc='lower right', frameon=False)

ax1 = f.add_subplot(244)
ax1.set_title(r'$R = (\pi,\pi,\pi)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.text(5,-0.10,'$Z=0.46, \gamma=0.006$', color='blue', fontsize=14)
ax1.set_xlim((0,20))
ax1.set_ylim((-0.9,0))
plt.yticks([-0.9,-0.7,-0.5,-0.3,-0.1])
p1 = plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].imag, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,10,10,10,240:240+nfreq].imag, label='$d_{xy},d_{xz},d_{yz}$', lw='2', marker='',c='blue')
plt.legend(loc='lower right', frameon=False)

ax1 = f.add_subplot(245)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{Re}[\Sigma(i\nu)]$', fontsize=14)
ax1.set_xlim((0,20))
plt.ylim(1.6,3.0)
plt.yticks([1.6,1.9,2.2,2.5,2.8])
plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].real, label=r'$\mathrm{DMFT}$', lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,0,0,240:240+nfreq].real, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='', c='blue')
plt.legend(loc=7, frameon=False)

ax1 = f.add_subplot(246)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_xlim((0,20))
plt.ylim(1.6,3.0)
plt.yticks([1.6,1.9,2.2,2.5,2.8])
plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].real, lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,10,0,240:240+nfreq].real, label='$d_{xy},d_{yz}$', lw='2', marker='', c='blue')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,0,10,0,240:240+nfreq].real, label='$d_{xz}$', lw='2', marker='', c='red')
plt.legend(loc=7, frameon=False)

ax1 = f.add_subplot(247)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_xlim((0,20))
plt.ylim(1.6,3.0)
plt.yticks([1.6,1.9,2.2,2.5,2.8])
plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].real, lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,10,10,0,240:240+nfreq].real, label='$d_{xz},d_{yz}$', lw='2', marker='', c='blue')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,10,10,0,240:240+nfreq].real, label='$d_{xy}$', lw='2', marker='', c='red')
plt.legend(loc=7, frameon=False)

ax1 = f.add_subplot(248)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_xlim((0,20))
plt.ylim(1.6,3.0)
plt.yticks([1.6,1.9,2.2,2.5,2.8])
plt.plot(fmats[:nfreq], asy['input/siw'][0,2000:2000+nfreq].real, lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,10,10,0,240:240+nfreq].real, label='$d_{xy},d_{xy},d_{xz}$', lw='2', marker='', c='blue')
plt.legend(loc=7, frameon=False)


plt.tight_layout()
plt.savefig('comparison_bugfix.eps', bbox_inches='tight', format='eps', dpi=300, rasterized=True)
plt.show()
