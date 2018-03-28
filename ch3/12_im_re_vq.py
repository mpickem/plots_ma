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

# qmc = h5py.File('../Output-SrVO3-calcgiw--bugfix--N60--qmconly/adga-20171124-162004.726-output.hdf5','r')
# asy = h5py.File('../Output-SrVO3-calcgiw--bugfix--N200/adga-20171124-025935.184-output.hdf5','r')

vq60 = h5py.File('../Output-SrVO3-readgiw--bugfix--N60-VQ-bugfix/adga-20171202-163625.597-output.hdf5','r')
# vq60g = h5py.File('../Output-SrVO3-readgiw--bugfix--N60-VQ-gamma/adga-20171213-204018.437-output.hdf5','r')
vq60g = h5py.File('../Output-SrVO3-readgiw--bugfix--N60-VQ-gamma-new/adga-20171214-121454.127-output.hdf5','r')
vq120 = h5py.File('../Output-SrVO3-readgiw--bugfix--N120-VQ-bugfix/adga-20171202-221617.720-output.hdf5','r')

vq60_hf = np.zeros((3,3,20,20,20,120),dtype=np.complex128)
vq60['selfenergy/nonloc/hartree_fock'].read_direct(vq60_hf)
vq60g_hf = np.zeros((3,3,20,20,20,120),dtype=np.complex128)
vq60g['selfenergy/nonloc/hartree_fock'].read_direct(vq60g_hf)
vq120_hf = np.zeros((3,3,20,20,20,240),dtype=np.complex128)
vq120['selfenergy/nonloc/hartree_fock'].read_direct(vq120_hf)

# fmats = np.loadtxt('siw_dmft.dat', skiprows=2, usecols=(4,), unpack=True)
nfreq = 16
freq = 60
beta = 10.0
fmats = np.arange(np.pi/beta,100,2*np.pi/beta)

## Polynomial fit to extract Z and gamma -- we use a polynomial third order for 5 frequencies
polynomial = poly.polyfit(fmats[:5], vq60g['input/siw'][0,2000:2005].imag, deg=3)
print('DMFT - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:5], vq60g['selfenergy/nonloc/dga'][0,0,0,0,0,freq:freq+5].imag, deg=3)
print('DGA-G-xy-xz-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:5], vq60g['selfenergy/nonloc/dga'][0,0,0,10,0,freq:freq+5].imag, deg=3)
print('DGA-X-xy-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:5], vq60g['selfenergy/nonloc/dga'][1,1,0,10,0,freq:freq+5].imag, deg=3)
print('DGA-X-xz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:5], vq60g['selfenergy/nonloc/dga'][0,0,10,10,0,freq:freq+5].imag, deg=3)
print('DGA-M-xy - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:5], vq60g['selfenergy/nonloc/dga'][1,1,10,10,0,freq:freq+5].imag, deg=3)
print('DGA-M-xz-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:5], vq60g['selfenergy/nonloc/dga'][1,1,10,10,10,freq:freq+5].imag, deg=3)
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
ax1.set_ylim(-0.95,-0.25)
ax1.text(0.5,-0.39,'$Z=0.56, \gamma=0.55$', color='purple', fontsize=14)
ax1.text(0.5,-0.34,'$Z=0.49, \gamma=0.36$', color='gray', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].imag, label=r'$\mathrm{DMFT}$', lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][0,0,0,0,0,freq:freq+nfreq].imag, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o', c='purple')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(242)
ax1.set_title(r'$X = (0,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.95,-0.25)
ax1.text(0.5,-0.39,'$Z=0.56, \gamma=0.56$', color='green', fontsize=14)
ax1.text(0.5,-0.34,'$Z=0.50, \gamma=0.42$', color='purple', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].imag, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][0,0,0,10,0,freq:freq+nfreq].imag, label=r'$d_{xy},d_{yz}$', lw='2', marker='o', c='purple')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][1,1,0,10,0,freq:freq+nfreq].imag, label=r'$d_{xz}$', lw='2', marker='s', c='green')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(243)
ax1.set_title(r'$M = (\pi,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.95,-0.25)
ax1.text(0.5,-0.39,'$Z=0.51, \gamma=0.41$', color='green', fontsize=14)
ax1.text(0.5,-0.34,'$Z=0.51, \gamma=0.42$', color='purple', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].imag, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][1,1,10,10,0,freq:freq+nfreq].imag, label=r'$d_{xz},d_{yz}$', lw='2', marker='o', c='purple')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][0,0,10,10,0,freq:freq+nfreq].imag, label=r'$d_{xy}$', lw='2', marker='s', c='green')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(244)
ax1.set_title(r'$R = (\pi,\pi,\pi)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(-0.95,-0.25)
ax1.text(0.5,-0.34,'$Z=0.51, \gamma=0.39$', color='purple', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].imag, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][0,0,10,10,10,freq:freq+nfreq].imag, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o', c='purple')
plt.legend(loc=4, frameon=False)


ax1 = f.add_subplot(245)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{Re}[\Sigma(i\nu)]$', fontsize=14, labelpad=10)
ax1.set_ylim(1.7,2.9)
# ax1.text(0.5,-0.46,'$Z=0.61, \gamma=0.55$', color='red', fontsize=14)
# ax1.text(0.5,-0.40,'$Z=0.49, \gamma=0.36$', color='gray', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].real, label=r'$\mathrm{DMFT}$', lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][0,0,0,0,0,freq:freq+nfreq].real + vq60g['selfenergy/nonloc/hartree_fock'][0,0,0,0,0,freq:freq+nfreq].real \
  - np.mean(vq60g['selfenergy/nonloc/hartree_fock'], axis=(2,3,4))[0,0,freq:freq+nfreq].real, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o', c='purple')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(246)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(1.7,2.9)
# ax1.text(0.5,-0.46,'$Z=0.58, \gamma=0.46$', color='red', fontsize=14)
# ax1.text(0.5,-0.40,'$Z=0.62, \gamma=0.56$', color='blue', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].real, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][0,0,0,10,0,freq:freq+nfreq].real + vq60g['selfenergy/nonloc/hartree_fock'][0,0,0,10,0,freq:freq+nfreq].real \
  - np.mean(vq60g['selfenergy/nonloc/hartree_fock'], axis=(2,3,4))[0,0,freq:freq+nfreq].real, label=r'$d_{xy},d_{yz}$', lw='2', marker='o', c='purple')
p3 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][1,1,0,10,0,freq:freq+nfreq].real + vq60g['selfenergy/nonloc/hartree_fock'][1,1,0,10,0,freq:freq+nfreq].real \
  - np.mean(vq60g['selfenergy/nonloc/hartree_fock'], axis=(2,3,4))[1,1,freq:freq+nfreq].real, label=r'$d_{xz}$', lw='2', marker='s', c='green')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(247)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(1.7,2.9)
# ax1.text(0.5,-0.46,'$Z=0.57, \gamma=0.44$', color='red', fontsize=14)
# ax1.text(0.5,-0.40,'$Z=0.58, \gamma=0.46$', color='blue', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].real, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][1,1,10,10,0,freq:freq+nfreq].real + vq60g['selfenergy/nonloc/hartree_fock'][1,1,10,10,0,freq:freq+nfreq].real \
  - np.mean(vq60g['selfenergy/nonloc/hartree_fock'], axis=(2,3,4))[1,1,freq:freq+nfreq].real, label=r'$d_{xz},d_{yz}$', lw='2', marker='o', c='purple')
p3 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][0,0,10,10,0,freq:freq+nfreq].real + vq60g['selfenergy/nonloc/hartree_fock'][0,0,10,10,0,freq:freq+nfreq].real \
  - np.mean(vq60g['selfenergy/nonloc/hartree_fock'], axis=(2,3,4))[0,0,freq:freq+nfreq].real, label=r'$d_{xy}$', lw='2', marker='s', c='green')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(248)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylim(1.7,2.9)
# ax1.text(0.5,-0.46,'$Z=0.57, \gamma=0.44$', color='red', fontsize=14)
# ax1.text(0.5,-0.40,'$Z=0.58, \gamma=0.46$', color='blue', fontsize=14)
p1 = plt.plot(fmats[:nfreq], vq60g['input/siw'][0,2000:2000+nfreq].real, lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], vq60g['selfenergy/nonloc/dga'][1,1,10,10,10,freq:freq+nfreq].real + vq60g['selfenergy/nonloc/hartree_fock'][1,1,10,10,10,freq:freq+nfreq].real \
  - np.mean(vq60g['selfenergy/nonloc/hartree_fock'], axis=(2,3,4))[1,1,freq:freq+nfreq].real, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o', c='purple')
plt.legend(loc=4, frameon=False)


plt.tight_layout()
plt.savefig('2x4_hf_gamma.eps', bbox_inches='tight', format='eps', dpi=600, rasterized=True)
# plt.show()


# g = plt.figure(figsize=(8,12))
# g.add_subplot(311)
# plt.plot(np.linspace(-60,59,120), vq60['selfenergy/nonloc/dga'][0,0,0,0,0,:].imag, label='$N=60, \Gamma$')
# plt.plot(np.linspace(-60,59,120), vq60g['selfenergy/nonloc/dga'][0,0,0,0,0,:].imag, label='$N=60, \Gamma$ with $\gamma^\omega$')
# plt.plot(np.linspace(-120,119,240), vq120['selfenergy/nonloc/dga'][0,0,0,0,0,:].imag, label='$N=120, \Gamma$')
# plt.xlim(-120,119)
# plt.legend(frameon=False)
# g.add_subplot(312)
# plt.plot(np.linspace(-60,59,120), vq60['selfenergy/nonloc/dga'][0,0,0,10,0,:].imag, label='$N=60, X$')
# plt.plot(np.linspace(-60,59,120), vq60g['selfenergy/nonloc/dga'][0,0,0,10,0,:].imag, label='$N=60, X$ with $\gamma^\omega$')
# plt.plot(np.linspace(-120,119,240), vq120['selfenergy/nonloc/dga'][0,0,0,10,0,:].imag, label='$N=120, X$')
# plt.xlim(-120,119)
# plt.legend(frameon=False)
# g.add_subplot(313)
# plt.plot(np.linspace(-60,59,120), vq60['selfenergy/nonloc/dga'][0,0,10,10,0,:].imag, label='$N=60, M$')
# plt.plot(np.linspace(-60,59,120), vq60g['selfenergy/nonloc/dga'][0,0,10,10,0,:].imag, label='$N=60, M$ with $\gamma^\omega$')
# plt.plot(np.linspace(-120,119,240), vq120['selfenergy/nonloc/dga'][0,0,10,10,0,:].imag, label='$N=120, M$')
# plt.xlim(-120,119)
# plt.legend(frameon=False)
# plt.tight_layout()
# plt.savefig('asymptotitcs.eps', bbox_inches='tight', format='eps', dpi=300, rasterized=True)


# g = plt.figure(figsize=(8,12))
# g.add_subplot(311)
# plt.plot(np.linspace(-60,59,120), vq60['selfenergy/nonloc/dga'][0,0,0,0,0,:].real, label='$N=60, \Gamma$')
# plt.plot(np.linspace(-60,59,120), vq60g['selfenergy/nonloc/dga'][0,0,0,0,0,:].real, label='$N=60, \Gamma$ with $\gamma^\omega$')
# plt.plot(np.linspace(-120,119,240), vq120['selfenergy/nonloc/dga'][0,0,0,0,0,:].real, label='$N=120, \Gamma$')
# plt.xlim(-120,119)
# plt.legend(frameon=False)
# g.add_subplot(312)
# plt.plot(np.linspace(-60,59,120), vq60['selfenergy/nonloc/dga'][0,0,0,10,0,:].real, label='$N=60, X$')
# plt.plot(np.linspace(-60,59,120), vq60g['selfenergy/nonloc/dga'][0,0,0,10,0,:].real, label='$N=60, X$ with $\gamma^\omega$')
# plt.plot(np.linspace(-120,119,240), vq120['selfenergy/nonloc/dga'][0,0,0,10,0,:].real, label='$N=120, X$')
# plt.xlim(-120,119)
# plt.legend(frameon=False)
# g.add_subplot(313)
# plt.plot(np.linspace(-60,59,120), vq60['selfenergy/nonloc/dga'][0,0,10,10,0,:].real, label='$N=60, M$')
# plt.plot(np.linspace(-60,59,120), vq60g['selfenergy/nonloc/dga'][0,0,10,10,0,:].real, label='$N=60, M$ with $\gamma^\omega$')
# plt.plot(np.linspace(-120,119,240), vq120['selfenergy/nonloc/dga'][0,0,10,10,0,:].real, label='$N=120, M$')
# plt.xlim(-120,119)
# plt.legend(frameon=False)
# plt.tight_layout()

plt.show()
