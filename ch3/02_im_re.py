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

qmc = h5py.File('../Output-SrVO3-calcgiw--bugfix--N60--qmconly/adga-20171124-162004.726-output.hdf5','r')
asy = h5py.File('../Output-SrVO3-calcgiw--bugfix--N200/adga-20171124-025935.184-output.hdf5','r')

# fmats = np.loadtxt('siw_dmft.dat', skiprows=2, usecols=(4,), unpack=True)
# use genfromtxt ... more sophisticated
fmats, siw_dmft_r, siw_dmft_i = np.genfromtxt('siw_dmft.dat', skip_header=2, usecols=(4,5,6), unpack=True)

nfreq = 16

# gamma -- xy = xz = yz
siw_old_g_r, siw_old_g_i = np.genfromtxt('siwk_old_gamma.dat', usecols=(1,2), unpack=True)
siw_new_g_r, siw_new_g_i = np.genfromtxt('siwk_new_gamma.dat', usecols=(1,2), unpack=True)

# X -- xy = yz
siw_old_x_0_r, siw_old_x_0_i = np.genfromtxt('siwk_old_X.dat', usecols=(1,2), unpack=True)
siw_new_x_0_r, siw_new_x_0_i = np.genfromtxt('siwk_new_X.dat', usecols=(1,2), unpack=True)
# X -- xz
siw_old_x_1_r, siw_old_x_1_i = np.genfromtxt('siwk_old_X.dat', usecols=(7,8), unpack=True)
siw_new_x_1_r, siw_new_x_1_i = np.genfromtxt('siwk_new_X.dat', usecols=(7,8), unpack=True)

# M -- xy
siw_old_m_0_r, siw_old_m_0_i = np.genfromtxt('siwk_old_M.dat', usecols=(1,2), unpack=True)
siw_new_m_0_r, siw_new_m_0_i = np.genfromtxt('siwk_new_M.dat', usecols=(1,2), unpack=True)
# M -- xz=yz
siw_old_m_1_r, siw_old_m_1_i = np.genfromtxt('siwk_old_M.dat', usecols=(7,8), unpack=True)
siw_new_m_1_r, siw_new_m_1_i = np.genfromtxt('siwk_new_M.dat', usecols=(7,8), unpack=True)

freq, aw_g_dmft = np.genfromtxt('svo-spectra/dmft/aw_gkiw_dmft_0-0-0_0.dat', usecols=(0,1),skip_header=1,  unpack=True)
aw_x_0_dmft = np.genfromtxt('svo-spectra/dmft/aw_gkiw_dmft_0-10-0_0.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_x_1_dmft = np.genfromtxt('svo-spectra/dmft/aw_gkiw_dmft_0-10-0_1.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_m_0_dmft = np.genfromtxt('svo-spectra/dmft/aw_gkiw_dmft_10-10-0_0.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_m_1_dmft = np.genfromtxt('svo-spectra/dmft/aw_gkiw_dmft_10-10-0_1.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_r_0_dmft = np.genfromtxt('svo-spectra/dmft/aw_gkiw_dmft_10-10-10_0.dat', usecols=(1,),skip_header=1, unpack=True)

aw_g_dga = np.genfromtxt('svo-spectra/dga/aw_gkiw_dga_0-0-0_0.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_x_0_dga = np.genfromtxt('svo-spectra/dga/aw_gkiw_dga_0-10-0_0.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_x_1_dga = np.genfromtxt('svo-spectra/dga/aw_gkiw_dga_0-10-0_1.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_m_0_dga = np.genfromtxt('svo-spectra/dga/aw_gkiw_dga_10-10-0_0.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_m_1_dga = np.genfromtxt('svo-spectra/dga/aw_gkiw_dga_10-10-0_1.dat', usecols=(1,),skip_header=1,  unpack=True)
aw_r_0_dga = np.genfromtxt('svo-spectra/dga/aw_gkiw_dga_10-10-10_0.dat', usecols=(1,), skip_header=1, unpack=True)

## Polynomial fit to extract Z and gamma -- we use a polynomial third order for 5 frequencies
polynomial = poly.polyfit(fmats[:5], siw_dmft_i[:5], deg=3)
print('DMFT - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
# polynomial = poly.polyfit(fmats[:5], siw_new_g_i[:5], deg=3)
polynomial = poly.polyfit(fmats[:5], asy['selfenergy/nonloc/dga'][0,0,0,0,0,200:205].imag, deg=3)
print('DGA-G-xy-xz-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
# polynomial = poly.polyfit(fmats[:5], siw_new_x_0_i[:5], deg=3)
polynomial = poly.polyfit(fmats[:5], asy['selfenergy/nonloc/dga'][0,0,0,10,0,200:205].imag, deg=3)
print('DGA-X-xy-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
# polynomial = poly.polyfit(fmats[:5], siw_new_x_1_i[:5], deg=3)
polynomial = poly.polyfit(fmats[:5], asy['selfenergy/nonloc/dga'][1,1,0,10,0,200:205].imag, deg=3)
print('DGA-X-xz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
# polynomial = poly.polyfit(fmats[:5], siw_new_m_0_i[:5], deg=3)
polynomial = poly.polyfit(fmats[:5], asy['selfenergy/nonloc/dga'][0,0,10,10,0,200:205].imag, deg=3)
print('DGA-M-xy - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
# polynomial = poly.polyfit(fmats[:5], siw_new_m_1_i[:5], deg=3)
polynomial = poly.polyfit(fmats[:5], asy['selfenergy/nonloc/dga'][1,1,10,10,0,200:205].imag, deg=3)
print('DGA-M-xz-yz - gamma - Z: ', -polynomial[0], (1-polynomial[1])**(-1))
polynomial = poly.polyfit(fmats[:5], asy['selfenergy/nonloc/dga'][1,1,10,10,10,200:205].imag, deg=3)
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
ax1.text(0.5,-0.41,'$Z=0.49, \gamma=0.38$', color='blue', fontsize=14)
ax1.text(0.5,-0.35,'$Z=0.49, \gamma=0.36$', color='gray', fontsize=14)
plt.ylim(-0.95,-0.25)
plt.yticks([-0.9,-0.7,-0.5,-0.3])
# ax1.axes().get_xaxis().set_visible(False)
# plt.xticks([])
p1 = plt.plot(fmats[:nfreq], siw_dmft_i[:nfreq], label=r'$\mathrm{DMFT}$', lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,0,0,200:216].imag, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o', c='blue')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(242)
ax1.set_title(r'$X = (0,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.text(0.5,-0.41,'$Z=0.50, \gamma=0.38$', color='red', fontsize=14)
ax1.text(0.5,-0.35,'$Z=0.50, \gamma=0.40$', color='blue', fontsize=14)
plt.yticks([-0.9,-0.7,-0.5,-0.3])
plt.ylim(-0.95,-0.25)
p1 = plt.plot(fmats[:nfreq], siw_dmft_i[:nfreq], lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,10,0,200:216].imag, label='$d_{xy},d_{yz}$', lw='2', marker='o',c='blue')
p3 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,0,10,0,200:216].imag, label='$d_{xz}$', lw='2', marker='s',c='red')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(243)
ax1.set_title(r'$M = (\pi,\pi,0)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.text(0.5,-0.41,'$Z=0.50, \gamma=0.42$', color='red', fontsize=14)
ax1.text(0.5,-0.35,'$Z=0.50, \gamma=0.40$', color='blue', fontsize=14)
plt.yticks([-0.9,-0.7,-0.5,-0.3])
plt.ylim(-0.95,-0.25)
p1 = plt.plot(fmats[:nfreq], siw_dmft_i[:nfreq], lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,10,10,0,200:216].imag, label='$d_{xz},d_{yz}$', lw='2', marker='o',c='blue')
p3 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,10,10,0,200:216].imag, label='$d_{xy}$', lw='2', marker='s',c='red')
plt.legend(loc=4, frameon=False)

ax1 = f.add_subplot(244)
ax1.set_title(r'$R = (\pi,\pi,\pi)$', fontsize=14)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.text(0.5,-0.35,'$Z=0.50, \gamma=0.42$', color='blue', fontsize=14)
plt.yticks([-0.9,-0.7,-0.5,-0.3])
plt.ylim(-0.95,-0.25)
p1 = plt.plot(fmats[:nfreq], siw_dmft_i[:nfreq], lw='4', ls='--', c='gray')
p2 = plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,10,10,10,200:216].imag, label='$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o',c='blue')
plt.legend(loc=4, frameon=False)


# next row

ax1 = f.add_subplot(245)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
ax1.set_ylabel(r'$\mathrm{Re}[\Sigma(i\nu)]$', fontsize=14, labelpad=10)
plt.ylim(1.9,2.9)
plt.plot(fmats[:nfreq], siw_dmft_r[:nfreq], label=r'$\mathrm{DMFT}$', lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,0,0,200:216].real, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o', c='blue')
plt.legend(loc=7, frameon=False)

ax1 = f.add_subplot(246)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
plt.ylim(1.9,2.9)
plt.plot(fmats[:nfreq], siw_dmft_r[:nfreq], lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,0,10,0,200:216].real, label='$d_{xy},d_{yz}$', lw='2', marker='o', c='blue')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,0,10,0,200:216].real, label='$d_{xz}$', lw='2', marker='s', c='red')
plt.legend(loc=7, frameon=False)

ax1 = f.add_subplot(247)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
plt.ylim(1.9,2.9)
plt.plot(fmats[:nfreq], siw_dmft_r[:nfreq], lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][1,1,10,10,0,200:216].real, label='$d_{xz},d_{yz}$', lw='2', marker='o', c='blue')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,10,10,0,200:216].real, label='$d_{xy}$', lw='2', marker='s', c='red')
plt.legend(loc=7, frameon=False)

ax1 = f.add_subplot(248)
ax1.set_xlabel(r'$i\nu$', fontsize=14)
plt.ylim(1.9,2.9)
plt.plot(fmats[:nfreq], siw_dmft_r[:nfreq], lw='4', ls='--', c='gray')
plt.plot(fmats[:nfreq], asy['selfenergy/nonloc/dga'][0,0,10,10,10,200:216].real, label=r'$d_{xy},d_{xz},d_{yz}$', lw='2', marker='o', c='blue')
plt.legend(loc=7, frameon=False)

# next row

# ax1 = f.add_subplot(349)
# ax1.set_xlabel(r'$\nu$', fontsize=14)
# ax1.set_ylabel(r'$A(\nu)$', fontsize=14, labelpad=10)
# plt.xlim(-3,4)
# plt.ylim(0,1.6)
# plt.yticks([0,0.4,0.8,1.2,1.6])
# plt.plot(freq, aw_g_dmft, label=r'$\mathrm{DMFT}$', lw='4', ls='--', c='gray')
# plt.plot(freq, aw_g_dga, label='$d_{xy},d_{xz},d_{yz}$', lw='2', ls='-', c='blue')
# plt.legend(loc=1, frameon=False)

# ax1 = f.add_subplot(3,4,10)
# ax1.set_xlabel(r'$\nu$', fontsize=14)
# plt.xlim(-3,4)
# plt.ylim(0,1.6)
# plt.yticks([0,0.4,0.8,1.2,1.6])
# plt.plot(freq,aw_x_0_dmft, lw='4', ls='--', c='gray')
# plt.plot(freq,aw_x_1_dmft, lw='4', ls='--', c='gray')
# plt.plot(freq,aw_x_0_dga, label='$d_{xy},d_{yz}$', lw='2', ls='-', c='blue')
# plt.plot(freq,aw_x_1_dga, label='$d_{xz}$', lw='2', ls='-', c='red')
# plt.legend(loc=1, frameon=False)

# ax1 = f.add_subplot(3,4,11)
# plt.xlim(-3,4)
# plt.ylim(0,1.6)
# plt.yticks([0,0.4,0.8,1.2,1.6])
# ax1.set_xlabel(r'$\nu$', fontsize=14)
# plt.plot(freq,aw_m_0_dmft, lw='4', ls='--', c='gray')
# plt.plot(freq,aw_m_1_dmft, lw='4', ls='--', c='gray')
# plt.plot(freq,aw_m_1_dga, label='$d_{xz},d_{yz}$', lw='2', ls='-', c='blue')
# plt.plot(freq,aw_m_0_dga, label='$d_{xy}$', lw='2', ls='-', c='red')
# plt.legend(loc=1, frameon=False)

# ax1 = f.add_subplot(3,4,12)
# plt.xlim(-3,4)
# plt.ylim(0,1.6)
# plt.yticks([0,0.4,0.8,1.2,1.6])
# ax1.set_xlabel(r'$\nu$', fontsize=14)
# plt.plot(freq,aw_r_0_dmft, lw='4', ls='--', c='gray')
# plt.plot(freq,aw_r_0_dga, label='$d_{xy},d_{xz},d_{yz}$', lw='2', ls='-', c='blue')
# plt.legend(loc=1, frameon=False)


plt.tight_layout()
plt.savefig('2x4.eps', bbox_inches='tight', format='eps', dpi=600, rasterized=True)
plt.show()
