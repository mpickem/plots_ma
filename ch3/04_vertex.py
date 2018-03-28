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

asy = h5py.File('out-plot.hdf5','r')

beta = 38.0
fmats = np.linspace(-119*np.pi/beta, 121*np.pi/beta, 121)

# for symmetric interpolation and plotting
qmc_d_plt = np.zeros((121,121),dtype = np.float64)
qmc_m_plt = np.zeros((121,121),dtype = np.float64)
asy_d_plt = np.zeros((121,121),dtype = np.float64)
asy_m_plt = np.zeros((121,121),dtype = np.float64)
mix_d_plt = np.zeros((121,121),dtype = np.float64)
mix_m_plt = np.zeros((121,121),dtype = np.float64)

qmc_d_plt[:120,:120] = asy['qmc_dens']
qmc_m_plt[:120,:120] = asy['qmc_magn']
asy_d_plt[:120,:120] = asy['asy_dens']
asy_m_plt[:120,:120] = asy['asy_magn']
mix_d_plt[:120,:120] = asy['mix_dens']
mix_m_plt[:120,:120] = asy['mix_magn']


f = plt.figure(figsize=(8,5)) # this looks good

a = asy_d_plt[:-1,:-1].min()
b = asy_d_plt[:-1,:-1].max()

ax1 = f.add_subplot(231)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_title(r"$\mathrm{(a)\enspace} F_{\mathrm{FULL}}^{\mathrm{dens}}(i\nu,i\nu',i\omega=0)$")
ax1.set_ylabel(r"$i\nu'$", labelpad=-5, fontsize=14)
ax1.set_xlim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_ylim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_aspect('equal')
plt.xticks([-8,-4,0,4,8])
plt.yticks([-8,-4,0,4,8])
plt.pcolormesh(fmats,fmats,qmc_d_plt, cmap='inferno', rasterized=True, vmin=a, vmax=b)

ax1 = f.add_subplot(232)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_title(r"$\mathrm{(b)\enspace} F_{\mathrm{ASY}}^{\mathrm{dens}}(i\nu,i\nu',i\omega=0)$")
ax1.set_xlim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_ylim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_aspect('equal')
plt.xticks([-8,-4,0,4,8])
plt.yticks([-8,-4,0,4,8])
plt.pcolormesh(fmats,fmats,asy_d_plt, cmap='inferno', rasterized=True, vmin=a, vmax=b)

ax1 = f.add_subplot(233)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_title(r"$\mathrm{(c)\enspace} F_{\mathrm{REP}}^{\mathrm{dens}}(i\nu,i\nu',i\omega=0)$")
ax1.set_xlim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_ylim(-119*np.pi/beta,119*np.pi/beta)
# ax1.set_aspect('equal')
plt.xticks([-8,-4,0,4,8])
plt.yticks([-8,-4,0,4,8])
plt.pcolormesh(fmats,fmats,mix_d_plt, cmap='inferno', rasterized=True, vmin=a, vmax=b)
cb1 = plt.colorbar(ticks=[a,a+(b-a)/5, a+(b-a)*2/5, a+(b-a)*3/5, a+(b-a)*4/5, b], format='%.1f', shrink=0.995)
cb1.ax.tick_params(labelsize=14)

a = asy_m_plt[:-1,:-1].min()
b = asy_m_plt[:-1,:-1].max()

ax1 = f.add_subplot(234)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_title(r"$\mathrm{(d)\enspace} F_{\mathrm{FULL}}^{\mathrm{magn}}(i\nu,i\nu',i\omega=0)$")
ax1.set_xlabel(r'$i\nu$', labelpad=2, fontsize=14)
ax1.set_ylabel(r"$i\nu'$", labelpad=-5, fontsize=14)
ax1.set_xlim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_ylim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_aspect('equal')
plt.xticks([-8,-4,0,4,8])
plt.yticks([-8,-4,0,4,8])
plt.pcolormesh(fmats,fmats,qmc_m_plt, cmap='inferno', rasterized=True, vmin=a, vmax=b)

ax1 = f.add_subplot(235)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_title(r"$\mathrm{(e)\enspace} F_{\mathrm{ASY}}^{\mathrm{magn}}(i\nu,i\nu',i\omega=0)$")
ax1.set_xlabel(r'$i\nu$', labelpad=2, fontsize=14)
ax1.set_xlim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_ylim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_aspect('equal')
plt.xticks([-8,-4,0,4,8])
plt.yticks([-8,-4,0,4,8])
plt.pcolormesh(fmats,fmats,asy_m_plt, cmap='inferno', rasterized=True, vmin=a, vmax=b)

ax1 = f.add_subplot(236)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_title(r"$\mathrm{(f)\enspace} F_{\mathrm{REP}}^{\mathrm{magn}}(i\nu,i\nu',i\omega=0)$")
ax1.set_xlabel(r'$i\nu$', labelpad=2, fontsize=14)
ax1.set_xlim(-119*np.pi/beta,119*np.pi/beta)
ax1.set_ylim(-119*np.pi/beta,119*np.pi/beta)
# ax1.set_aspect('equal')
plt.xticks([-8,-4,0,4,8])
plt.yticks([-8,-4,0,4,8])
plt.pcolormesh(fmats,fmats,mix_m_plt, cmap='inferno', rasterized=True, vmin=a, vmax=b)
cb1 = plt.colorbar(ticks=[a,a+(b-a)/5, a+(b-a)*2/5, a+(b-a)*3/5, a+(b-a)*4/5, b], format='%.1f', shrink=0.995)
cb1.ax.tick_params(labelsize=14)

plt.tight_layout()
plt.savefig('qmc_asy_mix.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()
