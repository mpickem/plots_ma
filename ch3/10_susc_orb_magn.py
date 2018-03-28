#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import h5py
import sys

# make matplot math look like latex
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

asy = h5py.File('output-srvo3-beta38-n240/adga-20171209-114252.312-output.hdf5','r')

frequencies = 240
band = 0

beta = 38.0
fmats = np.linspace(np.pi/beta, 401*np.pi/beta, 201)
susc_dens = np.zeros((3,3,20,20,20,2*frequencies+1), dtype=np.complex128)
asy['susceptibility/nonloc/dens'].read_direct(susc_dens)
susc_magn = np.zeros((3,3,20,20,20,2*frequencies+1), dtype=np.complex128)
asy['susceptibility/nonloc/magn'].read_direct(susc_magn)
asy.close()

print(fmats)

dens_plt = np.zeros((3,3,22,22),dtype = np.float64)
magn_plt = np.zeros((3,3,22,22),dtype = np.float64)
dens_bz = np.zeros((3,3,21,21),dtype = np.float64)
magn_bz = np.zeros((3,3,21,21),dtype = np.float64)

dens_bz[:,:,:20,:20] = np.copy(susc_dens[...,0,frequencies].real)
magn_bz[:,:,:20,:20] = np.copy(susc_magn[...,0,frequencies].real)

# Extending the BZ to go from 0 to 2pi
dens_bz[:,:,:,20] = dens_bz[:,:,:,0]
dens_bz[:,:,20,:] = dens_bz[:,:,0,:]
magn_bz[:,:,:,20] = magn_bz[:,:,:,0]
magn_bz[:,:,20,:] = magn_bz[:,:,0,:]

plotx   = np.linspace(0,2*np.pi,22) # one more to shift the pixels
ploty   = np.linspace(0,2*np.pi,22) # one more to shift the pixels

dens_plt[:,:,:11,:11] = np.copy(dens_bz[:,:,10:,10:])
dens_plt[:,:,:11,10:21] = np.copy(dens_bz[:,:,10:,:11])
dens_plt[:,:,10:21,:11] = np.copy(dens_bz[:,:,:11,10:])
dens_plt[:,:,10:21,10:21] = np.copy(dens_bz[:,:,:11,:11])

magn_plt[:,:,:11,:11] = np.copy(magn_bz[:,:,10:,10:])
magn_plt[:,:,:11,10:21] = np.copy(magn_bz[:,:,10:,:11])
magn_plt[:,:,10:21,:11] = np.copy(magn_bz[:,:,:11,10:])
magn_plt[:,:,10:21,10:21] = np.copy(magn_bz[:,:,:11,:11])

# A4 - paper size = 8.3 by 11.7 inches
f = plt.figure(figsize=(8,8)) # this looks good
gs1 = gridspec.GridSpec(3,3)
gs1.update(hspace=0, wspace=0.25)

for i in xrange(3):
  for j in xrange(3):
    ax1 = plt.subplot(gs1[i,j])
    ax1.tick_params(axis='both', which='major', labelsize=14)
    if i==0 and j==0:
      ax1.set_title(r'$xy$', fontsize=14)
    if i==0 and j==1:
      ax1.set_title(r'$xz$', fontsize=14)
    if i==0 and j==2:
      ax1.set_title(r'$yz$', fontsize=14)
    if j==0:
      ax1.set_ylabel(r'$q_y$', labelpad=-5, fontsize=14)
    else:
      plt.yticks(visible=False)
    if i==2:
      ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
    else:
      plt.xticks(visible=False)
    ax1.set_xlim(-np.pi,np.pi)
    ax1.set_ylim(-np.pi,np.pi)
    ax1.set_aspect('equal')
    plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
    plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
    plt.pcolormesh(plotx-np.pi,ploty-np.pi,magn_plt[i,j,...], cmap='magma', rasterized=True)
    a = magn_plt[i,j,:-1,:-1].min()
    b = magn_plt[i,j,:-1,:-1].max()
    cb1 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.2f', shrink=0.650)
    cb1.ax.tick_params(labelsize=14)


# plt.tight_layout(w_pad=0.5, h_pad=0.2)
plt.savefig('srvo3-b38-susc-magn.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()
