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

asy = h5py.File('../Output-u43-wo-gamma-newgiw/adga-20180125-211135.125-output.hdf5','r')

frequencies = 60
band = 0

beta = 38.0
fmats = np.linspace(np.pi/beta, 401*np.pi/beta, 201)
susc_dens = np.zeros((6,6,80,80,1,2*frequencies+1), dtype=np.complex128)
asy['susceptibility/nonloc/dens'].read_direct(susc_dens)
susc_magn = np.zeros((6,6,80,80,1,2*frequencies+1), dtype=np.complex128)
asy['susceptibility/nonloc/magn'].read_direct(susc_magn)
asy.close()

print(fmats)

dens_plt = np.zeros((6,6,82,82),dtype = np.float64)
magn_plt = np.zeros((6,6,82,82),dtype = np.float64)
dens_bz = np.zeros((6,6,81,81),dtype = np.float64)
magn_bz = np.zeros((6,6,81,81),dtype = np.float64)

dens_bz[:,:,:80,:80] = np.copy(susc_dens[...,0,frequencies].real)
magn_bz[:,:,:80,:80] = np.copy(susc_magn[...,0,frequencies].real)

# Extending the BZ to go from 0 to 2pi
dens_bz[:,:,:,80] = dens_bz[:,:,:,0]
dens_bz[:,:,80,:] = dens_bz[:,:,0,:]
magn_bz[:,:,:,80] = magn_bz[:,:,:,0]
magn_bz[:,:,80,:] = magn_bz[:,:,0,:]

plotx   = np.linspace(0,2*np.pi,82) # one more to shift the pixels
ploty   = np.linspace(0,2*np.pi,82) # one more to shift the pixels

# this is for 0 to 2pi plots
# dens_plt[:,:,:81,:81] = dens_bz[()]
# magn_plt[:,:,:81,:81] = magn_bz[()]

# this is for -pi to pi plots
dens_plt[:,:,:41,:41] = np.copy(dens_bz[:,:,40:,40:])
dens_plt[:,:,:41,40:81] = np.copy(dens_bz[:,:,40:,:41])
dens_plt[:,:,40:81,:41] = np.copy(dens_bz[:,:,:41,40:])
dens_plt[:,:,40:81,40:81] = np.copy(dens_bz[:,:,:41,:41])

magn_plt[:,:,:41,:41] = np.copy(magn_bz[:,:,40:,40:])
magn_plt[:,:,:41,40:81] = np.copy(magn_bz[:,:,40:,:41])
magn_plt[:,:,40:81,:41] = np.copy(magn_bz[:,:,:41,40:])
magn_plt[:,:,40:81,40:81] = np.copy(magn_bz[:,:,:41,:41])

# A4 - paper size = 8.3 by 11.7 inches
f = plt.figure(figsize=(16,16)) # this looks good
gs1 = gridspec.GridSpec(6,6)
gs1.update(hspace=0, wspace=0.25)


for i in xrange(6):
  for j in xrange(i,6):
    ax1 = plt.subplot(gs1[i,j])
    ax1.tick_params(axis='both', which='major', labelsize=14)
    if (i==j):
      ax1.set_ylabel(r'$q_y$', labelpad=-5, fontsize=14)
      ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
    else:
      plt.xticks(visible=False)
      plt.yticks(visible=False)
    if (i==0):
      if (j==0):
        ax1.set_title(r'$xy\;(\mathrm{top})$')
      elif (j==1):
        ax1.set_title(r'$xz\;(\mathrm{top})$')
      elif (j==2):
        ax1.set_title(r'$yz\;(\mathrm{top})$')
      elif (j==3):
        ax1.set_title(r'$xy\;(\mathrm{bottom})$')
      elif (j==4):
        ax1.set_title(r'$xz\;(\mathrm{bottom})$')
      elif (j==5):
        ax1.set_title(r'$yz\;(\mathrm{bottom})$')
    ax1.set_xlim(-np.pi,np.pi)
    ax1.set_ylim(-np.pi,np.pi)
    ax1.set_aspect('equal')
    plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
    plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
    plt.pcolormesh(plotx-np.pi,ploty-np.pi,dens_plt[i,j,...], cmap='magma', rasterized=True)
    a = dens_plt[i,j,:-1,:-1].min()
    b = dens_plt[i,j,:-1,:-1].max()
    cb1 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.2f', shrink=0.650)
    cb1.ax.tick_params(labelsize=14)

# ax = f.add_axes([0,0.5,0.5,0.5])
# ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
# ax.patch.set_alpha(0.2)
# ax.patch.set_color('red')

# ax = f.add_axes([0.5,0,0.5,0.5])
# ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
# ax.patch.set_alpha(0.2)
# ax.patch.set_color('blue')

# ax = f.add_axes([0.5,0.5,0.5,0.5])
# ax.xaxis.set_visible(False)
# ax.yaxis.set_visible(False)
# ax.patch.set_alpha(0.2)
# ax.patch.set_color('black')

# plt.tight_layout(w_pad=0.5, h_pad=0.2)
plt.savefig('layer-cond-3x6-susc-dens.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()
