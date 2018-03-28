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

bulk = h5py.File('output-srvo3-beta38-n240/adga-20171209-114252.312-output.hdf5','r')
layer = h5py.File('/home/lv70585/pickemERC/MA_ADGA_runs/Output-u43-wo-gamma-newgiw/adga-20180125-211135.125-output.hdf5','r')

bulk_f = 240
layer_f = 60

beta = 38.0
fmats = np.linspace(np.pi/beta, 401*np.pi/beta, 201)
bulk_dens = np.zeros((3,3,20,20,20,2*bulk_f+1), dtype=np.complex128)
bulk['susceptibility/nonloc/dens'].read_direct(bulk_dens)
bulk_magn = np.zeros((3,3,20,20,20,2*bulk_f+1), dtype=np.complex128)
bulk['susceptibility/nonloc/magn'].read_direct(bulk_magn)
bulk.close()

layer_dens = np.zeros((6,6,80,80,1,2*layer_f+1), dtype=np.complex128)
layer['susceptibility/nonloc/dens'].read_direct(layer_dens)
layer_magn = np.zeros((6,6,80,80,1,2*layer_f+1), dtype=np.complex128)
layer['susceptibility/nonloc/magn'].read_direct(layer_magn)
layer.close()

#  BULK
bulk_dens_plt = np.zeros((3,3,22,22),dtype = np.float64)
bulk_magn_plt = np.zeros((3,3,22,22),dtype = np.float64)
bulk_dens_bz = np.zeros((3,3,21,21),dtype = np.float64)
bulk_magn_bz = np.zeros((3,3,21,21),dtype = np.float64)

bulk_dens_bz[:,:,:20,:20] = np.copy(bulk_dens[...,0,bulk_f].real)
bulk_magn_bz[:,:,:20,:20] = np.copy(bulk_magn[...,0,bulk_f].real)

# Extending the BZ to go from 0 to 2pi
bulk_dens_bz[:,:,:,20] = bulk_dens_bz[:,:,:,0]
bulk_dens_bz[:,:,20,:] = bulk_dens_bz[:,:,0,:]
bulk_magn_bz[:,:,:,20] = bulk_magn_bz[:,:,:,0]
bulk_magn_bz[:,:,20,:] = bulk_magn_bz[:,:,0,:]

bulk_plotx   = np.linspace(0,2*np.pi,22) # one more to shift the pixels
bulk_ploty   = np.linspace(0,2*np.pi,22) # one more to shift the pixels

bulk_dens_plt[:,:,:11,:11] = np.copy(bulk_dens_bz[:,:,10:,10:])
bulk_dens_plt[:,:,:11,10:21] = np.copy(bulk_dens_bz[:,:,10:,:11])
bulk_dens_plt[:,:,10:21,:11] = np.copy(bulk_dens_bz[:,:,:11,10:])
bulk_dens_plt[:,:,10:21,10:21] = np.copy(bulk_dens_bz[:,:,:11,:11])

bulk_magn_plt[:,:,:11,:11] = np.copy(bulk_magn_bz[:,:,10:,10:])
bulk_magn_plt[:,:,:11,10:21] = np.copy(bulk_magn_bz[:,:,10:,:11])
bulk_magn_plt[:,:,10:21,:11] = np.copy(bulk_magn_bz[:,:,:11,10:])
bulk_magn_plt[:,:,10:21,10:21] = np.copy(bulk_magn_bz[:,:,:11,:11])

#  LAYER
layer_dens_plt = np.zeros((6,6,82,82),dtype = np.float64)
layer_magn_plt = np.zeros((6,6,82,82),dtype = np.float64)
layer_dens_bz = np.zeros((6,6,81,81),dtype = np.float64)
layer_magn_bz = np.zeros((6,6,81,81),dtype = np.float64)

layer_dens_bz[:,:,:80,:80] = np.copy(layer_dens[...,0,layer_f].real)
layer_magn_bz[:,:,:80,:80] = np.copy(layer_magn[...,0,layer_f].real)

# Extending the BZ to go from 0 to 2pi
layer_dens_bz[:,:,:,80] = layer_dens_bz[:,:,:,0]
layer_dens_bz[:,:,80,:] = layer_dens_bz[:,:,0,:]
layer_magn_bz[:,:,:,80] = layer_magn_bz[:,:,:,0]
layer_magn_bz[:,:,80,:] = layer_magn_bz[:,:,0,:]

layer_plotx   = np.linspace(0,2*np.pi,82) # one more to shift the pixels
layer_ploty   = np.linspace(0,2*np.pi,82) # one more to shift the pixels

layer_dens_plt[:,:,:41,:41] = np.copy(layer_dens_bz[:,:,40:,40:])
layer_dens_plt[:,:,:41,40:81] = np.copy(layer_dens_bz[:,:,40:,:41])
layer_dens_plt[:,:,40:81,:41] = np.copy(layer_dens_bz[:,:,:41,40:])
layer_dens_plt[:,:,40:81,40:81] = np.copy(layer_dens_bz[:,:,:41,:41])

layer_magn_plt[:,:,:41,:41] = np.copy(layer_magn_bz[:,:,40:,40:])
layer_magn_plt[:,:,:41,40:81] = np.copy(layer_magn_bz[:,:,40:,:41])
layer_magn_plt[:,:,40:81,:41] = np.copy(layer_magn_bz[:,:,:41,40:])
layer_magn_plt[:,:,40:81,40:81] = np.copy(layer_magn_bz[:,:,:41,:41])

# A4 - paper size = 8.3 by 11.7 inches
f = plt.figure(figsize=(14,6)) # this looks good
gs1 = gridspec.GridSpec(2,4,hspace=0.0,wspace=0.4)
# gs1.update(hspace=0.5, wspace=0.75)

ax1 = plt.subplot(gs1[0,0])
ax1.set_title(r'$\mathrm{bulk}$', fontsize=14)
# ax1.set_title(r'$\mathrm{(a)}\;\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_ylabel(r'$q_y$', labelpad=-5, fontsize=14)
# ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(bulk_plotx-np.pi,bulk_ploty-np.pi,np.sum(bulk_dens_plt[()], axis=(0,1)), cmap='magma', rasterized=True)
a = np.sum(bulk_dens_plt[:,:,:-1,:-1], axis=(0,1)).min()
b = np.sum(bulk_dens_plt[:,:,:-1,:-1], axis=(0,1)).max()
cb1 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb1.ax.tick_params(labelsize=14)

ax1 = plt.subplot(gs1[0,1])
ax1.set_title(r'$\mathrm{top\;layer}$', fontsize=14)
# ax1.set_title(r'$\mathrm{(a)}\;\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
# ax1.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
# ax1.set_xlabel(r'$k_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(layer_plotx-np.pi,layer_ploty-np.pi,np.sum(layer_dens_plt[:3,:3,...], axis=(0,1)), cmap='magma', rasterized=True)
a = np.sum(layer_dens_plt[:3,:3,:-1,:-1], axis=(0,1)).min()
b = np.sum(layer_dens_plt[:3,:3,:-1,:-1], axis=(0,1)).max()
cb2 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb2.ax.tick_params(labelsize=14)

ax1 = plt.subplot(gs1[0,2])
ax1.set_title(r'$\mathrm{bottom\;layer}$', fontsize=14)
# ax1.set_title(r'$\mathrm{(a)}\;\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
# ax1.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
# ax1.set_xlabel(r'$k_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(layer_plotx-np.pi,layer_ploty-np.pi,np.sum(layer_dens_plt[3:,3:,...], axis=(0,1)), cmap='magma', rasterized=True)
a = np.sum(layer_dens_plt[3:,3:,:-1,:-1], axis=(0,1)).min()
b = np.sum(layer_dens_plt[3:,3:,:-1,:-1], axis=(0,1)).max()
cb3 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb3.ax.tick_params(labelsize=14)

ax1 = plt.subplot(gs1[0,3])
ax1.set_title(r'$\mathrm{off\;impurity}$', fontsize=14)
# ax1.set_title(r'$\mathrm{(a)}\;\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
# ax1.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
# ax1.set_xlabel(r'$k_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(layer_plotx-np.pi,layer_ploty-np.pi,np.sum(layer_dens_plt[3:,:3,...], axis=(0,1)), cmap='magma', rasterized=True)
a = np.sum(layer_dens_plt[3:,:3,:-1,:-1], axis=(0,1)).min()
b = np.sum(layer_dens_plt[3:,:3,:-1,:-1], axis=(0,1)).max()
cb4 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb4.ax.tick_params(labelsize=14)

ax1 = plt.subplot(gs1[1,0])
# ax1.set_title(r'$\mathrm{(b)}\;\mathrm{Re}[\chi_M(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_ylabel(r'$q_y$', labelpad=-5, fontsize=14)
ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(bulk_plotx-np.pi,bulk_ploty-np.pi,2*np.sum(bulk_magn_plt[()], axis=(0,1)), cmap='magma', rasterized=True)
a = 2*np.sum(bulk_magn_plt[:,:,:-1,:-1], axis=(0,1)).min()
b = 2*np.sum(bulk_magn_plt[:,:,:-1,:-1], axis=(0,1)).max()
cb5 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb5.ax.tick_params(labelsize=14)

ax1 = plt.subplot(gs1[1,1])
# ax1.set_title(r'$\mathrm{(a)}\;\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
# ax1.set_ylabel(r'$q_y$', labelpad=-5, fontsize=14)
ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(layer_plotx-np.pi,layer_ploty-np.pi,2*np.sum(layer_magn_plt[:3,:3,...], axis=(0,1)), cmap='magma', rasterized=True)
a = 2*np.sum(layer_magn_plt[:3,:3,:-1,:-1], axis=(0,1)).min()
b = 2*np.sum(layer_magn_plt[:3,:3,:-1,:-1], axis=(0,1)).max()
cb6 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb6.ax.tick_params(labelsize=14)

ax1 = plt.subplot(gs1[1,2])
# ax1.set_title(r'$\mathrm{bottom\;layer}$')
# ax1.set_title(r'$\mathrm{(a)}\;\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
# ax1.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(layer_plotx-np.pi,layer_ploty-np.pi,2*np.sum(layer_magn_plt[3:,3:,...], axis=(0,1)), cmap='magma', rasterized=True)
a = 2*np.sum(layer_magn_plt[3:,3:,:-1,:-1], axis=(0,1)).min()
b = 2*np.sum(layer_magn_plt[3:,3:,:-1,:-1], axis=(0,1)).max()
cb7 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb7.ax.tick_params(labelsize=14)

ax1 = plt.subplot(gs1[1,3])
# ax1.set_title(r'$\mathrm{off\;impurity}$')
# ax1.set_title(r'$\mathrm{(a)}\;\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$')
ax1.tick_params(axis='both', which='major', labelsize=14)
# ax1.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(layer_plotx-np.pi,layer_ploty-np.pi,2*np.sum(layer_magn_plt[3:,:3,...], axis=(0,1)), cmap='magma', rasterized=True)
a = 2*np.sum(layer_magn_plt[3:,:3,:-1,:-1], axis=(0,1)).min()
b = 2*np.sum(layer_magn_plt[3:,:3,:-1,:-1], axis=(0,1)).max()
cb8 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.705)
cb8.ax.tick_params(labelsize=14)



plt.gcf().text(0.93,0.41, r'$\mathrm{Re}[\chi_M(\mathbf{q},\omega=0)]$', fontsize=18, rotation=90)
plt.gcf().text(0.93,0.81, r'$\mathrm{Re}[\chi_D(\mathbf{q},\omega=0)]$', fontsize=18, rotation=90)


# plt.tight_layout(w_pad=0.5, h_pad=0.2)
plt.savefig('susc_comp_bulk_layer.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()
