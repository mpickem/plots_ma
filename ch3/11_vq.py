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


vq_file = h5py.File('V_q_40-40-40_20-20-20.hdf5','r')
vq = vq_file['00001'][:400].real # 1111 at k_z = 0
vq_ordered = vq.reshape((20,20))

vq = vq_file['00001'][()].real # 1111 at k_z = 0
vq_odd = vq_file['00011'][()].real
vq_ordered3d = vq.reshape((20,20,20))

vq_bz = np.zeros((21,21), dtype=np.float64)
vq_bz[:20,:20] = vq_ordered[()]
vq_bz[20,:] = vq_bz[0,:]
vq_bz[:,20] = vq_bz[:,0]

bulk_plotx   = np.linspace(0,2*np.pi,22) # one more to shift the pixels
bulk_ploty   = np.linspace(0,2*np.pi,22) # one more to shift the pixels

vq_plot = np.zeros((22,22), dtype=np.float64)

vq_plot[:11,:11] = np.copy(vq_bz[10:,10:])
vq_plot[:11,10:21] = np.copy(vq_bz[10:,:11])
vq_plot[10:21,:11] = np.copy(vq_bz[:11,10:])
vq_plot[10:21,10:21] = np.copy(vq_bz[:11,:11])

# A4 - paper size = 8.3 by 11.7 inches
f = plt.figure(figsize=(8,3)) # this looks good
gs1 = gridspec.GridSpec(1,2)
gs1.update(hspace=4, wspace=0.75)

ax1 = plt.subplot(gs1[0,0])
ax1.set_title(r'$V(q_x,q_y,q_z=0)$')
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_ylabel(r'$q_y$', labelpad=-5, fontsize=14)
ax1.set_xlabel(r'$q_x$', labelpad=-1, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.pcolormesh(bulk_plotx-np.pi,bulk_ploty-np.pi,vq_plot, cmap='BuPu', rasterized=True)
a = vq_plot[:21,:21].min()
b = vq_plot[:21,:21].max()
cb1 = plt.colorbar(ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.1f', shrink=0.750)
cb1.ax.tick_params(labelsize=14)

# ax2 = plt.subplot(gs1[0,1])
ax2 = f.add_axes([0.53,0.20,0.4,0.60])
plt.plot(np.arange(11), vq_ordered3d[0,np.arange(11),0].real, lw=2, color='red')
plt.plot(np.arange(10,21), vq_ordered3d[np.arange(11),10,0].real, lw=2, color='red')
plt.plot(np.arange(20,31), vq_ordered3d[np.arange(10,-1,-1),np.arange(10,-1,-1),0].real, lw=2, color='red')
plt.plot(np.arange(30,41), vq_ordered3d[np.arange(11),np.arange(11),np.arange(11)].real, lw=2, color='red')
plt.xticks([0,10,20,30,40],[r'$\Gamma$','$X$', r'$M$', r'$\Gamma$', r'$R$'])
ax2.set_ylim(-10,150)
plt.yticks([0,40,80,120])
plt.axhline(y=0, lw=1, color='black')
ax2.set_ylabel(r'$V(\mathbf{q})$', labelpad=1, fontsize=14)



# plt.tight_layout(w_pad=0.5, h_pad=0.2)
plt.savefig('vq_kz_0.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()


print(np.sum(vq_ordered3d, axis=(0,1,2)))
