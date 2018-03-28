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

band = 0
frequencies = 240


# asy = h5py.File('output-srvo3-beta38-n60/adga-20171207-114802.202-output.hdf5','r')
# asy = h5py.File('output-srvo3-beta38-n60-gamma/adga-20171213-174024.449-output.hdf5','r')
# asy = h5py.File('output-srvo3-beta38-n120/adga-20171207-215254.224-output.hdf5','r')
asy = h5py.File('output-srvo3-beta38-n240/adga-20171209-114252.312-output.hdf5','r')
beta = 38.0
fmats = np.linspace(np.pi/beta, 1001*np.pi/beta, 501)
siwdga = np.zeros((3,3,20,20,20,2*frequencies), dtype=np.complex128)
siwdmft = np.zeros((3,4000), dtype=np.complex128)
hk = np.zeros((3,3,20,20,20), dtype=np.complex128)
asy['selfenergy/nonloc/dga'].read_direct(siwdga)
asy['input/siw'].read_direct(siwdmft)
asy['input/hk'].read_direct(hk)
asy.close()



print(fmats[:50])

# x,y = np.argwhere(np.abs(hk[0,0,:,:,0]) < 0.3).T
# plt.scatter(x,y)
# plt.show()
# sys.exit()

polynomial = poly.polyfit(fmats[:2], siwdmft[band,2000:2002].imag, deg=1)
print("gamma_dmft: ", -polynomial[0], "   Z_dmft: ", (1-polynomial[1])**(-1))

ReSi0_bz = np.zeros((21,21),dtype = np.float64)
ImSi0_bz = np.zeros((21,21),dtype = np.float64)
gamma_bz = np.zeros((21,21),dtype = np.float64)
Z_bz     = np.zeros((21,21),dtype = np.float64)

# interpolated from 20 to 100 intervals -> 101 values + 1 additional for pcolormesh
ReSi0_plt = np.zeros((102,102),dtype = np.float64)
ImSi0_plt = np.zeros((102,102),dtype = np.float64)
gamma_plt = np.zeros((102,102),dtype = np.float64)
Z_plt     = np.zeros((102,102),dtype = np.float64)

nfreq = 16

for ikx in xrange(20):
    for iky in xrange(20):
        for ikz in xrange(1):
            # this can be made more elegantly but
            # since we do the loop anyways ...
            ReSi0_bz[ikx,iky] = siwdga[band,band,ikx,iky,ikz,frequencies].real
            ImSi0_bz[ikx,iky] = siwdga[band,band,ikx,iky,ikz,frequencies].imag
            polynomial = poly.polyfit(fmats[:2], siwdga[band,band,ikx,iky,ikz,frequencies:frequencies+2].imag, deg=1) # linear interp
            gamma_bz[ikx,iky] = -polynomial[0]
            Z_bz[ikx,iky] = (1-polynomial[1])**(-1)

# Extending the BZ to go from 0 to 2pi
ReSi0_bz[:,20] = ReSi0_bz[:,0]
ReSi0_bz[20,:] = ReSi0_bz[0,:]
ImSi0_bz[:,20] = ImSi0_bz[:,0]
ImSi0_bz[20,:] = ImSi0_bz[0,:]
gamma_bz[:,20] = gamma_bz[:,0]
gamma_bz[20,:] = gamma_bz[0,:]
Z_bz[:,20] = Z_bz[:,0]
Z_bz[20,:] = Z_bz[0,:]

origx   = np.linspace(0,2*np.pi,21)
origy   = np.linspace(0,2*np.pi,21)
interpx = np.linspace(0,2*np.pi,101)
interpy = np.linspace(0,2*np.pi,101)
plotx   = np.linspace(0,2*np.pi,102) # one more to shift the pixels
ploty   = np.linspace(0,2*np.pi,102) # one more to shift the pixels

# Interpolation
fint = interp2d(origx,origy, ReSi0_bz[()], kind='cubic')
data = fint(interpx,interpy)
ReSi0_plt[:51,:51] = np.copy(data[50:,50:])
ReSi0_plt[:51,50:101] = np.copy(data[50:,:51])
ReSi0_plt[50:101,:51] = np.copy(data[:51,50:])
ReSi0_plt[50:101,50:101] = np.copy(data[:51,:51])

fint = interp2d(origx,origy, ImSi0_bz[()], kind='cubic')
data = fint(interpx,interpy)
ImSi0_plt[:51,:51] = np.copy(data[50:,50:])
ImSi0_plt[:51,50:101] = np.copy(data[50:,:51])
ImSi0_plt[50:101,:51] = np.copy(data[:51,50:])
ImSi0_plt[50:101,50:101] = np.copy(data[:51,:51])

fint = interp2d(origx,origy, gamma_bz[()], kind='cubic')
data = fint(interpx,interpy)
gamma_plt[:51,:51] = np.copy(data[50:,50:])
gamma_plt[:51,50:101] = np.copy(data[50:,:51])
gamma_plt[50:101,:51] = np.copy(data[:51,50:])
gamma_plt[50:101,50:101] = np.copy(data[:51,:51])

fint = interp2d(origx,origy, Z_bz[()], kind='cubic')
data = fint(interpx,interpy)
Z_plt[:51,:51] = np.copy(data[50:,50:])
Z_plt[:51,50:101] = np.copy(data[50:,:51])
Z_plt[50:101,:51] = np.copy(data[:51,50:])
Z_plt[50:101,50:101] = np.copy(data[:51,:51])

# A4 - paper size = 8.3 by 11.7 inches
f = plt.figure(figsize=(11,2.5)) # this looks good
# f = plt.figure(figsize=(6.3,2)) # this is the correct size but looks like garbage

ax1 = f.add_subplot(141)
ax1.tick_params(axis='both', which='major', labelsize=14)
ax1.set_title(r'$\mathrm{(a)\enspace Re}[\Sigma(\mathbf{k},i\nu_0)]$')
ax1.set_xlabel(r'$k_x$', labelpad=-1, fontsize=14)
ax1.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
ax1.set_xlim(-np.pi,np.pi)
ax1.set_ylim(-np.pi,np.pi)
ax1.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
p1 = plt.pcolormesh(plotx-np.pi,ploty-np.pi,ReSi0_plt, cmap='rainbow', rasterized=True)
a = ReSi0_plt[:-1,:-1].min()
b = ReSi0_plt[:-1,:-1].max()
cb1 = plt.colorbar(p1, ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.2f', shrink=0.715)
cb1.ax.tick_params(labelsize=14)

ax2 = f.add_subplot(142)
ax2.tick_params(axis='both', which='major', labelsize=14)
ax2.set_title(r'$\mathrm{(b)\enspace Im}[\Sigma(\mathbf{k},i\nu_0)]$')
ax2.set_xlabel(r'$k_x$', labelpad=-1, fontsize=14)
ax2.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
ax2.set_xlim(-np.pi,np.pi)
ax2.set_ylim(-np.pi,np.pi)
ax2.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
p2 = plt.pcolormesh(plotx-np.pi,ploty-np.pi,ImSi0_plt, cmap='rainbow_r', rasterized=True)
a = ImSi0_plt[:-1,:-1].min()
b = ImSi0_plt[:-1,:-1].max()
cb2 = plt.colorbar(p2, ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.715)
cb2.ax.tick_params(labelsize=14)

ax3 = f.add_subplot(143)
ax3.tick_params(axis='both', which='major', labelsize=14)
ax3.set_title(r'$\mathrm{(c)}\enspace\gamma_{\mathbf{k}}$')
ax3.set_xlabel(r'$k_x$', labelpad=-1, fontsize=14)
ax3.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
ax3.set_xlim(-np.pi,np.pi)
ax3.set_ylim(-np.pi,np.pi)
ax3.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
p3 = plt.pcolormesh(plotx-np.pi,ploty-np.pi,gamma_plt, cmap='rainbow', rasterized=True)
a = gamma_plt[:-1,:-1].min()
b = gamma_plt[:-1,:-1].max()
cb3 = plt.colorbar(p3, ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.715)
cb3.ax.tick_params(labelsize=14)

ax4 = f.add_subplot(144)
ax4.tick_params(axis='both', which='major', labelsize=14)
ax4.set_title(r'$\mathrm{(d)}\enspace Z_{\mathbf{k}}$')
ax4.set_xlabel(r'$k_x$', labelpad=-1, fontsize=14)
ax4.set_ylabel(r'$k_y$', labelpad=-5, fontsize=14)
ax4.set_xlim(-np.pi,np.pi)
ax4.set_ylim(-np.pi,np.pi)
ax4.set_aspect('equal')
plt.xticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
plt.yticks([-np.pi,0,np.pi],[r'$-\pi$','$0$', r'$\pi$'])
p4 = plt.pcolormesh(plotx-np.pi,ploty-np.pi,Z_plt, cmap='rainbow_r', rasterized=True)
a = Z_plt[:-1,:-1].min()
b = Z_plt[:-1,:-1].max()
cb4 = plt.colorbar(p4, ticks=[a,a+(b-a)/3,a+(b-a)*2/3, b], format='%.3f', shrink=0.715)
cb4.ax.tick_params(labelsize=14)


plt.tight_layout()
plt.savefig('1x4_gamma_n240.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()
