#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import numpy.polynomial.polynomial as poly
from scipy.interpolate import interp2d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import h5py
import sys

# make matplot math look like latex
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

asy = h5py.File('adga-20180126-013536.333-output.hdf5','r')

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

# A4 - paper size = 8.3 by 11.7 inches
f = plt.figure(figsize=(7,8)) # this looks good

ax1 = f.add_subplot(421)
ax1.set_xlim(0,120)
ax1.set_ylim(0.12,0.20)
ax1.set_ylabel(r'$\chi_{D\;1111}$', labelpad=15, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
plt.yticks([0.12,0.14,0.16,0.18,0.20])

plt.plot(np.arange(41), susc_dens[0,0,0,np.arange(41),0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(40,81), susc_dens[0,0,np.arange(0,41),40,0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(80,121), susc_dens[0,0,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'blue', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')
ax1.get_yaxis().set_label_coords(-0.2,0.5)

ax1 = f.add_subplot(422)
ax1.set_xlim(0,120)
# ax1.set_ylim(-400,400)
ax1.set_ylabel(r'$\chi_{M\;1111}$', labelpad=8, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
# plt.yticks([-400,-200,0,200,400], fontsize=14)
ax1.get_yaxis().set_label_coords(-0.2,0.5)

plt.plot(np.arange(41), susc_magn[0,0,0,np.arange(41),0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(40,81), susc_magn[0,0,np.arange(0,41),40,0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(80,121), susc_magn[0,0,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'red', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')

ax1 = f.add_subplot(423)
ax1.set_xlim(0,120)
# ax1.set_ylim(0.12,0.20)
ax1.set_ylabel(r'$\chi_{D\;2222}$', labelpad=8, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
# plt.yticks([0.12,0.14,0.16,0.18,0.20], fontsize=14)

plt.plot(np.arange(41), susc_dens[1,1,0,np.arange(41),0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(40,81), susc_dens[1,1,np.arange(0,41),40,0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(80,121), susc_dens[1,1,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'blue', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')
ax1.get_yaxis().set_label_coords(-0.2,0.5)

axins = zoomed_inset_axes(ax1, 6, loc='lower right')  # zoom = 3
axins.plot(np.arange(68,81), susc_dens[1,1,np.arange(28,41),40,0,frequencies].real, 'black', lw=2,)
# sub region of the original image
x1, x2, y1, y2 = 68, 74, 0.0275, 0.031
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# fix the number of ticks on the inset axes
plt.xticks([])
plt.yticks([])
mark_inset(ax1, axins, loc1=3, loc2=1, fc="none", ec="0.5")

ax1 = f.add_subplot(424)
ax1.set_xlim(0,120)
# ax1.set_ylim(-400,400)
ax1.set_ylabel(r'$\chi_{M\;2222}$', labelpad=8, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
# plt.yticks([-400,-200,0,200,400], fontsize=14)
ax1.get_yaxis().set_label_coords(-0.2,0.5)

plt.plot(np.arange(41), susc_magn[1,1,0,np.arange(41),0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(40,81), susc_magn[1,1,np.arange(0,41),40,0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(80,121), susc_magn[1,1,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'red', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')

ax1 = f.add_subplot(425)
ax1.set_xlim(0,120)
# ax1.set_ylim(0.12,0.20)
ax1.set_ylabel(r'$\chi_{D\;4444}$', labelpad=8, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
# plt.yticks([0.12,0.14,0.16,0.18,0.20], fontsize=14)

plt.plot(np.arange(41), susc_dens[3,3,0,np.arange(41),0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(40,81), susc_dens[3,3,np.arange(0,41),40,0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(80,121), susc_dens[3,3,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'blue', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')
ax1.get_yaxis().set_label_coords(-0.2,0.5)

ax1 = f.add_subplot(426)
ax1.set_xlim(0,120)
# ax1.set_ylim(-400,400)
ax1.set_ylabel(r'$\chi_{M\;4444}$', labelpad=8, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
# plt.yticks([-400,-200,0,200,400], fontsize=14)

plt.plot(np.arange(41), susc_magn[3,3,0,np.arange(41),0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(40,81), susc_magn[3,3,np.arange(0,41),40,0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(80,121), susc_magn[3,3,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'red', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')
ax1.get_yaxis().set_label_coords(-0.2,0.5)

ax1 = f.add_subplot(427)
ax1.set_xlim(0,120)
# ax1.set_ylim(0.12,0.20)
ax1.set_ylabel(r'$\chi_{D\;5555}$', labelpad=8, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
# plt.yticks([0.12,0.14,0.16,0.18,0.20], fontsize=14)
ax1.get_yaxis().set_label_coords(-0.2,0.5)

plt.plot(np.arange(41), susc_dens[4,4,0,np.arange(41),0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(40,81), susc_dens[4,4,np.arange(0,41),40,0,frequencies].real, 'blue', lw='2')
plt.plot(np.arange(80,121), susc_dens[4,4,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'blue', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')

ax1 = f.add_subplot(428)
ax1.set_xlim(0,120)
# ax1.set_ylim(-400,400)
ax1.set_ylabel(r'$\chi_{M\;5555}$', labelpad=8, fontsize=14)
plt.xticks([0,40,80,120], [r'$\Gamma$', r'$X$', r'$M$', r'$\Gamma$'], fontsize=14)
# plt.yticks([-400,-200,0,200,400], fontsize=14)

plt.plot(np.arange(41), susc_magn[4,4,0,np.arange(41),0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(40,81), susc_magn[4,4,np.arange(0,41),40,0,frequencies].real, 'red', lw='2')
plt.plot(np.arange(80,121), susc_magn[4,4,np.arange(40,-1,-1),np.arange(40,-1,-1),0,frequencies].real , 'red', lw='2')
plt.axhline(y=0, xmin=0, xmax=120, color='black')
ax1.get_yaxis().set_label_coords(-0.2,0.5)

axins = zoomed_inset_axes(ax1, 6, loc='upper left')  # zoom = 3
axins.plot(np.arange(58,68), susc_magn[4,4,np.arange(18,28),40,0,frequencies].real, 'black', lw=2,)
# sub region of the original image
x1, x2, y1, y2 = 58, 65, 2.6, 3.1
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# fix the number of ticks on the inset axes
plt.xticks([])
plt.yticks([])
mark_inset(ax1, axins, loc1=3, loc2=1, fc="none", ec="0.5")

plt.tight_layout()
plt.savefig('layer-diverging-susc.eps', bbox_inches='tight', format='eps', dpi=600, Rasterized=True)
plt.show()
