#! /usr/bin/env python

from __future__ import print_function, division, absolute_import
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import h5py
import sys

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

dmft_cond = h5py.File('svo2-kanamori-run00-left2right-FULL-2017-12-18-Mon-14-33-56.hdf5','r')
dmft_ins = h5py.File('svo2-run03-wormmeas-full-2017-11-06-Mon-13-55-47.hdf5','r')


dmft_0 = []
dmft_0.append(dmft_cond['dmft-001/ineq-001/gtau/value'])
dmft_0.append(dmft_cond['dmft-001/ineq-002/gtau/value'])

dmft_1 = []
dmft_1.append(dmft_ins['stat-001/ineq-001/gtau/value'])
dmft_1.append(dmft_ins['stat-001/ineq-002/gtau/value'])


f = plt.figure(figsize=(8,7))
beta = 38

ax1 = f.add_subplot(211)
ax1.set_title(r"$U' = 2.8 eV \;(\mathrm{conducting})$", fontsize=14)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_ylabel(r'$G(\tau)$', fontsize=14)
ax1.set_xlabel(r'$\tau$', fontsize=14)
ax1.set_ylim(0,1)
ax1.set_xlim(0,beta)
ax1.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.xticks([0,beta/4,beta/2, beta*3/4, beta])
plt.plot(np.linspace(0,beta,1000), dmft_0[0][0,0,:], 'r', lw=2, label='top layer')
plt.plot(np.linspace(0,beta,1000), dmft_0[0][1,0,:], 'r-.', lw=2)
plt.plot(np.linspace(0,beta,1000), dmft_0[1][0,0,:], 'b', lw=2, label='bottom layer')
plt.plot(np.linspace(0,beta,1000), dmft_0[1][1,0,:], 'b-.', lw=2)
plt.legend(loc = 'upper right', frameon=False, fontsize=16)


axins = zoomed_inset_axes(ax1, 15, loc='upper center')  # zoom = 15
axins.plot(np.linspace(0,beta,1000), dmft_0[0][0,0,:], 'r', lw=2, label='top layer')
axins.plot(np.linspace(0,beta,1000), dmft_0[0][1,0,:], 'r-.', lw=2)
axins.plot(np.linspace(0,beta,1000), dmft_0[1][0,0,:], 'b', lw=2, label='bottom layer')
axins.plot(np.linspace(0,beta,1000), dmft_0[1][1,0,:], 'b-.', lw=2)
plt.axvline(x=19.0, linestyle='--', color='gray')

# sub region of the original image
x1, x2, y1, y2 = 18.5, 19.5, 0.02, 0.05
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# fix the number of ticks on the inset axes
axins.yaxis.get_major_locator().set_params(nbins=3)
axins.xaxis.get_major_locator().set_params(nbins=3)

# plt.xticks(visible=False)
# plt.yticks(visible=False)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax1, axins, loc1=3, loc2=4, fc="none", ec="0.5")



ax2 = f.add_subplot(212)
ax2.set_title(r"$U' = 3.5 eV \;(\mathrm{insulating})$", fontsize=14)
ax2.tick_params(axis='both', which='major', labelsize=16)
ax2.set_ylabel(r'$G(\tau)$', fontsize=14)
ax2.set_xlabel(r'$\tau$', fontsize=14)
ax2.set_ylim(0,1)
ax2.set_xlim(0,beta)
ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
plt.xticks([0,beta/4,beta/2, beta*3/4, beta])
plt.plot(np.linspace(0,beta,2000), dmft_1[0][0,0,:], 'r', lw=2, label='top layer')
plt.plot(np.linspace(0,beta,2000), dmft_1[0][1,0,:], 'r-.', lw=2)
plt.plot(np.linspace(0,beta,2000), dmft_1[1][0,0,:], 'b', lw=2, label='bottom layer')
plt.plot(np.linspace(0,beta,2000), dmft_1[1][1,0,:], 'b-.', lw=2)

plt.legend(loc = 'upper right', frameon=False, fontsize=16)

axins = zoomed_inset_axes(ax2, 15, loc='upper center')  # zoom = 15
axins.plot(np.linspace(0,beta,2000), dmft_1[0][0,0,:], 'r', lw=2, label='top layer')
axins.plot(np.linspace(0,beta,2000), dmft_1[0][1,0,:], 'r-.', lw=2)
axins.plot(np.linspace(0,beta,2000), dmft_1[1][0,0,:], 'b', lw=2, label='bottom layer')
axins.plot(np.linspace(0,beta,2000), dmft_1[1][1,0,:], 'b-.', lw=2)
plt.axvline(x=19.0, linestyle='--', color='gray')

# sub region of the original image
x1, x2, y1, y2 = 18.5, 19.5, 0.00, 0.03
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
# fix the number of ticks on the inset axes
axins.yaxis.get_major_locator().set_params(nbins=3)
axins.xaxis.get_major_locator().set_params(nbins=3)

# plt.xticks(visible=False)
# plt.yticks(visible=False)

# draw a bbox of the region of the inset axes in the parent axes and
# connecting lines between the bbox and the inset axes area
mark_inset(ax2, axins, loc1=3, loc2=4, fc="none", ec="0.5")


plt.tight_layout()
plt.savefig('dmft_gtau.eps', format='eps', dpi=600, rasterized=True)
plt.show()
