#!/usr/bin/env python
import numpy as np
import h5py
import matplotlib.pyplot as plt

# make matplot math look like latex
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

f=h5py.File('giwk_dmft.hdf5','r')
w=f['axes/w'][()]
nw=w.shape[0]

kpoints=f['data'].keys()
nk=len(kpoints)

spec_top=np.zeros((nw,nk+1))
spec_bot=np.zeros((nw,nk+1))

for i,kpoint in enumerate(kpoints):
    #spec[:,i]=np.sum(f['data/'+kpoint+'/giw/maxent'][()],axis=0)
    spec_top[:,i]=np.sum(f['data/'+kpoint+'/giw/maxent'][:3],axis=0)
    spec_bot[:,i]=np.sum(f['data/'+kpoint+'/giw/maxent'][3:],axis=0)
    # spec[:,i]=f['data/'+kpoint+'/giw/maxent'][0]
    #spec[:,i]=f['data/'+kpoint+'/giwtr/maxent'][()]

spec_top[:,-1] = spec_top[:,0]
spec_bot[:,-1] = spec_bot[:,0]

f1=plt.figure(figsize=(8,6.2))
ax1=f1.add_subplot(211)
ax1.set_title(r'$\mathrm{top\;layer}$')
ax1.set_ylabel(r'$E\;[\mathrm{eV}]$', labelpad=-1, fontsize=14)
plt.pcolormesh(np.arange(nk+1),w,spec_top,vmin=0, vmax=4, rasterized=True)
plt.xlim(0,120)
plt.xticks([0,40,80,120],[r'$\Gamma$','$X$', '$M$', r'$\Gamma$'], fontsize=14)
plt.ylim(-1,1.5)
cb1 = plt.colorbar(ticks=[0,1,2,3,4], format='%.1f', pad=0.02)
cb1.set_label(r'$A(\mathbf{k},E)\;[\mathrm{arb.\;unit}]$', fontsize=14)

ax2=f1.add_subplot(212)
ax2.set_title(r'$\mathrm{bottom\;layer}$')
ax2.set_ylabel(r'$E\;[\mathrm{eV}]$', labelpad=2, fontsize=14)
plt.pcolormesh(np.arange(nk+1),w,spec_bot,vmin=0, vmax=4, rasterized=True)
plt.xlim(0,120)
plt.xticks([0,40,80,120],[r'$\Gamma$','$X$', '$M$', r'$\Gamma$'], fontsize=14)
plt.ylim(-1,1.5)
cb2 = plt.colorbar(ticks=[0,1,2,3,4], format='%.1f', pad=0.02)
cb2.set_label(r'$A(\mathbf{k},E)\;[\mathrm{arb.\;unit}]$', fontsize=14)

plt.tight_layout()
plt.savefig('dmft_spectrum.eps', format='eps', dpi=600, bbox_inches='tight', Rasterized=True)
plt.show()
