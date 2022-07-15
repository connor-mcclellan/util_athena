import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
from scipy import interpolate

mp = 1.6726e-24

def loadout(fname):
    try:
        f = np.loadtxt(fname, skiprows=2)
    except Exception as e:
        rowfail = int(str(e).split('line ')[-1])
        f = np.loadtxt(fname, skiprows=0, max_rows=rowfail-5)
    return f

fig, ax = plt.subplots(1, 1)
r, rho, nC, n_H, ionfrac, column = loadout('cycle400_halfion.txt').T

ion_cross_idx = np.argwhere(np.diff(np.sign(ionfrac - 0.5))).flatten()
nCH_cross_idx = np.argwhere(np.diff(np.sign(n_H - nC))).flatten()

ax.plot(r, n_H, 'b--', label=r'$n_H$')
ax.plot(r, nC, 'r--', label=r'$n_c$')
#ax.axvline(r[nCH_cross_idx], c='m', alpha=0.5)
#ax.axhline(nC[nCH_cross_idx], c='m', alpha=0.5)
ax.set_ylabel("number density", c='b')
ax.set_xlabel("r")
ax.set_yscale('log')
ax.set_xscale('log')
ax.legend(loc='lower left')

ax2 = ax.twinx()
ax2.plot(r, ionfrac, 'g-', label='ionization fraction')
#ax2.axvline(r[ion_cross_idx], c='limegreen', alpha=0.5)
ax2.axhline(0.5, c='limegreen', alpha=0.5)
ax2.set_ylabel("ionization fraction", c='g')
ax2.set_ylim((0, 1))
ax2.set_xlim((1e10, 1e11))
ax2.set_xscale('log')
ax2.set_yscale('linear')
ax2.legend(loc='lower right')
plt.title('last timestep')
plt.show()
