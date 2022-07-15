import numpy as np
import matplotlib.pyplot as plt
import argparse
import pdb
from scipy import interpolate

mp = 1.6726e-24
alpha = 4.18e-13

def loadout(fname):
    try:
        f = np.loadtxt(fname, skiprows=2)
    except Exception as e:
        rowfail = int(str(e).split('line ')[-1])
        f = np.loadtxt(fname, skiprows=0, max_rows=rowfail-5)
    return f

fig, ax = plt.subplots(1, 1)
r, rho, nC, n_H, ionfrac, column = loadout('cycle400_halfion.txt').T

Gamma = alpha * nC

ax.plot(column, Gamma, 'k--')
ax.set_ylabel(r"ionization rate (s$^{-1}$)")
ax.set_xlabel(r"$N_H$ (cm$^{-2}$)")
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlim((1e16, 1e21))
ax.set_ylim((5e-11, 9e-5))
plt.show()
