import numpy as np
import matplotlib.pyplot as plt
import argparse


def loadout(fname):
    try:
        f = np.loadtxt(fname, skiprows=4)
    except Exception as e:
        rowfail = int(str(e).split('line ')[-1])
        f = np.loadtxt(fname, skiprows=4, max_rows=rowfail-5)
    return f

f = loadout('out_taunu.txt')
logbins = 10**np.linspace(np.log10(min(f)), np.log10(max(f))) 
n1, bins1, patches = plt.hist(f, bins=logbins, histtype='step', label=r'$\tau(\nu)$')

f = loadout('out_taunu0.txt')
logbins = 10**np.linspace(np.log10(min(f)), np.log10(max(f))) 
n2, bins2, patches = plt.hist(f, bins=logbins, histtype='step', label=r'$\tau(\nu_0)$')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()


