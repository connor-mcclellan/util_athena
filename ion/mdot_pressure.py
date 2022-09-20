import numpy as np
import matplotlib.pyplot as plt
import pdb

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150

# Constants
mp = 1.6726e-24
kb = 1.380649e-16

def mdot_pres(GM, temp, rbase, pbase, mmw=1):
    rhobase = pbase * mmw * mp / kb / temp
    a = np.sqrt(kb * temp / mmw / mp)
    rs = GM / 2 / a**2
    mdot = 4 * np.pi * rs**2 * a * rhobase * np.exp(3./2 - 2*rs/rbase)
    print("SONIC POINT: ", rs)
    pdb.set_trace()
    return mdot


GM = 1.2668653e23
rbase = 1e10
pbase = 10
temp = 1e4


fig, axs = plt.subplots(2, 2, figsize=(8, 6))
axs = np.ndarray.flatten(axs)
colors = ['r', 'orange', 'g', 'c']
xvars = [GM, temp, rbase, pbase]
varnames = ['GM', 'T', r'$r_{base}$', r'$P_{base}$']

for i, xvar in enumerate(xvars):
    x_orig = xvar
    y_orig = mdot_pres(GM, temp, rbase, pbase)

    xvar *= np.linspace(0.5, 2, 100)
    xvars[i] = xvar
    yvar = mdot_pres(*xvars)

    xvars = [GM, temp, rbase, pbase]
    axs[i].plot(xvar, yvar, c=colors[i])
    axs[i].set_xlabel(varnames[i])
    axs[i].set_ylabel('$\dot M$')

    axs[i].axvline(x_orig, c=colors[i], alpha=0.5)

    axs[i].set_yscale('log')

plt.tight_layout()
plt.show()
