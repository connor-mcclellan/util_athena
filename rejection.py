import numpy as np
import matplotlib.pyplot as plt
import pdb

from scipy.integrate import quad

import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 150

c = 2.99792458e10
tgas = 1e4
nu0 = 2.468e15
mass = 1.660538782e-24
kb = 1.380649e-16

vth = np.sqrt( 2 * kb * tgas / mass)
doppwidth = nu0 * vth / c;
lorwidth = 6.265e8/(4.*np.pi)
a = lorwidth / doppwidth

tcoeff = 1.698161839733523

def p(t): # Original function
    return np.exp(-(tdiff/t)**2) * lo_efreq * np.exp(-lo_efreq * t)

def f(t): # Comparison function
    return lo_efreq * np.exp(-lo_efreq * t)

def sample():
    reject = True
    n_rejections = 0
    while reject:
        # Sample an area under f(t)
        B_samp = np.random.random()

        # Find t for which the area under f(t) to the left of t is equal to B_samp
        # Left bound for t integral is light crossing time
        t_samp = - 1. / lo_efreq * np.log(np.exp(-lo_efreq * tlc) - B_samp)

        # Sample a value between 0 and f(t_samp)
        f_samp = f(t_samp) * np.random.random()

        # Accept if f_samp < p(t_samp), reject if f_samp > p(t_samp)
        if f_samp <= p(t_samp):
            reject = False
        else: 
            n_rejections += 1
    return t_samp, n_rejections

def rejection_method(r0, tau0):
    global tlc, tdiff, lo_efreq
    tlc = r0 / c
    tdiff = tlc * (a * tau0)**(1./3.)
    lo_efreq = 1./(tcoeff * (a * tau0)**(1./3.))

    Nsamps = 10000
    tsamps = np.zeros(Nsamps)
    rejections = np.zeros(Nsamps)
    for i in range(Nsamps):
        tsamps[i], rejections[i] = sample()
    
    # Normalize and plot original distribution
    tvals = np.logspace(np.log10(tlc), 3, 200)
    pnorm, _ = quad(p, tlc, np.inf)
    pdata = np.array([tvals, p(tvals)/pnorm])

    return tvals, pdata, tsamps, rejections


r0 = 1e11
tau0 = 1e7

tvals, pdata, tsamps, rejections = rejection_method(r0, tau0)
plt.hist(rejections, bins=25, histtype='step')
plt.xlabel('number of rejections')
plt.ylabel('N')
plt.title('R = 1e11, tau0 = 1e7')
plt.show()

'''
r0=[1e8, 1e9, 1e10, 1e11]
tau0=1e7

if type(r0) == list:
    varytau = False
elif type(tau0) == list:
    varytau = True

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(8, 6))
for i, ax in enumerate([ax1, ax2, ax3, ax4]):
    if varytau:
        tauval = tau0[i]
        r0val = r0
    else:
        tauval = tau0
        r0val = r0[i]
    tvals, pdata, tsamps = rejection_method(r0val, tauval)
    bins = np.logspace(np.log10(tlc), 3, 75)
    ax.hist(tsamps, bins=bins, density=True, log=True, histtype='step', alpha=0.75, label='sampled times')
    ax.plot(*pdata, alpha=0.75, label='$P(t)$')
#    ax.plot(tvals, f(tvals)/pnorm, 'r--', alpha=0.75, label='$f(t)$, comparison function')
    ax.set_xscale('log')
    ax.set_ylim((1e-6, 0.5))
    ax.text(0.7, 0.9, r"$\tau_0 = ${:.0e}".format(tauval), transform=ax.transAxes)
    ax.text(0.7, 0.83, r"$R = ${:.0e}".format(r0val), transform=ax.transAxes)

    ax.set_xlabel('$t$')
    ax.set_ylabel('$P(t)$')

ax4.legend(bbox_to_anchor=(0.0, -0.2), loc='upper center', fontsize='small', ncol=3, frameon=False)
plt.tight_layout()
plt.subplots_adjust(bottom=0.15, wspace=0, hspace=0)
plt.show()
'''
