import numpy as np
import matplotlib.pyplot as plt
import pdb

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
Pnm = 2.45

def midpoint_diff(t):
  midpoints = 0.5*(t[1:]+t[:-1])
  dt = np.diff(midpoints)
  dt = np.insert(dt, 0, midpoints[1] - midpoints[0])
  dt = np.append(dt, midpoints[-1] - midpoints[-2])
  return dt

def p(t):
    return np.exp(-(tdiff/t)**2) * lo_efreq * Pnm * np.exp(-lo_efreq * t)

def f(t):
    return lo_efreq * Pnm * np.exp(-lo_efreq * t)

def f_int(llim, t):
    return Pnm * (np.exp(-lo_efreq * llim) - np.exp(-lo_efreq * t))

def sample():
    reject = True
    while reject:
        # Sample an area under f(t)
        B_samp = A * np.random.random()

        # Find t for which the area under f(t) to the left of t is equal to B_samp
        t_samp = - 1. / lo_efreq * np.log(np.exp(-lo_efreq * llim) - B_samp/Pnm)

        # Sample a value between 0 and f(t_samp)
        f_samp = f(t_samp) * np.random.random()

        # Accept if f_samp < p(t_samp), reject if f_samp > p(t_samp)
        if f_samp <= p(t_samp):
            reject = False
        else: 
            pass
    return t_samp

def rejection_method(r0, tau0):
    global tlc, tdiff, lo_efreq, A, pnorm, llim, rlim
    tlc = r0 / c
    tdiff = tlc * (a * tau0)**(1./3.)
    lo_efreq = 1./(tcoeff * (a * tau0)**(1./3.))

    tvals = np.logspace(0, 7, 200)

    # Find effective limits of the p(t) function
    test_pvals = p(tvals)
    thresh_inds = np.where(test_pvals >= 1e-6)[0]
    llim = tvals[thresh_inds[0]]
    rlim = tvals[thresh_inds[-1]]

    tvals = np.logspace(np.log10(llim), np.log10(rlim), 200)

    A = f_int(llim, rlim)
    pnorm = np.sum(p(tvals)*midpoint_diff(tvals))

    Nsamps = 100000
    tsamps = np.zeros(Nsamps)
    for i in range(Nsamps):
        tsamps[i] = sample()
    
    pdata = np.array([tvals, p(tvals)/pnorm])
    fdata = np.array([tvals, f(tvals)/pnorm])

    return tvals, pdata, fdata, tsamps

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
    tvals, pdata, fdata, tsamps = rejection_method(r0val, tauval)
    bins = np.logspace(np.log10(llim), np.log10(rlim), 75)
    ax.hist(tsamps, bins=bins, density=True, log=True, histtype='step', alpha=0.75, label='sampled times')
    ax.plot(tvals, p(tvals)/pnorm, alpha=0.75, label='$P(t)$')
    ax.plot(tvals, f(tvals)/pnorm, 'r--', alpha=0.75, label='$f(t)$, comparison function')
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
