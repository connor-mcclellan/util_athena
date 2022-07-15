from util import Line, read_bin, Params
import numpy as np
from plot_lya import plot_spec_lya
import astropy.constants as c
import astropy.units as u
import argparse
import pdb
import athena_mc_spec as mcspec
import matplotlib.pyplot as plt
from scipy.special import wofz
from scipy.interpolate import interp1d
from scipy.integrate import quad

'''
Must be run from /LyraShared/bcm2vn/athena/outputs/accel with athena mc
plotting modules in the same directory.
'''


def voigtx(a, x, full=False):
    # Full
    if full:
        z = x + a*1j
        H = np.real(wofz(z))
        line_profile = H/np.sqrt(np.pi)
    else:
        line_profile = a / np.pi / (0.01 + x**2)
    return line_profile

def sigma_func(x, a):
    phix = voigtx(a, x, full=True)
    integrand = np.sqrt(2.0/3.0)/phix
    return integrand

def get_sigma(x, a):
    try:
        result = quad(sigma_func,0.0,x, args=(a,))[0]
    except:
        result = []
        for xval in x:
            result.append(quad(sigma_func,0.0,xval, args=(a,))[0])
        result = np.array(result)
    return result

def monte_carlo(x, bins=64):

    # Plot histogram to get n and bin positions, then clear it
    n, bins, patches = plt.hist(x, bins=bins)
    plt.cla()

    # Calculate centers of bins and normalization factor
    bincenters = 0.5 * (bins[1:] + bins[:-1])
    N = np.sum(n)
    norm = 1. / N / (bins[1:] - bins[:-1])
    err = 2. * np.sqrt(n / N**2.) * N

    return bincenters, n * norm, err * norm

def hsp_analytic(x, line, p, temp=1e4, radius=1e11,
                 L=1., tau0=1e4, xi=0.0, normed=True, core=True):

    # Quantities derived from constants
    vth = np.sqrt(2.0 * c.k_B.cgs.value * temp / c.m_p.cgs.value)

    # Quantities derived from arguments
    delta = line.nu0 * vth / c.c.cgs.value
    a = line.gamma / (4.0 * np.pi * delta)
    print('a: ', a)
    if core:
        phix = voigtx(a, x, full=True)    
        sigma = get_sigma(x, a)
        sigma0 = line.strength / (np.sqrt(np.pi) * delta)
        sigmai = get_sigma(xi, a)
    else:
        phix = voigtx(a, x)
        sigma = p.beta * x**3. / a
        sigma0 = line.strength / (np.sqrt(np.pi) * delta)
        sigmai = p.beta * xi**3. / a
    numden = tau0 / (sigma0 * radius)
    kx = numden * line.strength / delta
    z = (sigma - sigmai) / (kx * radius)

    # Equations
    Jprefac = np.sqrt(6.0) / (16.0 * np.pi**3.) * kx**2. * L / delta

    print('Jprefac: ', Jprefac)
    print('sigma: ', sigma)
    print('z: ', z)

    Hsp_analytic = Jprefac / (3.0 * phix) / (kx * radius)**3. * \
        (np.pi**2. / 2.0) / (1.0 + np.cosh(np.pi * z))
    print('cosh: ', 1 + np.cosh(np.pi * z))
    print('Hsp_analytic: ', Hsp_analytic)
    return Hsp_analytic

# Execute main function
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile',
        help='input photon list filename')
    args = parser.parse_args()

    # read spectrum as dict from infile
    spectrum = mcspec.read_spectrum(args.infile)


    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    x, ax = plot_spec_lya(spectrum, 0, ax)
    x_theory = np.linspace(x[0], x[-1], 200)

    # Set up line and get data from bin files
    lya = Line(1215.6701, 0.4164, 6.265e8)
    p = Params(line=lya, temp=1e4, tau0=1e4, num_dens=1701290465.5139434, 
           energy=1., R=1e11, sigma_source=0., n_points=1e4)

    norm = 4.0 * np.pi * p.R**2. * p.delta * 4.0 * np.pi / p.energy
    print('norm: ', norm)

    # H_0 from analytic solution
    sol = hsp_analytic(x_theory, lya, p, core=False) * norm
    sol2 = hsp_analytic(x_theory, lya, p, core=True) * norm
    ax.plot(x_theory, sol, 'g-', linewidth=2, alpha=0.5, label=r'$H_0$')
    ax.plot(x_theory, sol2, 'm-', linewidth=2, alpha=0.5, label=r'$H_0$, full voigt')
    plt.legend()
    plt.show()
    print("sol: ", sol)
