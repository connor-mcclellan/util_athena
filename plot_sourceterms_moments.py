import numpy as np
from athena_read import athdf
import pdb
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"text.usetex": True})
import argparse
from pathlib import Path
from glob import glob

def extract_input_val(string):
    return string.split("=")[1].split("#")[0].strip()

def plot_frame(prim, mcmom, mcsrc, logpath, basename, cycle):
    radius = prim['x1v']

    # Flux radiative acceleration data
    radaccelf = mcsrc['RadforceF1'] / mcmom['rho'] #prim['r0'] * prim['rho'] / mcmom['rho']
    radaccelf = np.median(radaccelf, axis=(0))
    radaccelf = radaccelf[1:]
    radaccelf_avg = np.median(radaccelf, axis=(0))

    # Scattering radiative acceleration data
    radaccels = mcsrc['RadforceS1'] / mcmom['rho']
    radaccels = np.median(radaccels, axis=(0))
    radaccels = radaccels[1:]
    radaccels_avg = np.median(radaccels, axis=(0))

    # Density
    rho = prim['rho']
    rho = np.median(rho, axis=0)

    # Ionization fraction
    ionfrac = prim['r0'] / mcmom['rho'] * prim['rho']
    ionfrac = np.median(ionfrac, axis=0)

    # Radial velocity
    vel1 = prim['vel1']
    vel1 = np.median(vel1, axis=0)

    # Flux and energy density
    # do volume average - x1f gives you faces, construct exact volume of every cell and then 
    # weight every cell by its local volume and then average
    er = np.median(mcmom['Ermc'],axis=(0))*2.99792458e10 
    fr = np.median(mcmom['Frmc1'],axis=(0))

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=True, figsize=(12, 8))
    array_dims = np.shape(radaccels)
    for j in range(array_dims[0]):
        ax1.scatter(radius, radaccelf[j, :], color='c', marker='o', s=3, alpha=0.25)
        ax1.scatter(radius, radaccels[j, :], color='orange', marker='o', s=3, alpha=0.25)
        ax2.scatter(radius, er[j, :], color='k', marker='o', s=3, alpha=0.25)
        ax2.scatter(radius, fr[j, :], color='r', marker='o', s=3, alpha=0.25)
        ax3.scatter(radius, rho[j, :], color='g', marker='o', s=3, alpha=0.25)
        ax4.scatter(radius, ionfrac[j, :], color='m', marker='o', s=3, alpha=0.25)
    ax1.plot(radius, np.abs(radaccelf_avg), color='c', marker='o', ms=3)
    ax1.plot(radius, np.abs(radaccels_avg), color='orange', marker='o', ms=3)
    ax1.plot(radius, np.sqrt(np.median(radaccelf**2, axis=(0))), color='c', marker='o', ms=3, ls='--')
    ax1.plot(radius, np.sqrt(np.median(radaccels**2, axis=(0))), color='orange', marker='o', ms=3, ls='--')

    ax3.set_xlabel('r')
    ax4.set_xlabel('r')

    ax1.set_ylabel('acceleration (cm s-2)')
    ax2.set_ylabel(r'$F_r$, $cE_r$')
    ax3.set_ylabel('density (g/cm3)')
    ax4.set_ylabel('ionization fraction')

    ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax3.set_xscale('log')
    ax4.set_xscale('log')

    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')

    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position("right")

    ax1.set_ylim((1e-8, 1e0))
    ax2.set_ylim((1e-2, 1e4))
    ax4.set_ylim((0., 1.1))

    nphot, nx1, nx2, nx3, pbase, temp, accel, runtime = [np.nan]*8
    with open(logpath, 'r') as f:
        log = f.readlines()
        for l in log:
            if "nphot" in l:
                nphot = int(extract_input_val(l))
            if "nx1" in l and "Number of zones" in l:
                nx1 = int(extract_input_val(l))
            if "nx2" in l and "Number of zones" in l:
                nx2 = int(extract_input_val(l))
            if "nx3" in l and "Number of zones" in l:
                nx3 = int(extract_input_val(l))
            if "pbase" in l:
                pbase = float(extract_input_val(l))
            if "temp" in l:
                temp = float(extract_input_val(l))
            if "acceleration" in l:
                accel = str(extract_input_val(l) == 'true')
            if "cpu time used" in l:
                runtime = float(extract_input_val(l))

    captiontext = r"\noindent \textbf{{Number of photons:}} {:.0e}. \textbf{{Mesh:}} {} $\times$ {} $\times$ {}.\\ \textbf{{Base pressure:}} {:.2e} dyne/cm$^2$. \textbf{{Temperature:}} {:.0e} K.\\ \textbf{{Code acceleration:}} {}. \textbf{{Wall time:}} {:.2f} hrs.".format(nphot, nx1, nx2, nx3, pbase, temp, accel, runtime/60./60.)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0, bottom=0.15, wspace=0.)
    plt.suptitle(args.title)

    formatlegend = [Line2D([1], [0], color='c', marker='o', linestyle="None", label=r'Path-length radiative acceleration'), 
                    Line2D([1], [0], color='orange', marker='o', linestyle="None", label=r'Scattering radiative acceleration'), 
                    Line2D([1], [0], color='k', marker='o', linestyle='-', label='Mean value'), 
                    Line2D([1], [0], color='k', marker='o', linestyle='--', label='RMS value')]
    fmtlegend = ax1.legend(handles=formatlegend, loc='upper right', bbox_to_anchor=[1.0, 0.98], frameon=False)
    formatlegend = [Line2D([1], [0], color='k', marker='o', linestyle="None", label=r'$c E_r$'), 
                    Line2D([1], [0], color='r', marker='o', linestyle="None", label=r'$F_r$')]
    fmtlegend = ax2.legend(handles=formatlegend, loc='upper right', bbox_to_anchor=[1.0, 0.98], frameon=False)
    formatlegend = [Line2D([1], [0], color='g', marker='o', linestyle="None", label=r'Density')]
    fmtlegend = ax3.legend(handles=formatlegend, loc='upper right', bbox_to_anchor=[1.0, 0.98], frameon=False)
    formatlegend = [Line2D([1], [0], color='m', marker='o', linestyle="None", label=r'Ionization fraction')]
    fmtlegend = ax4.legend(handles=formatlegend, loc='upper right', bbox_to_anchor=[1.0, 0.98], frameon=False)
    ax4.text(-1, -1.15, captiontext, transform=ax2.transAxes, fontsize=12)
    plt.savefig(basename+'den-vel.{:05d}.png'.format(cycle))
    plt.close()
#    plt.show()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("basename", type=str)
    parser.add_argument("-t", "--title", type=str)
    parser.add_argument("-d", "--dynamic", action='store_true', default=False)
    args = parser.parse_args()

    if args.dynamic:
        nfiles = len(glob(args.basename+'.out1.*.athdf'))
    else:
        nfiles = 1

    for cycle in range(nfiles):
        cyclestr = "{:05d}.athdf".format(cycle)
        prim = athdf(args.basename+'.out1.'+cyclestr, quantities=["rho", 'r0', 'vel1', 'press'])
        mcmom = athdf(args.basename+'.out2.'+cyclestr, quantities=["Ermc", "Frmc1", "rho", "tgas"])
        mcsrc = athdf(args.basename+'.out3.'+cyclestr, quantities=["RadforceF1", "RadforceS1"])

        # Load run info from log
        logpath = Path(args.basename+'.txt')
        if not logpath.exists():
            logpath = Path("output.txt")

        plot_frame(prim, mcmom, mcsrc, logpath, args.basename, cycle)
