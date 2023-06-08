from athena_read import athdf
from glob import glob
import argparse
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
#from matplotlib import rc
import numpy as np
import astropy.constants as c
import astropy.units as u
import multiprocessing as mp
from pathlib import Path
import pdb

#rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
#rc("text", usetex=True)

# ========================================================================================
#                                                                                        -
# Use the following command to produce a .mp4 from this file's outputs:                  -
#                                                                                        -
# > cat $(find . -name "[name]*.png" -print | sort) | ffmpeg -framerate 60 -i - out.mp4  -
#                                                                                        -
# Replace [name] with your files' names.                                                 -
#                                                                                        -
# ========================================================================================

temp = 1e4

# fundamental constants
cgrav=6.67429e-8
hplanck=6.62606896e-27
clight=2.997924589e10
kboltz=1.3806504e-16
charge=4.80320427e-10
abohr=0.52917720859e-8
melectron=9.10938215e-28
mproton=1.672621637e-24
amu=1.660538782e-24

# Lyman alpha stuff
lambda0 = 1215.6701 * 1.e-8
osc_strength = 0.4164
Gamma= 6.265e8
nu0 = clight / lambda0  
line_strength = np.pi*charge**2/(melectron*clight) * osc_strength
vth = np.sqrt( 2.0 * kboltz * temp / amu )
delta = nu0*vth/clight



def midpoint_diff(t):
    midpoints = 0.5 * (t[1:] + t[:-1])
    dt = np.diff(midpoints)
    dt = np.insert(dt, 0, midpoints[1] - midpoints[0])
    dt = np.append(dt, midpoints[-1] - midpoints[-2])
    return dt


def find_filenames(path):
    files = sorted(glob(path+'/*'))
    files = [str(Path(f).resolve()) for f in files if "xdmf" not in f and "txt" not in f]
    basename = files[0].split(".")[0]

    uovfiles = sorted(glob(basename + ".uov0.*.athdf"))
    primfiles = sorted(glob(basename + ".out1.*.athdf"))
    momfiles = sorted(glob(basename + ".out2.*.athdf"))
    mcsrcfiles = sorted(glob(basename + ".out3.*.athdf"))
    return uovfiles, primfiles, momfiles, mcsrcfiles


def get_times(prim, istart, iend):
    times = []
    for i in range(istart, iend):
        # Load in data from each output block
        primdata = athdf(prim[i])
        times.append(primdata["Time"])
    return np.array(times)


def plot_frame(istart, iend, times1, times2):
    for i in range(istart, iend):

        uovdata = athdf(uov1[i])
        primdata = athdf(prim1[i])
        momdata = athdf(mom1[i])
        srcdata = athdf(src1[i])

        this_time = primdata["Time"]
        time_match_ind = np.argwhere(times2 > this_time)[0][0]
        
        uov2before = athdf(uov2[time_match_ind-1])
        prim2before = athdf(prim2[time_match_ind-1])
        mom2before = athdf(mom2[time_match_ind-1])
        src2before = athdf(src2[time_match_ind-1])

        uov2after = athdf(uov2[time_match_ind])
        prim2after = athdf(prim2[time_match_ind])
        mom2after = athdf(mom2[time_match_ind])
        src2after = athdf(src2[time_match_ind])

        # Apply a linear interpolation to the output2 data to match the time of output1

        # Get relevant quantities from data
        radius = primdata["x1v"]
        column = uovdata['column']
        radial_column = np.mean(np.mean(column, axis=0), axis=0)
        sigma = line_strength / np.sqrt(np.pi) / delta
        tau = sigma * radial_column
        # temp = primdata['temp']

        # reconstruct volume from cell face coordinates
        r = primdata["x1f"][np.newaxis, np.newaxis, :]
        th = primdata["x2f"][np.newaxis, :, np.newaxis]
        ph = primdata["x3f"][:, np.newaxis, np.newaxis]
        vol = (
            1 / 3
            * (r[:, :, 1:] ** 3 - r[:, :, :-1] ** 3)
            * (np.cos(th[:, :-1, :]) - np.cos(th[:, 1:, :]))
            * (ph[1:, :, :] - ph[:-1, :, :])
        )
        vol = vol.astype(np.float)
        volshells = np.sum(np.sum(vol, axis=0), axis=0)
        totvol = np.sum(volshells)

        # Energy density and flux
        er = np.sum(np.sum(momdata["Ermc"]*vol, axis=0), axis=0) / volshells
        fr = np.sum(np.sum(momdata["Frmc1"]*vol, axis=0), axis=0) / volshells

        # Flux radiative acceleration data
        radaccelf = srcdata["RadforceF1"] / momdata["rho"] * vol
        radaccelf = np.sum(radaccelf, axis=(0))
        radaccelf_avg = np.sum(radaccelf, axis=(0)) / volshells

        # Scattering radiative acceleration dataz
        radaccels = srcdata["RadforceS1"] / momdata["rho"] * vol
        radaccels = np.sum(radaccels, axis=(0))
        radaccels_avg = np.sum(radaccels, axis=(0)) / volshells

        # Density
        rho = primdata["rho"]
        rho = np.mean(np.mean(rho, axis=0), axis=0)

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 6), sharex=True, dpi=100)

        # PLOT 1: Acceleration & Hydro quantities
        # =======================================

        # Panel 1.5: Acceleration for first output
        # -----------------------
        ax1.plot(tau, radaccelf_avg, **kwa, c="c", label="rad accel flux")
        ax1.plot(tau, -radaccelf_avg, **kwaneg, c="c")
        ax1.plot(tau, radaccels_avg, **kwa, c="orange", label="rad accel scattering")
        ax1.plot(tau, -radaccels_avg, **kwaneg, c="orange")

        # Panel formatting
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylim((1e-9, 1e5))
        #ax1.set_ylim((1e-9, 1e9))
        #ax1.set_xlim((1.1999, 1.2015))
        ax1.grid(which='both', alpha=0.75, lw=0.3)
        ax1.set_ylabel("a [dyne cm${^-3}$]")
        ax1.legend(**lega)
        ax1.set_title((prim1[i]).split("/")[-1])
        ax1.set_xlabel("optical depth")
        ax1.invert_xaxis()

        # Panel 1.2: Acceleration for second output
        # -----------------------

        # Flux radiative acceleration data
        radaccelf = src2before["RadforceF1"] / mom2before["rho"] * vol
        radaccelf = np.sum(radaccelf, axis=(0))
        radaccelf_avg = np.sum(radaccelf, axis=(0)) / volshells

        # Scattering radiative acceleration dataz
        radaccels = src2before["RadforceS1"] / mom2before["rho"] * vol
        radaccels = np.sum(radaccels, axis=(0))
        radaccels_avg = np.sum(radaccels, axis=(0)) / volshells

        # Radiation acceleration
        ax2.plot(tau, radaccelf_avg, **kwa, c="c", label="rad accel flux")
        ax2.plot(tau, -radaccelf_avg, **kwaneg, c="c")
        ax2.plot(tau, radaccels_avg, **kwa, c="orange", label="rad accel scattering")
        ax2.plot(tau, -radaccels_avg, **kwaneg, c="orange")

        # Panel formatting
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_ylim((1e-9, 1e5))
        ax2.grid(which='both', alpha=0.75, lw=0.3)
        ax2.yaxis.tick_right()
        ax2.set_ylabel("a [dyne cm${^-3}$]")
        ax2.legend(**lega)
        ax2.set_title((prim2[time_match_ind-1]).split("/")[-1])
        ax2.set_xlabel("optical depth")
        ax2.invert_xaxis()

        # Panel 1.3: Radiation moments for first output
        # ----------------------------
        # Energy density and flux
        er = np.sum(np.sum(momdata["Ermc"]*vol, axis=0), axis=0) / volshells
        fr = np.sum(np.sum(momdata["Frmc1"]*vol, axis=0), axis=0) / volshells

        ax3.plot(tau, c.c.cgs * er, **kwa, c='k', label=r'$cE_r$')
        ax3.plot(tau, -c.c.cgs * er, **kwaneg, c='k')
        ax3.plot(tau, fr, **kwa, c='r', label='$F_r$')
        ax3.plot(tau, -fr, **kwaneg, c='r')
        ax3.xaxis.set_major_formatter(ScalarFormatter())

        # Panel formatting
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_ylim((1e-4, 1e4))
        ax3.grid(which='both', alpha=0.75, lw=0.3)
        ax3.set_ylabel("[erg cm$^{-2}$ s$^{-1}$]")
        ax3.yaxis.set_label_position("right")
        ax3.legend(**lega)
        ax3.set_xlabel("optical depth")
        ax3.tick_params(axis='x', which='minor', labelsize=7)

        captiontext = "t={:.4e}".format(primdata["Time"])
        ax3.text(0, -.25, captiontext, transform=ax3.transAxes, fontsize=12)

        # Panel 1.4: Radiation moments for second output
        # ----------------------------

        # Energy density and flux
        er = np.sum(np.sum(mom2before["Ermc"]*vol, axis=0), axis=0) / volshells
        fr = np.sum(np.sum(mom2before["Frmc1"]*vol, axis=0), axis=0) / volshells

        ax4.plot(tau, c.c.cgs * er, **kwa, c='k', label=r'$cE_r$')
        ax4.plot(tau, -c.c.cgs * er, **kwaneg, c='k')
        ax4.plot(tau, fr, **kwa, c='r', label='$F_r$')
        ax4.plot(tau, -fr, **kwaneg, c='r')

        # Panel formatting
        ax4.set_xscale("log")
        ax4.set_yscale("log")
        ax4.set_ylim((1e-4, 1e4))
        #ax4.set_ylim((1e23, 1e31))
        ax4.grid(which='both', alpha=0.75, lw=0.3)
        ax4.set_ylabel("[erg cm$^{-2}$ s$^{-1}$]")
        ax4.yaxis.set_label_position("right")
        ax4.yaxis.tick_right()
        ax4.legend(**lega)
        ax4.set_xlabel("optical depth")
        ax4.xaxis.set_major_formatter(ScalarFormatter())
        ax4.xaxis.set_minor_formatter(ScalarFormatter())
        ax4.tick_params(axis='x', which='minor', labelsize=7)
        
        captiontext = "t={:.4e}".format(prim2before["Time"])
        ax4.text(0, -0.25, captiontext, transform=ax4.transAxes, fontsize=12)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0, bottom=0.15, wspace=0.0)
        filename = uov1[i]
        filename2 = (uov2[time_match_ind-1]).split("/")[-1].replace(".uov0", "").replace(".athdf", "")
        framename = filename.replace(".athdf", "_"+filename2+".png").replace(".uov0", "")

        if not args.show:
            plt.savefig(framename, dpi=150)
            plt.close()
        else:
            plt.show()


if __name__ == "__main__":

#    if args.cores > 1:
#        batch = int(len(uovfiles)/args.cores)
#        pool = mp.Pool(processes=args.cores)
#        for j in range(args.cores):
#            istart = int(j*batch)
#            iend = int((j+1)*batch)
#            if j == args.cores-1:
#                iend = len(uovfiles)
#            print(istart, iend)
#            pool.apply_async(plot_frame, args=(istart, iend))
#        pool.close()
#        pool.join()
#    else:

    parser = argparse.ArgumentParser()
    parser.add_argument("path1", type=str)
    parser.add_argument("path2", type=str)
    parser.add_argument("--show", "-w", action='store_true', default=False)
    parser.add_argument("--cores", "-c", type=int, default=1)
    args = parser.parse_args()

    uov1, prim1, mom1, src1 = find_filenames(args.path1)
    uov2, prim2, mom2, src2 = find_filenames(args.path2)

    ns_mass = 1.4 * u.M_sun.cgs

    # Characteristic quantities
    temp0 = (c.c.cgs**2 * c.m_p.cgs / c.k_B.cgs).to(u.K)
    crat = 1.0
    arad = 4 * c.sigma_sb.cgs / c.c.cgs
    kes = 0.2 * u.cm**2 * u.g**(-1)
    rg = (2 * c.G.cgs * ns_mass / (c.c.cgs**2)).to(u.cm)
    prat = (arad * temp0**4 * kes * rg / c.c.cgs**2).decompose()
    rho0 = (1 / kes / rg).to(u.g / u.cm**3)

    # Default plot kwargs and legend kwargs
    kwa = {'alpha': 0.75, 'lw': 1}
    kwaneg = {'alpha': 0.5, 'lw': 1, 'ls': '--'}
    #lega = {'bbox_to_anchor': (1.02, 0.8), 'loc': "upper left", 'fontsize': "x-small", 'frameon': False}
    lega = {'frameon': False, "fontsize": "x-small"}

    times1 = get_times(prim1, 0, len(uov1))
    times2 = get_times(prim2, 0, len(uov2))

    if len(times1) > len(times2):
        times1, times2 = times2, times1
        prim1, prim2 = prim2, prim1
        uov1, uov2 = uov2, uov1
        mom1, mom2 = mom2, mom1
        src1, src2 = src2, src1

    plot_frame(0, len(uov1), times1, times2)
