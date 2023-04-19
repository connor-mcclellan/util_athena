from athena_read import tab
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

def midpoint_diff(t):
    midpoints = 0.5 * (t[1:] + t[:-1])
    dt = np.diff(midpoints)
    dt = np.insert(dt, 0, midpoints[1] - midpoints[0])
    dt = np.append(dt, midpoints[-1] - midpoints[-2])
    return dt


def find_filenames(path):
    files = sorted(glob(path+'/*'))
    files = [str(Path(f).resolve()) for f in files if "rst" not in f and "txt" not in f]
    basename = files[0].split(".")[0]

    uovfiles = sorted(glob(basename + ".uov.*.tab"))
    primfiles = sorted(glob(basename + ".prim.*.tab"))
    momfiles = sorted(glob(basename + ".rad_mom.*.tab"))
    return uovfiles, primfiles, momfiles


def get_times(prim, istart, iend):
    times = []
    for i in range(istart, iend):
        # Load in data from each output block
        primdata = tab(prim[i])
        times.append(primdata["time"])
    return np.array(times)


def plot_frame(istart, iend, times1, times2):
    for i in range(istart, iend):

        uovdata = tab(uov1[i])
        primdata = tab(prim1[i])
        momdata = tab(mom1[i])

        this_time = primdata["time"]
        time_match_ind = np.argwhere(times2 > this_time)[0][0]
        
        uov2before = tab(uov2[time_match_ind-1])
        prim2before = tab(prim2[time_match_ind-1])
        mom2before = tab(mom2[time_match_ind-1])

        uov2after = tab(uov2[time_match_ind])
        prim2after = tab(prim2[time_match_ind])
        mom2after = tab(mom2[time_match_ind])

        # Apply a linear interpolation to the output2 data to match the time of output1


        # Get relevant quantities from data
        radius = uovdata["x1v"] * rg / 1e6
        pressure = primdata["press"]
        rho = primdata["rho"]
        mu = 1.
        # temp = primdata['temp']

        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(8, 6), sharex=True, dpi=100)

        # PLOT 1: Acceleration & Hydro quantities
        # =======================================

        # Panel 1.5: Acceleration for first output
        # -----------------------

        # Gravitational Acceleration
        grav_accel = uovdata["GravSrc_IM1"]
        ax1.plot(radius, grav_accel, **kwa, c="limegreen", label="grav accel")

        # Radiation acceleration - Prat * sigma_s * F_com / rho
        rad_accel = (prat * primdata["Sigma_s_0"] / rho * momdata["Fr01"])
        code_accel = uovdata["RadAccel_IM1"]
        ax1.plot(radius, rad_accel, **kwa, c="darkred", label="rad accel")
        ax1.plot(radius, -rad_accel, **kwaneg, c="darkred")
        ax1.plot(radius, code_accel, **kwa, c="c", label="code rad accel")
        ax1.plot(radius, -code_accel, **kwaneg, c="c")

        # Pressure gradient acceleration
        grad_P = midpoint_diff(pressure) / midpoint_diff(radius / rg) / rho
        ax1.plot(radius, np.abs(grad_P), **kwa, c="purple", label=r"$|\nabla P|$")

        # Panel formatting
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylim((1e-9, 1e9))
        ax1.set_xlim((1.1999, 1.2015))
        ax1.grid(which='both', alpha=0.75, lw=0.3)
        ax1.set_ylabel("a [dyne cm${^-3}$]")
        ax1.legend(**lega)
        ax1.set_title((prim1[i]).split("/")[-1])
        ax1.set_xlabel("radius [cm] / 1e6")

        # Panel 1.2: Acceleration for second output
        # -----------------------

        # Gravitational Acceleration
        grav_accel = uov2before["GravSrc_IM1"]
        ax2.plot(radius, grav_accel, **kwa, c="limegreen", label="grav accel")

        # Radiation acceleration - Prat * sigma_s * F_com / rho
        rad_accel = (prat * prim2before["Sigma_s_0"] / rho * mom2before["Fr01"])
        code_accel = uov2before["RadAccel_IM1"]
        ax2.plot(radius, rad_accel, **kwa, c="darkred", label="rad accel")
        ax2.plot(radius, -rad_accel, **kwaneg, c="darkred")
        ax2.plot(radius, code_accel, **kwa, c="c", label="code rad accel")
        ax2.plot(radius, -code_accel, **kwaneg, c="c")

        # Pressure gradient acceleration
        grad_P = midpoint_diff(pressure) / midpoint_diff(radius / rg) / rho
        ax2.plot(radius, np.abs(grad_P), **kwa, c="purple", label=r"$|\nabla P|$")

        # Panel formatting
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_ylim((1e-9, 1e9))
        ax2.grid(which='both', alpha=0.75, lw=0.3)
        ax2.set_ylabel("a [dyne cm${^-3}$]")
        ax2.yaxis.tick_right()
        ax2.legend(**lega)
        ax2.set_title((prim2[time_match_ind-1]).split("/")[-1])
        ax2.set_xlabel("radius [cm]")

        # Panel 1.3: Radiation moments for first output
        # ----------------------------
        # Energy density and flux
        er = arad * temp0**4 * primdata['Er']
        fr = c.c.cgs * arad * temp0**4 * primdata['Fr1']

        print("flux0: {:.5e}".format(primdata['Fr1'][2]))
        ax3.plot(radius, c.c.cgs * er, **kwa, c='k', label=r'$cE_r$')
        ax3.plot(radius, -c.c.cgs * er, **kwaneg, c='k')
        ax3.plot(radius, fr, **kwa, c='r', label='$F_r$')
        ax3.plot(radius, -fr, **kwaneg, c='r')
        ax3.plot(radius, (c.c.cgs * c.G.cgs * ns_mass / kes / radius**2).cgs, **kwa, c='limegreen', label=r'$F_{\rm Edd}$')
        ax3.xaxis.set_major_formatter(ScalarFormatter())

        # Panel formatting
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_ylim((1e23, 1e31))
        ax3.grid(which='both', alpha=0.75, lw=0.3)
        ax3.set_ylabel("[erg cm$^{-2}$ s$^{-1}$]")
        ax3.yaxis.set_label_position("right")
        ax3.legend(**lega)
        ax3.set_xlabel("radius [cm] / 1e6")
        ax3.tick_params(axis='x', which='minor', labelsize=7)

        captiontext = "t={:.4e}".format(primdata["time"])
        ax3.text(0, -.25, captiontext, transform=ax3.transAxes, fontsize=12)

        # Panel 1.4: Radiation moments for second output
        # ----------------------------
        # Energy density and flux
        er = arad * temp0**4 * prim2before['Er']
        fr = c.c.cgs * arad * temp0**4 * prim2before['Fr1']

        print("flux0: {:.5e}".format(prim2before['Fr1'][2]))
        ax4.plot(radius, c.c.cgs * er, **kwa, c='k', label=r'$cE_r$')
        ax4.plot(radius, -c.c.cgs * er, **kwaneg, c='k')
        ax4.plot(radius, fr, **kwa, c='r', label='$F_r$')
        ax4.plot(radius, -fr, **kwaneg, c='r')
        ax4.plot(radius, (c.c.cgs * c.G.cgs * ns_mass / kes / radius**2).cgs, **kwa, c='limegreen', label=r'$F_{\rm Edd}$')

        # Panel formatting
        ax4.set_xscale("log")
        ax4.set_yscale("log")
        ax4.set_ylim((1e23, 1e31))
        ax4.grid(which='both', alpha=0.75, lw=0.3)
        ax4.set_ylabel("[erg cm$^{-2}$ s$^{-1}$]")
        ax4.yaxis.set_label_position("right")
        ax4.yaxis.tick_right()
        ax4.legend(**lega)
        ax4.set_xlabel("radius [cm] / 1e6")
        ax4.xaxis.set_major_formatter(ScalarFormatter())
        ax4.xaxis.set_minor_formatter(ScalarFormatter())
        ax4.tick_params(axis='x', which='minor', labelsize=7)
        

        captiontext = "t={:.4e}".format(prim2before["time"])
        ax4.text(0, -0.25, captiontext, transform=ax4.transAxes, fontsize=12)
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0, bottom=0.15, wspace=0.0)
        filename = uov1[i]
        filename2 = (uov2[time_match_ind-1]).split("/")[-1].replace(".uov", "").replace(".tab", "")
        framename = filename.replace(".tab", "_"+filename2+"_CLOSEUP.png").replace(".uov", "")

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

    uov1, prim1, mom1 = find_filenames(args.path1)
    uov2, prim2, mom2 = find_filenames(args.path2)

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

    plot_frame(0, len(uov1), times1, times2)
