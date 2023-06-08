from athena_read import tab
from glob import glob
import argparse
import matplotlib.pyplot as plt
#from matplotlib import rc
import numpy as np
import astropy.constants as c
import astropy.units as u
import multiprocessing as mp
import pdb
from pathlib import Path

#rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
#rc("text", usetex=True)

# ==========================================================================
#                                                                          -
# Use the following command to produce a .mp4 from this file's outputs:    -
#                                                                          -
# > ffmpeg -framerate 12 -i [basename].%05d.png out.mp4                    -
#                                                                          -
# Replace [basename] with your files' basenames.                           -
#                                                                          -
# ==========================================================================


parser = argparse.ArgumentParser()
parser.add_argument("path", type=str)
parser.add_argument("--start", "-s", type=int, default=0)
parser.add_argument("--show", "-w", action='store_true', default=False)
parser.add_argument("--cores", "-c", type=int, default=1)
args = parser.parse_args()


def midpoint_diff(t):
    midpoints = 0.5 * (t[1:] + t[:-1])
    dt = np.diff(midpoints)
    dt = np.insert(dt, 0, midpoints[1] - midpoints[0])
    dt = np.append(dt, midpoints[-1] - midpoints[-2])
    return dt


# Hydro outputs
basename = glob(str(Path(args.path).resolve())+"/*")[0].split(".")[0]
uovfiles = sorted(glob(basename + ".block0.uov.*.tab"))
primfiles = sorted(glob(basename + ".block0.prim.*.tab"))
momfiles = sorted(glob(basename + ".block0.rad_mom.*.tab"))
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


# For each timestep of the outputs
def plot_frame(istart, iend):
    for i in range(istart, iend):
        # Load in data from each output block
        uovdata = tab(uovfiles[i])
        primdata = tab(primfiles[i])
        momdata = tab(momfiles[i])

        # Get relevant quantities from data
        radius = uovdata["x1v"] * rg
        pressure = primdata["press"]
        rho = primdata["rho"]
        mu = 1.
        # temp = primdata['temp']

        temp0 = c.c.cgs**2 * mu * c.m_p.cgs / c.k_B.cgs

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(8, 6), sharex=True, dpi=100)

        # PLOT 1: Acceleration & Hydro quantities
        # =======================================

        # Panel 1.1: Density
        # ------------------
        ax1.plot(radius, rho, **kwa, c='darkblue', label="density")

        # Panel formatting
        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1.set_ylabel(r'$\rho$ [g cm$^{-3}$]')
        ax1.set_ylim((1e-2, 1e8))
        ax1.grid(which='both', alpha=0.75, lw=0.3)
        ax1.legend(**lega)
        #ax1.set_xlim(1.19997e6, 1.2005e6)

        # Panel 1.2: Radial Velocity and Sound Speed
        # ------------------------------------------
        cs = np.sqrt(pressure / rho)
        ax2.plot(radius, primdata["vel1"], **kwa, c='k', label="velocity")
        ax2.plot(radius, -primdata["vel1"], **kwaneg, c='k')
        ax2.plot(radius, cs, **kwa, c="orange", label="sound speed")
        ax2.plot(radius, np.ones(len(radius)) * np.sqrt(5e11 / temp0), **kwa, c='c', label='isothermal sound speed')

        # Panel formatting
        ax2.set_yscale("log")
        ax2.set_ylim((1e-9, 10))
        ax2.set_ylabel(r'$v / c$')
        ax2.yaxis.set_label_position("right")
        ax2.yaxis.tick_right()
        ax2.grid(which='both', alpha=0.75, lw=0.3)
        ax2.legend(**lega)

        # Panel 1.3: Pressure
        # -------------------
        ax3.plot(radius, pressure, **kwa, c="m", label="pressure")

        # Panel formatting
        ax3.set_yscale("log")
        ax3.set_ylim((1e-6, 1e7))
        ax3.set_ylabel('P [dyne cm$^{-2}$]')
        ax3.grid(which='both', alpha=0.75, lw=0.3)
        ax3.legend(**lega)

        # Panel 1.4: Temperature
        # ----------------------

        arad = 4 * c.sigma_sb.cgs / c.c.cgs
        ax4.plot(radius, temp0 * pressure / rho, label="temperature", **kwa, c="teal")
        ax4.plot(radius, temp0 * (primdata['Er'])**(0.25), label="radiation temperature", **kwa, c="darkred")

        # Panel formatting
        ax4.set_xscale("log")
        ax4.set_yscale("log")
        ax4.set_ylim((1e2, 1e15))
        ax4.grid(which='both', alpha=0.75, lw=0.3)
        ax4.set_ylabel('T [K]')
        ax4.yaxis.set_label_position('right')
        ax4.yaxis.tick_right()
        ax4.legend(**lega)


        # Panel 1.5: Acceleration
        # -----------------------

        # Gravitational Acceleration
        grav_accel = uovdata["GravSrc_IM1"]
        ax5.plot(radius, grav_accel, **kwa, c="limegreen", label="grav accel")

        # Radiation acceleration - Prat * sigma_s * F_com / rho
        rad_accel = (prat * primdata["Sigma_s_0"] / rho * primdata["Fr01"])
        code_accel = uovdata["RadAccel_IM1"] / rho
#        ax5.plot(radius, prat*primdata["Sigma_s_0"]/rho)
        ax5.plot(radius, rad_accel, **kwa, c="darkred", label="rad accel")
        ax5.plot(radius, -rad_accel, **kwaneg, c="darkred")
        ax5.plot(radius, code_accel, **kwa, c="c", label="code rad accel")
        ax5.plot(radius, -code_accel, **kwaneg, c="c")

        # Pressure gradient acceleration
#        grad_P = -midpoint_diff(pressure) / midpoint_diff(radius / rg) / rho
#        ax5.plot(radius, grad_P, **kwa, c="purple", label=r"$\nabla P$")
#        ax5.plot(radius, -grad_P, **kwaneg, c="purple")


        # Panel formatting
        ax5.set_xscale("log")
        ax5.set_yscale("log")
        ax5.set_ylim((1e-8, 1e5))
        ax5.grid(which='both', alpha=0.75, lw=0.3)
        ax5.set_ylabel("a [dyne cm${^-3}$]")
        ax5.legend(**lega)
        ax5.tick_params(axis='x', which="minor", labelsize=7)

        ax5.set_xlabel("radius [cm]")

        # Panel 1.6: Radiation moments
        # ----------------------------
        # Energy density and flux
        er = arad * temp0**4 * primdata['Er']
        fr0 = c.c.cgs * arad * temp0**4 * primdata['Fr01']
        fr = c.c.cgs * arad * temp0**4 * primdata['Fr1']

        ax6.plot(radius, c.c.cgs * er, **kwa, c='k', label=r'$cE_r$')
        ax6.plot(radius, -c.c.cgs * er, **kwaneg, c='k')
        ax6.plot(radius, fr0, **kwa, c='c', label='$F_{r0}$')
        ax6.plot(radius, -fr0, **kwaneg, c='c')
        ax6.plot(radius, fr, **kwa, c='r', label='$F_r$')
        ax6.plot(radius, -fr, **kwaneg, c='r')
        ax6.plot(radius, (c.c.cgs * c.G.cgs * ns_mass / kes / radius**2).cgs, **kwa, c='limegreen', label=r'$F_{\rm Edd}$')

        # Panel formatting
        ax6.set_xscale("log")
        ax6.set_yscale("log")
        ax6.set_ylim((1e17, 1e37))
        ax6.grid(which='both', alpha=0.75, lw=0.3)
        ax6.set_ylabel("[erg cm$^{-2}$ s$^{-1}$]")
        ax6.yaxis.set_label_position("right")
        ax6.yaxis.tick_right()
        ax6.legend(**lega)
        ax6.set_xlabel("radius [cm]")
        ax6.tick_params(axis='x', which="minor", labelsize=7)
        captiontext = "t={:.4e}".format(primdata["time"])
        plt.tight_layout()
        plt.subplots_adjust(hspace=0.0, bottom=0.15, wspace=0.0)
        ax6.text(-1, -.5, captiontext, transform=ax6.transAxes, fontsize=12)
        filename = uovfiles[i]
        framename = filename.replace(".tab", "_CLOSEUP.png").replace(".uov", "")
        if not args.show:
            plt.savefig(framename, dpi=150)
            plt.close()
        else:
            plt.show()


if __name__ == "__main__":

    if args.cores > 1:
        batch = int(len(uovfiles)/args.cores)
        pool = mp.Pool(processes=args.cores)
        for j in range(args.cores):
            istart = int(j*batch)
            iend = int((j+1)*batch)
            if j == args.cores-1:
                iend = len(uovfiles)
            print(istart, iend)
            pool.apply_async(plot_frame, args=(istart, iend))
        pool.close()
        pool.join()
    else:
        plot_frame(0, len(uovfiles))
