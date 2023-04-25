import numpy as np
from athena_read import athdf
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse
from pathlib import Path
from glob import glob
import pdb

plt.rcParams.update({"text.usetex": True})

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

def arad(tau, F_s, tau_s, temp=1e4):
    # Handle both array-like and single-value arguments for tau
    if isinstance(tau, (list, tuple, np.ndarray)):
        tau_term = np.zeros(len(tau))
        tau_term[tau > tau_s] = (tau[tau > tau_s] + tau_s)**(-2/3) - np.abs(tau[tau > tau_s] - tau_s)**(-2/3)
        tau_term[tau < tau_s] = (tau[tau < tau_s] + tau_s)**(-2/3) + np.abs(tau[tau < tau_s] - tau_s)**(-2/3)
        tau_term = np.concatenate([tau_term[tau > tau_s], tau_term[tau < tau_s]])
        tau = np.concatenate([tau[tau > tau_s], tau[tau < tau_s]])
    else:
        if tau > tau_s:
            tau_term = (tau + tau_s)**(-2/3) - np.abs(tau - tau_s)**(-2/3)
        elif tau < tau_s:
            tau_term = (tau + tau_s)**(-2/3) + np.abs(tau - tau_s)**(-2/3)
        else:
            raise ValueError("Undefined result for tau = tau_s.")
    a_rad = 1e-2 * F_s * (1e4 / temp) * tau_term
    return tau, a_rad

def extract_input_val(string):
    return string.split("=")[1].split("#")[0].strip()


def midpoint_diff(t):
    midpoints = 0.5*(t[1:]+t[:-1])
    dt = np.diff(midpoints)
    dt = np.insert(dt, 0, midpoints[1] - midpoints[0])
    dt = np.append(dt, midpoints[-1] - midpoints[-2])
    return dt


def plot_frame(prim, mcmom, mcsrc, uov, logpath, basename, cycle, crd="1", ctend=np.median):

    # Convert spatial variable to optical depth
    radius = prim["x1v"]
    column = uov['column']
    radial_column = np.mean(np.mean(column, axis=0), axis=0)
    sigma = line_strength / np.sqrt(np.pi) / delta
    tau = sigma * radial_column

    # reconstruct volume from cell face coordinates
    r = prim["x1f"][np.newaxis, np.newaxis, :]
    th = prim["x2f"][np.newaxis, :, np.newaxis]
    ph = prim["x3f"][:, np.newaxis, np.newaxis]
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
    er = np.sum(np.sum(mcmom["Ermc"]*vol, axis=0), axis=0) / volshells
    fr = np.sum(np.sum(mcmom["Frmc" + crd]*vol, axis=0), axis=0) / volshells

    # Ionization fraction
    ionfrac = prim["r0"]
    ionfrac_3d = ionfrac.astype(np.float)
    ionfrac = ctend(ctend(ionfrac, axis=0), axis=0).astype(np.float)

    # Density
    rho = prim["rho"]
    rho = ctend(ctend(rho, axis=0), axis=0)

    # Flux radiative acceleration data
    radaccelf = mcsrc["RadforceF" + crd] / mcmom["rho"] * vol
    radaccelf = np.sum(radaccelf, axis=(0))
    radaccelf_avg = np.sum(radaccelf, axis=(0)) / volshells

    # Scattering radiative acceleration data
    radaccels = mcsrc["RadforceS" + crd] /mcmom["rho"] * vol
    radaccels = np.sum(radaccels, axis=(0))
    radaccels_avg = np.sum(radaccels, axis=(0)) / volshells

    # dU/dr
    dUdr = midpoint_diff(er) / midpoint_diff(radius)
    accel = dUdr / 3 / rho

    # Radial velocity
    vel1 = prim["vel" + crd]
    vel1 = np.sum(np.sum(vel1*vol, axis=0), axis=0)/volshells

    # Emissivity
    # ==========

    # Load debug outputs printed in code
    test = np.loadtxt("/LyraShared/bcm2vn/athena/outputs/emissivity_outputs.txt", delimiter=',')
    test_volume = test[:, 1]
    test_emissivity = test[:, 2] * hplanck * nu0

    em = np.sum(np.sum(mcmom['Emissivity'].astype(np.float), axis=0), axis=0)
    em_novol = np.sum(np.sum(mcmom['Emissivity'].astype(np.float), axis=0), axis=0) # * (128*32*32) / 100000
    em_theory = np.sum(np.sum(4.18e-13 * (ionfrac_3d * prim['rho'] / mproton)**2 * hplanck * nu0, axis=0),axis=0)

    # Hydrogen Column
    nh = ctend(ctend(uov['column'], axis=0), axis=0)

    # Find flux and optical depth near the source for Chase's analytic solution
    # This will be where the peak of the emissivity function is
    peak_ind = np.argmax(em) + 1
    tau_source = tau[peak_ind]
    fr_source = fr[peak_ind]

    # Get theory line from Chase's solution
    tau_plot, arad_plot = arad(tau, fr_source, tau_source)
    
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, sharex=True, figsize=(12, 8))

    # Panel 1 - Accelerations
    ax1.plot(tau, 1.2668653e23 / radius**2, color="pink", lw=2, label='grav')
    ax1.plot(tau, radaccelf_avg, color="c", marker="o", ms=3, label='path length')
    ax1.plot(tau, -radaccelf_avg, color="c", marker="^", ms=3, ls="--", alpha=0.5)
    ax1.plot(tau, radaccels_avg, color="orange", marker="o", ms=3, label='scattering')
    ax1.plot(tau, -radaccels_avg, color="orange", marker="^", ms=3, ls="--", alpha=0.5)
    ax1.plot(tau, -accel, color='purple', marker='o', lw=2, ms=3, label=r'$\frac{1}{3 \rho}\frac{dU}{dr}$')
    ax1.plot(tau, accel, color='purple', marker='o', lw=2, ls='--', ms=3, alpha=0.5)
    ax1.plot(tau_plot, arad_plot, color='limegreen', marker='o', lw=2, ms=3, label='Analytic solution')
    ax1.plot(tau_plot, -arad_plot, color='limegreen', marker='o', lw=2, ls='--', ms=3, label='Analytic solution', alpha=0.5)

    ax1.axhline(2*tau_source**(-2/3), color='green', ls='--')

    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylim((1e-12, 1e4))
    ax1.set_ylabel("acceleration (cm s-2)")
    ax1.invert_xaxis()
    ax1.legend(frameon=False)

    # Panel 2 - velocity
    #ax2.plot(radius, vel1, color="k", marker="o", ms=3)
    #ax2.plot(radius, -vel1, color="r", marker="^", ms=3, ls="--", alpha=0.5)
    ax2.set_ylim((1e-4, 3e3))
    ax2.plot(tau, er*clight, color='k', marker='o', ms=3, label='c * Er')
    ax2.plot(tau, -er*clight, 'k--', marker='o', ms=3)
    ax2.plot(tau, fr, color='r', marker='o', ms=3, label='Fr')
    ax2.plot(tau, -fr, 'r--', marker='o', ms=3)
    ax2.set_ylabel("er * c, fr")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")
    ax2.legend(frameon=False)

    # Panel 3 - density
    ax3.plot(tau, rho, color="g", marker="o", ms=3, label='Density')
    #ax3.set_ylim((1e-7, 1e1))
    ax3.set_xlabel(r"$\tau$")
    ax3.set_ylabel("density (g/cm3)")
    ax3.set_xscale("log")
    ax3.set_yscale("log")
    ax3.legend(frameon=False)

    # Panel 4 - ionization fraction
    ax4.plot(tau, 1 - ionfrac, color="m", marker="o", ms=3, label='1-ion frac')
    ax4.set_xlabel(r"$\tau$")
    ax4.set_ylabel("ionization fraction")
    ax4.set_ylim((0.000001, 2))
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position("right")
    ax4.legend(frameon=False)

    # Panel 5 - emissivity
    ax5.plot(tau, em, color='limegreen', marker='o', ms=3, label='emissivity')
    ax5.plot(tau, em_theory, color='green', ls='--', label=r'$\alpha n_p^2 h \nu_0$')
    ax5.set_xlabel(r"$\tau$")
    ax5.set_ylabel("emissivity (ergs s-1)")
    ax5.set_xscale("log")
    ax5.set_yscale("log")
    ax5.legend(frameon=False)

    # Panel 6 - Velocity
    ax6.plot(tau, vel1, color="darkred", marker="o", ms=3)
    ax6.plot(tau, -vel1, color="darkred", ls='--', marker="o", ms=3)
    ax6.set_xlabel(r"$\tau$")
    ax6.set_ylabel("velocity (cm/s)")
    ax6.set_ylim((1e1, 5e6))
    ax6.set_xscale("log")
    ax6.set_yscale("log")
    ax6.yaxis.tick_right()
    ax6.yaxis.set_label_position("right")

    # Figure caption generation
    nphot, nx1, nx2, nx3, pbase, temp, accel, runtime = [np.nan] * 8
    with open(logpath, "r") as f:
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
                accel = str(extract_input_val(l) == "true")
            if "cpu time used" in l:
                runtime = float(extract_input_val(l))

    captiontext = (r"\noindent \textbf{{Number of photons:}} {:.0e}. "
                   r"\textbf{{Mesh:}} {} $\times$ {} $\times$ {}.\\ "
                   r"\textbf{{Base pressure:}} {:.2e} dyne/cm$^2$. "
                   r"\textbf{{Temperature:}} {:.0e} K.\\ "
                   r"\textbf{{Code acceleration:}} {}. \textbf{{Wall time:}}"
                   r" {:.2f} hrs.".format(
                       nphot, nx1, nx2, nx3, pbase, temp, accel,
                       runtime / 60.0 / 60.0))

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.0, bottom=0.15, wspace=0.0)
    plt.suptitle(args.title)

    ax6.text(-1, -0.25, captiontext, transform=ax6.transAxes, fontsize=12)

    #plt.savefig(basename+'.{:05d}.png'.format(cycle))
    #plt.close()
    plt.show()
    
    # Panel 6 - Velocity
    #plt.figure()
    #plt.plot(radius, vel1, color="darkred", marker="o", ms=3)
    #plt.plot(radius, -vel1, color="darkred", ls='--', marker="o", ms=3)
    #plt.xlabel(r"$r$")
    #plt.ylabel("velocity (cm/s)")
    #ax6.set_ylim((0.000001, 2))
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.show()

    time = prim['Time']
    return er, time


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("basename", type=str)
    parser.add_argument("-t", "--title", type=str)
    parser.add_argument("-d", "--dynamic", action="store_true", default=False)
    args = parser.parse_args()

    crd = "1"
    ctend = np.mean

    if args.dynamic:
        nfiles = len(glob(args.basename + ".out1.*.athdf"))
    else:
        nfiles = 1

    times = []
    ratios = []
    median_er = []

    for cycle in range(nfiles):
        cyclestr = "{:05d}.athdf".format(cycle)
        prim = athdf(
            args.basename + ".out1." + cyclestr,
            quantities=["rho", "r0", "vel" + crd, "press"],
        )
        mcmom = athdf(
            args.basename + ".out2." + cyclestr,
#            quantities=["Ermc", "Frmc" + crd, "rho", "tgas", "Emissivity", "Eavemc", "Time"],
        )
        mcsrc = athdf(
            args.basename + ".out3." + cyclestr,
            quantities=["RadforceF" + crd, "RadforceS" + crd],
        )

        uov = athdf(args.basename+".uov0."+cyclestr)
        # import pdb
        # pdb.set_trace()
        # Load run info from log
        logpath = Path(args.basename + ".txt")
        if not logpath.exists():
            logpath = Path("output.txt")



        # COMPARE ENERGY DENSITY BETWEEN TIMESTEPS
        try:
            print(eri)

            # Convert spatial variable to optical depth
            radius = prim["x1v"]

            # reconstruct volume from cell face coordinates
            r = prim["x1f"][np.newaxis, np.newaxis, :]
            th = prim["x2f"][np.newaxis, :, np.newaxis]
            ph = prim["x3f"][:, np.newaxis, np.newaxis]
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
            er = np.sum(np.sum(mcmom["Ermc"]*vol, axis=0), axis=0) / volshells

            times.append(timei)
            ratios.append(np.median(er/eri))
            median_er.append(np.median(er))
            '''
            plt.plot(radius, er, label="i")
            plt.plot(radius, eri, label="i-1")
            plt.xlabel("radius")
            plt.ylabel("energy density")
            plt.xscale("log")
            plt.yscale("log")
            plt.legend()
            plt.show()

            plt.plot(radius, er/eri, label="i / i-1")
            plt.xlabel("radius")
            plt.ylabel("energy density ratio")
            plt.xscale("log")
            plt.legend()
            plt.title("time: {}   dt: {}".format(mcmom['Time'], mcmom['Time']-timei))
            plt.show()
            '''
        except:
            pass

        eri, timei = plot_frame(prim, mcmom, mcsrc, uov, logpath, args.basename, cycle, crd=crd, ctend=ctend)

    #plt.plot(times, 1/np.array(ratios), marker='o', markersize=3)
    #plt.xlabel('time')
    #plt.ylabel('$E_i$ / $E_{i+1}$')
    #plt.show()
