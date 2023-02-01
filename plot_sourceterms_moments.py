import numpy as np
from athena_read import athdf
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse
from pathlib import Path
from glob import glob

plt.rcParams.update({"text.usetex": True})


def extract_input_val(string):
    return string.split("=")[1].split("#")[0].strip()


def plot_frame(
    prim, mcmom, mcsrc, logpath, basename, cycle, crd="1", ctend=np.median
):

    radius = prim["x1v"]

    # reconstruct volume from cell face coordinates
    r = prim["x1f"][np.newaxis, np.newaxis, :]
    th = prim["x2f"][np.newaxis, :, np.newaxis]
    ph = prim["x3f"][:, np.newaxis, np.newaxis]
    vol = (
        1 / 3
        * (r[:, :, 1:] ** 3 - r[:, :, :-1] ** 3)
        * (np.cos(th[:, -1:, :]) - np.cos(th[:, 1:, :]))
        * (ph[1:, :, :] - ph[-1:, :, :])
    )

    # Flux radiative acceleration data
    radaccelf = mcsrc["RadforceF" + crd] / mcmom["rho"] * vol
    radaccelf = ctend(radaccelf, axis=(0))
    radaccelf = radaccelf[1:]
    radaccelf_avg = ctend(radaccelf, axis=(0))

    # Scattering radiative acceleration data
    radaccels = mcsrc["RadforceS" + crd] / mcmom["rho"] * vol
    radaccels = ctend(radaccels, axis=(0))
    radaccels = radaccels[1:]
    radaccels_avg = ctend(radaccels, axis=(0))

    # Density
    rho = prim["rho"]
    rho = ctend(ctend(rho, axis=0), axis=0)

    # Ionization fraction
    ionfrac = prim["r0"] / mcmom["rho"] * prim["rho"]
    ionfrac = ctend(ctend(ionfrac, axis=0), axis=0)

    # Radial velocity
    vel1 = prim["vel" + crd]
    vel1 = ctend(ctend(vel1, axis=0), axis=0)

    er = ctend(ctend(mcmom["Ermc"], axis=0), axis=0) * 2.99792458e10 * vol
    fr = ctend(ctend(mcmom["Frmc" + crd], axis=0), axis=0) * vol

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        2, 2, sharex=True, figsize=(12, 8)
    )

    # Panel 1 - Accelerations
    ax1.plot(radius, 1.2668653e23 / radius**2, color="pink", lw=2)
    ax1.plot(radius, radaccelf_avg, color="c", marker="o", ms=3)
    ax1.plot(
        radius, -radaccelf_avg, color="c", marker="^", ms=3, ls="--", alpha=0.5
    )
    ax1.plot(radius, radaccels_avg, color="orange", marker="o", ms=3)
    ax1.plot(
        radius,
        -radaccels_avg,
        color="orange",
        marker="^",
        ms=3,
        ls="--",
        alpha=0.5,
    )
    ax1.set_yscale("log")
    ax1.set_xscale("log")
    ax1.set_ylim((1e-8, 1e4))
    ax1.set_ylabel("acceleration (cm s-2)")

    # Panel 2 - velocity
    ax2.plot(radius, vel1, color="k", marker="o", ms=3)
    ax2.plot(radius, -vel1, color="r", marker="^", ms=3, ls="--", alpha=0.5)
    ax2.set_ylim((1e-2, 1e8))
    ax2.set_ylabel("velocity (cm/s)")
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.yaxis.tick_right()
    ax2.yaxis.set_label_position("right")

    # Panel 3 - density
    ax3.plot(radius, rho, color="g", marker="o", ms=3)
    ax3.set_ylim((1e-7, 1e1))
    ax3.set_xlabel("r")
    ax3.set_ylabel("density (g/cm3)")
    ax3.set_xscale("log")
    ax3.set_yscale("log")

    # Panel 4 - ionization fraction
    ax4.plot(radius, 1 - ionfrac, color="m", marker="o", ms=3)
    ax4.set_xlabel("r")
    ax4.set_ylabel("ionization fraction")
    ax4.set_ylim((0.000001, 2))
    ax4.set_xscale("log")
    ax4.set_yscale("log")
    ax4.yaxis.tick_right()
    ax4.yaxis.set_label_position("right")

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

    # Panel 1 legend
    formatlegend = [
        Line2D(
            [1],
            [0],
            color="pink",
            linestyle="-",
            label=r"Gravitational acceleration",
        ),
        Line2D(
            [1],
            [0],
            color="c",
            marker="o",
            linestyle="None",
            label=r"Path-length radiative acceleration",
        ),
        Line2D(
            [1],
            [0],
            color="orange",
            marker="o",
            linestyle="None",
            label=r"Scattering radiative acceleration",
        ),
        Line2D(
            [1],
            [0],
            color="k",
            marker="o",
            linestyle="-",
            label="Median value",
        ),
        Line2D(
            [1], [0], color="k", marker="o", linestyle="--", label="RMS value"
        ),
    ]

    fmtlegend = ax1.legend(
        handles=formatlegend,
        loc="upper right",
        bbox_to_anchor=[1.0, 0.98],
        frameon=False,
    )

    # Panel 3 legend
    formatlegend = [
        Line2D(
            [1], [0], color="g", marker="o", linestyle="None", label=r"Density"
        )
    ]
    fmtlegend = ax3.legend(
        handles=formatlegend,
        loc="upper right",
        bbox_to_anchor=[1.0, 0.98],
        frameon=False,
    )

    # Panel 4 legend
    formatlegend = [
        Line2D(
            [1],
            [0],
            color="m",
            marker="o",
            linestyle="None",
            label=r"1 - Ionization fraction",
        )
    ]
    fmtlegend = ax4.legend(
        handles=formatlegend,
        loc="upper right",
        bbox_to_anchor=[1.0, 0.98],
        frameon=False,
    )

    ax4.text(-1, -1.15, captiontext, transform=ax2.transAxes, fontsize=12)
    if crd == "1":
        ax1.set_title("acceleration in x1 direction")
    elif crd == "2":
        ax1.set_title("acceleration in x2 direction")
    elif crd == "3":
        ax1.set_title("acceleration in x3 direction")

    #    plt.savefig(basename+'.{:05d}.png'.format(cycle))
    #    plt.close()
    plt.show()


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

    for cycle in range(nfiles):
        cyclestr = "{:05d}.athdf".format(cycle)
        prim = athdf(
            args.basename + ".out1." + cyclestr,
            quantities=["rho", "r0", "vel" + crd, "press"],
        )
        mcmom = athdf(
            args.basename + ".out2." + cyclestr,
            quantities=["Ermc", "Frmc" + crd, "rho", "tgas"],
        )
        mcsrc = athdf(
            args.basename + ".out3." + cyclestr,
            quantities=["RadforceF" + crd, "RadforceS" + crd],
        )
        # import pdb
        # pdb.set_trace()
        # Load run info from log
        logpath = Path(args.basename + ".txt")
        if not logpath.exists():
            logpath = Path("output.txt")

        plot_frame(
            prim,
            mcmom,
            mcsrc,
            logpath,
            args.basename,
            cycle,
            crd=crd,
            ctend=ctend,
        )
