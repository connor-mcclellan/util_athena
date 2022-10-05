import numpy as np
from athena_read import athdf
import pdb
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"text.usetex": True})
import argparse

def extract_input_val(string):
    return string.split("=")[1].split("#")[0].strip()

parser = argparse.ArgumentParser()
parser.add_argument("basename", type=str)
parser.add_argument("-t", "--title", type=str)
parser.add_argument("-d", "--dynamic", action='store_true', default=False)
args = parser.parse_args()

if args.dynamic:
    cyclestr = "00001.athdf"
else:
    cyclestr = "00000.athdf"

prim = athdf(args.basename+'.out1.'+cyclestr, quantities=["rho"])
mcmom = athdf(args.basename+'.out2.'+cyclestr, quantities=["Ermc", "Frmc1"])
mcsrc = athdf(args.basename+'.out3.'+cyclestr, quantities=["RadforceF1", "RadforceS1"])
radius = prim['x1v']

# Flux radiative acceleration data
radaccelf = mcsrc['RadforceF1'] / prim['rho']
radaccelf = np.average(radaccelf, axis=(0))
radaccelf_avg = np.average(radaccelf, axis=(0))

# Scattering radiative acceleration data
radaccels = mcsrc['RadforceS1'] / prim['rho']
radaccels = np.average(radaccels, axis=(0))
radaccels_avg = np.average(radaccels, axis=(0))

# Flux and energy density
# do volume average - x1f gives you faces, construct exact volume of every cell and then 
# weight every cell by its local volume and then average
er = np.average(mcmom['Ermc'],axis=(0))*2.99792458e10 
fr = np.average(mcmom['Frmc1'],axis=(0))

fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 8))
array_dims = np.shape(radaccels)
for j in range(array_dims[0]):
    ax1.scatter(radius, radaccelf[j, :], color='c', marker='o', s=3, alpha=0.25)
    ax1.scatter(radius, radaccels[j, :], color='orange', marker='o', s=3, alpha=0.25)
    ax2.scatter(radius, er[j, :], color='k', marker='o', s=3, alpha=0.25)
    ax2.scatter(radius, fr[j, :], color='r', marker='o', s=3, alpha=0.25)

ax2.set_xlabel('r')
ax2.set_ylabel(r'$F_r$, $cE_r$')
ax2.set_xscale('log')
ax2.set_yscale('log')

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_ylabel('acceleration (cm s-2)')
ax1.set_ylim((1e-6, 1e1))

# Load run info from log
with open(args.basename+'.txt', 'r') as f:
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
plt.subplots_adjust(hspace=0.0, bottom=0.15)
plt.suptitle(args.title)


formatlegend = [Line2D([1], [0], color='c', marker='o', linestyle="None", label=r'Flux radiative acceleration'), Line2D([1], [0], color='orange', marker='o', linestyle="None", label=r'Scattering radiative acceleration')]
fmtlegend = ax1.legend(handles=formatlegend, loc='upper right', bbox_to_anchor=[1.0, 0.95], frameon=False)
formatlegend = [Line2D([1], [0], color='k', marker='o', linestyle="None", label=r'$c E_r$'), Line2D([1], [0], color='r', marker='o', linestyle="None", label=r'$F_r$')]
fmtlegend = ax2.legend(handles=formatlegend, loc='upper right', bbox_to_anchor=[1.0, 0.95], frameon=False)
ax2.text(0., -0.2, captiontext, transform=ax2.transAxes, fontsize=12)
plt.show()

