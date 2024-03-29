import numpy as np
from athena_read import athdf
import pdb
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex": True})
import argparse

GM = 1.2668653e23
mass = 1.660538782e-24
kb = 1.380649e-16
T = 1e4

isoth_sound = np.sqrt(kb * T / mass)
r_sonic = GM / 2. / isoth_sound**2.
print("SONIC POINT: {:.3e} CM".format(r_sonic))


parser = argparse.ArgumentParser()
parser.add_argument("basename", type=str)
parser.add_argument("-t", "--title", type=str)
args = parser.parse_args()



prim = athdf(args.basename+'.out1.00000.athdf', quantities=["rho"])
#mcmom = athdf(args.basename+'.out2.00000.athdf', quantities=[""])
mcsrc = athdf(args.basename+'.out3.00000.athdf', quantities=["RadforceF1", "RadforceS1"])
radius = prim['x1v']

# TODO: Temporary flip after Shane's updates to moments
radaccelf = mcsrc['RadforceS1'] / prim['rho']
radaccelf = np.average(radaccelf, axis=(0, 1))

radaccels = mcsrc['RadforceF1'] / prim['rho']
radaccels = np.average(radaccels, axis=(0, 1))

surf_g = GM / radius**2

fig, ax = plt.subplots()
ax.plot(radius, radaccelf, 'c', label='flux radiative acceleration', marker='o', ms=3)
ax.plot(radius, radaccels, 'orange', label='scattering radiative acceleration', marker='o', ms=3)
ax.set_xlabel('r')
ax.set_xscale('log')
ax.set_ylabel('acceleration (cm s-2)')
ax.set_yscale('log')
plt.title(args.title)

plt.legend()
plt.show()
