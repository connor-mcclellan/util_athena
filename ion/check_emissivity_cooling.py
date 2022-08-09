import numpy as np
from athena_read import athdf
import pdb
import matplotlib.pyplot as plt

nphot = 1e8
ncells = 128 * 64 * 64

kb = 1.380649e-16
mp = 1.6726e-24
temp = 1e4
vth = np.sqrt( 2. * kb * temp / mp)

c = 2.99792458e10
nu0 = 2.468e15
dopw = nu0 * vth / c

x0 = 0.
h = 6.62607015e-27;
energy0 = h * (nu0 + dopw * x0)

emission = np.load("emission.npy", allow_pickle=True) * energy0
volume = np.load("volume.npy", allow_pickle=True)
data = athdf("./ion.out2.00000.athdf", quantities=["Cooling"])
cooling = -data['Cooling'] * volume * ncells

# Radial average of cooling and emission
cooling = np.average(cooling, axis=(0, 1))
emission = np.average(emission, axis=(0, 1))

fig, ax = plt.subplots()
ax.plot(data['x1v'], cooling, 'c', lw=3, label='cooling')
ax.plot(data['x1v'], emission, 'r--', lw=3, label='emission')
ax.legend()
ax.set_xlabel("r")
ax.set_ylabel("spherically-averaged energy per second")
ax.set_yscale('log')
ax.set_xscale('log')
plt.show()


