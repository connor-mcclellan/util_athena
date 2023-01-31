import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure(figsize=(4.5, 4))
ax = fig.add_subplot(111, projection='3d')
import numpy as np
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.collections as mcoll
import matplotlib.path as mpath

def loadout(fname, **kwargs):
    try:
        f = np.loadtxt(fname, skiprows=4, **kwargs)
    except Exception as e:
        rowfail = int(str(e).split('line ')[-1])
        f = np.loadtxt(fname, skiprows=4, max_rows=rowfail-5, **kwargs)
    return f


def colorline(x, y, z, k=None, cmap=plt.get_cmap('bwr_r'), norm=plt.Normalize(0.0, 1.0), **kwargs):

    # Default colors equally spaced on [0,1]:
    if k is None:
        k = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(k, "__iter__"):  # to check for numerical input -- this is a hack
        k = np.array([k])

    k = np.asarray(k)

    segments = make_segments(x, y, z)
    lc = Line3DCollection(segments, array=k[1:], cmap=cmap, norm=norm, **kwargs)

    ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_segments(x, y, z):
    points = np.array([x, y, z]).T.reshape(-1, 1, 3)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments




R = 1e11

#x, y, z, accel = np.loadtxt("output.txt", skiprows=4, max_rows=144420, unpack=True)

x_ip, y_ip, z_ip, accel_ip = loadout("output_ip_acc100000.txt", unpack=True)
##x_ip_na, y_ip_na, z_ip_na, accel_ip_na = loadout("output_ip_noacc.txt", unpack=True)
#colorline(x_ip, y_ip, z_ip, k=accel_ip, lw=0.2)

#ax.set_xlim((-R, R))
#ax.set_ylim((-R, R))
#ax.set_zlim((-R, R))
#ax.plot(x, y, z, lw=0.1, c='b', label='Accelerated')
#ax.plot(x_na, y_na, z_na, lw=0.1, c='r', label='Not accelerated')
#ax.plot(x_na_nv, y_na_nv, z_na_nv, lw=0.1, c='r', label='Not accelerated')
#ax.plot(x_nv, y_nv, z_nv, lw=0.1, c='b', label='Accelerated')

ax.plot(x_ip, y_ip, z_ip, lw=0.1, c='b', label='Accelerated')
#ax.plot(x_ip_na, y_ip_na, z_ip_na, lw=0.1, c='r', label='Not accelerated')

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_ylabel("z")

def wireframesph(r, nth, nph):
    # Make data
    u = np.linspace(0, 2 * np.pi, nph*4)
    v = np.linspace(0, np.pi, nth*4)
    x = r * np.outer(np.cos(u), np.sin(v))
    y = r * np.outer(np.sin(u), np.sin(v))
    z = r * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(x, y, z, color='k', rstride=4, cstride=4, alpha=0.005 + 0.02*(r/R)**5, linewidth=0.1)

for r in np.linspace(0, R, 17):
    wireframesph(r, 8, 8)

# Make data for sphere
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = R * np.outer(np.cos(u), np.sin(v))
y = R * np.outer(np.sin(u), np.sin(v))
z = R * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
#ax.plot_surface(x, y, z, color='k', alpha=0.05)
ax.auto_scale_xyz(x, y, z)
ax.grid(False)
ax.tick_params(axis='both', which='major', labelsize=4)
ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
plt.legend()
plt.show()
