from athena_read import athdf
from glob import glob
import argparse
import matplotlib.pyplot as plt

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
parser.add_argument('--basename', '-b', type=str)
args = parser.parse_args()

def plot_frame(x, y, filename, xlabel=None, ylabel=None, ycol='r', xscale='linear',
               yscale='linear', xlim=None, ylim=None, y2=None, y2label=None, 
               y2scale='linear', y2lim=None, y2col='b', time=None, **kwargs):
    fig, ax = plt.subplots(1, 1)
    ax.plot(x, y, c=ycol, **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel, color=ycol)
    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    if y2 is not None:
        ax2 = ax.twinx()
        ax2.plot(x, y2, c=y2col, **kwargs)
        ax2.set_ylabel(y2label, color=y2col)
        ax2.set_yscale(y2scale)
        ax2.set_ylim(y2lim)

    if time is not None:
        ax.text(0.1, 0.85, 't={:.1f}'.format(time), transform=ax.transAxes)
    else:
        ax.text(0.1, 0.85, 'Frame {:d}'.format(int(f.split(['out1.'])[1].split('.athdf')[0])))

    framename = filename.replace('.athdf', '.png').replace('.out1', '')
    plt.savefig(framename)
    plt.close()

basename = args.basename
files = glob(basename+'*.athdf')

for f in files:
    data = athdf(f)
    ionization_frac =  data['r0'][0,0,:]
    plot_frame(data['x1v'], data['rho'][0, 0, :], f, yscale='log', xscale='log', xlabel='r (cm)',
               ylabel='rho (g cm-3)', ylim=(1e-18, 1e-11), ycol='r', time=data['Time'],
               y2=ionization_frac, y2label='ionization fraction', y2scale='linear', y2lim=(0., 1.),
               marker='o', ms='2', alpha=0.75)

