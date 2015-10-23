import sys
print(sys.version)

from pylab import *
from astropy.io import fits
from matplotlib.widgets import Slider, RadioButtons
import glob
import re

# Settings
fits_dir = "/Users/Tom/HIPE/plots/iter_test/"

obsids = [1342249237, 1342227726, 1342210936, 1342216940]
obsid = obsids[0]

band = 'PLW' #, 'PMW', 'PSW'

fig = figure()
plots = {}

styles = list(reversed(['-', '--', '-.', ':']))
render = plot

def update(_):
    hdulist = fits.open('/Users/Tom/HIPE/plots/spire/1342249237_NOMINAL_REGRID_PLW.fits')
    data = nan_to_num(hdulist[1].data)
    plots[0].set_ydata((data[int(sfreq.val)])[:len(plots[0].get_ydata())] )
    for n in files:
        hdulist = fits.open(files[n])
        data = nan_to_num(hdulist[1].data)
        plots[n].set_ydata( (data[int(sfreq.val)])[:len(plots[n].get_ydata())] )

files_orig = {}
for fname in glob.glob(fits_dir + str(obsid) + '_' + band + '*.fits'):
    n = int(re.findall('\d+', fname.split('_')[-1])[0])
    files_orig[n] = fname

plot_iters = [5, 10, 15, 20, 30, 50, 80]
files = {k: v for k, v in files_orig.items() if k in plot_iters}

hdulist = fits.open('/Users/Tom/HIPE/plots/spire/%d_NOMINAL_REGRID_PLW.fits' % obsid)
data = nan_to_num(hdulist[1].data)
h = data.shape[1]
plots[0], = render(data[h // 2], label='Nominal', linewidth=1, color='k')

for i, n in enumerate(sort([k for k, v in files.items()])):
    hdulist = fits.open(files[n])
    data = nan_to_num(hdulist[1].data)
    h = data.shape[1]
    plots[n], = render(data[h // 2], label=n, linewidth=1, linestyle=styles[i % len(styles)])

xlim(h / 4, 3 * h / 4)
ylim(0, 150)
ylabel('Flux MJy/sr')
xlabel('X Pixel Coordinate')
legend(loc=1)

sfreq = Slider(axes([0.25, 0.85, 0.45, 0.02]), 'Y Pixel Coordinate', 0, h, valinit=h // 2, valfmt='%d')
sfreq.on_changed(update)

figure()
n = sort([k for k, v in files_orig.items()])

hdulist = fits.open('/Users/Tom/HIPE/plots/spire/%d_NOMINAL_REGRID_PLW.fits' % obsid)
max_image = nan_to_num(hdulist[1].data)

diffs = zeros_like(n)

for i in n:
    hdulist = fits.open(files_orig[i])
    data = nan_to_num(hdulist[1].data)
    diffs[i-1] = sqrt(((max_image - data)**2).mean())

ylabel('RMS Pixel Difference from nominal')
xlabel('Number of Iterations')

plot(n, diffs)
plot(n, diffs, 'kx', markersize=8)
savefig('doc/iter-rms-%d.pdf' % obsid)

show()
