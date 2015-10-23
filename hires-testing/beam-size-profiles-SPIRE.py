from common import *
from astropy import wcs
from scipy.optimize import curve_fit
from scipy.ndimage.filters import convolve
from matplotlib.widgets import Slider, RadioButtons
import aplpy

######### Settings ###################

# Directory to find files
fits_dir = "/Users/Tom/HIPE/plots/SPIRE/"

# Observation id to use
obsid = 1342216940

# Band to plot
band = 'PLW' #, 'PMW', 'PSW'

# Beam sizes to plot (this should be the same as in beam-size-SPIRE.py, 0 uses nominal map)
beam_sizes = {
    'PLW' : [0, 20, 50, 80, 100, 150],
    'PMW' : [0, 15, 35, 60, 80, 100],
    'PSW' : [0, 10, 25, 40, 60, 80],
}

fig = figure()
plots = {}

styles = list(reversed(['-', '--', '-.', ':', '-.']))
render = plot

def update(_):
    for b in beam_sizes[band]:
        if b == 0:
            hdulist = fits.open(fits_dir + '{}_NOMINAL_REGRID_{}.fits'.format(obsid, band))
        else:
            hdulist = fits.open(fits_dir + '{}_HIRES_{}_BEAMHSIZE_{}.fits'.format(obsid, band, b))
        data = nan_to_num(hdulist[1].data)
        plots[b].set_ydata( (data[int(sfreq.val)])[:len(plots[b].get_ydata())] )
    draw_idle()

for i, b in enumerate(beam_sizes[band]):
    if b == 0:
        hdulist = fits.open(fits_dir + '{}_NOMINAL_REGRID_{}.fits'.format(obsid, band))
    else:
        hdulist = fits.open(fits_dir + '{}_HIRES_{}_BEAMHSIZE_{}.fits'.format(obsid, band, b))
    data = nan_to_num(hdulist[1].data)
    h = data.shape[1]
    if b == 0:
        plots[b],  = render(data[h // 2], label='Nominal', linewidth=1, color='k')
    else:
        plots[b],  = render(data[h // 2], label=b, linewidth=1, linestyle=styles[i % len(styles)])

xlim(h / 4, 3 * h / 4)
ylim(0, 150)
ylabel('Flux MJy/sr')
xlabel('X Pixel Coordinate')
legend(loc=1)

sfreq = Slider(axes([0.25, 0.85, 0.45, 0.02]), 'Y Pixel Coordinate', 0, h, valinit=h // 2, valfmt='%d')
sfreq.on_changed(update)

show()
