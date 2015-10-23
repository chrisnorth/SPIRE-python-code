from common import *
from astropy import wcs
from scipy.optimize import curve_fit
from scipy.ndimage.filters import convolve
from matplotlib.ticker import FormatStrFormatter
import aplpy

obsids = [1342189427, 1342188588, 1342185538]

beams = {
    'PLW' : [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 200],
    'PMW' : [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125],
    'PSW' : [5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 75, 100],
}

bands = ['PLW', 'PMW', 'PSW']
fits_dir_base = "/Users/Tom/HIPE/plots/"

for obsid in obsids:
    fits_dir = fits_dir_base + str(obsid) + '/'

    figure()
    plot_rows = len(beams['PLW'])
    plot_cols = len(bands)
    for k, band in enumerate(bands):
        truth_hdulist = fits.open(fits_dir + 'TRUTH_{}.fits'.format(band))

        fdata = nan_to_num(truth_hdulist[1].data)

        amax = fdata.argmax()
        y = int(floor(amax / fdata.shape[1]))

        for i, b in enumerate(beams[band]):
            subplot_index = i*plot_cols + k + 1
            #print(subplot_index, k, i)
            hdulist = fits.open(fits_dir + 'HIRES_{}_BEAMHSIZE_{}.fits'.format(band, b))
            subplot(plot_rows, plot_cols, subplot_index)
            plot(nan_to_num(hdulist[1].data[y]))
            xticks([])
            yticks([])
            if i == 0:
                tstr = band + ' ' + str(b)
            else:
                tstr = str(b)
            title(tstr, x=0.2, y=0.3)
            xs = xlim()
            xlim(xs[1] / 4, 3 * xs[1] / 4)

    tight_layout()
    subplots_adjust(hspace=0, wspace=0)
show()
