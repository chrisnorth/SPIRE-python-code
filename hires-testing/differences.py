from common import *
from astropy import wcs
from scipy.optimize import curve_fit
from scipy.ndimage.filters import convolve
from matplotlib.ticker import FormatStrFormatter
import aplpy

y_format =  FormatStrFormatter('%.2e')
fits_dir_base = "/Users/Tom/HIPE/plots/"

obsids = [1342189427, 1342188588, 1342185538]

beams = {
    'PLW' : [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 200],
    'PMW' : [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125],
    'PSW' : [5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 75, 100],
}

bands = ['PLW', 'PMW', 'PSW']

weight_r = {
    'PLW' : 12,
    'PMW' : 8,
    'PSW' : 6
}

minima_locs = {
    'PLW' : [41, 82, 123],
    'PMW' : [31, 62, 93],
    'PSW' : [22, 45, 68]
}

def diffs_for_obs(obsid):
    fits_dir = fits_dir_base + str(obsid) + '/'
    figure(figsize=(12,6))

    for band in bands:
        truth_hdulist = fits.open(fits_dir + 'TRUTH_{}.fits'.format(band))

        truth_image = truth_hdulist[1].data
        truth_wcs = wcs.WCS(truth_hdulist[1].header)

        truth_image = nan_to_num(truth_image)
        truth_peak = truth_image.max()

        weight_image = truth_image.copy()

        # Gernerate a weighting based on proximity to bright pixels
        # Create a disc as filter
        r = weight_r[band]
        y,x = np.ogrid[-r: r+1, -r: r+1]
        mask = x**2+y**2 <= r**2

        weight_image = convolve(truth_image, mask)
        weight_image[weight_image > 1] = 1
        weight_image[weight_image < 0] = 0

        weight_image = nan_to_num(weight_image)

        # imshow(weight_image)
        # show()
        # return

        #weight_image = nan_to_num(weight_image) / weight_image.max()

        if band == 'PLW':
            subplot('131')
        elif band == 'PMW':
            subplot('132')
        else:
            subplot('133')

        title(band)
        render = plot
        hires_differences = []
        hires_diff_maps = []

        for b in beams[band]:
            if b == 0:
                hdulist = fits.open(fits_dir + 'ORIGINAL_{}.fits'.format(band, b))
            else:
                hdulist = fits.open(fits_dir + 'HIRES_{}_BEAMHSIZE_{}.fits'.format(band, b))

            hires_image = hdulist[1].data
            truth_image_temp = truth_image.copy()

            # Propogate NaNs into both images so we arent summing overlapping pixels
            hires_image[isnan(truth_image_temp)] = nan
            truth_image_temp[isnan(hires_image)] = nan

            weight_image_temp = weight_image.copy()
            weight_image_temp[isnan(hires_image)] = nan

            # Normalise image to have same total flux as truth image
            flux_correction_factor = (nansum(truth_image_temp) / nansum(hires_image))
            print(flux_correction_factor)
            hires_image *= flux_correction_factor

            # Take differences and weight based on brightness in truth
            im_diff = (truth_image - hires_image) * weight_image
            hires_diff_maps.append(im_diff)
            #hires_differences.append(sqrt(nanfilter(truth_image - hires_image)**2).sum() / len(nanfilter(truth_image_temp)))
            hires_differences.append(sqrt(nanfilter(truth_image - hires_image)**2).sum() / nansum(weight_image))

        render(beams[band], hires_differences, label=band)
        gca().yaxis.set_major_formatter(y_format)

        print('\t'.join([str(x) for x in beams]))
        print('\t'.join([str(x) for x in hires_differences]))

        if band == 'PMW':
            xlabel('Beam Half Pixel Size')
        if band == 'PLW':
            ylabel('RMS Image Pixel Difference (Relative)')

        yl = ylim()

        vlines(minima_locs[band], *yl)
        ylim(yl)

    tight_layout()
    savefig('doc/beam-size-{}.pdf'.format(obsid))

if __name__ == '__main__':
    for obsid in obsids:
        diffs_for_obs(obsid)
