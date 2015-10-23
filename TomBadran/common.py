from pylab import *
from astropy.io import fits
import aplpy
from scipy import optimize


def gaussian2D(height, cx, cy, w, h, c):
    return lambda x, y: height * exp(-(((cx - x) / w)**2 + ((cy - y) / h)**2) / 2) + c


def fitgaussian2D(data):
    def errorfunction(p):
        return ravel(gaussian2D(*p)(*indices(data.shape)) - data)
    max_pos = unravel_index(data.argmax(), data.shape)
    params = data.max(), max_pos[0], max_pos[1], 10, 10, median(data)
    p, success = optimize.leastsq(errorfunction, params)
    return p


def nanfilter(data):
     return data[~isnan(data)]


def nansum(data):
    return nanfilter(data).sum()


def SNR(file, r1=36, r2=42, log=False):
        source = fits.open(file)
        signal = nanfilter(source[1].data).max()
        w, h = source[1].data.shape
        cx, cy = 86, 95

        bgs = []
        for i in range(0, h):
            for j in range(0, w):
                r = sqrt(((cx - j))**2 + (cy - i)**2)
                if r>= r1 and r <= r2:
                    bgs.append(source[1].data[i][j])

        bg = sqrt((array(bgs)**2).mean())

        if log:
            print(signal, bg)

        return signal / bg
