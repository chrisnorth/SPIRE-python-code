# Program to crop and recenter PSF images
# written by Matthew Smith Apr 2013

# import modules
import os
import numpy
import math
import scipy
import sys
from os.path import join as pj
import pyfits
from scipy.optimize import fmin

# select folder file files

#folder="/data/Herschel/Calibration/Inputs/Mars_50004e4b"
#folder="/data/Herschel/Calibration/Inputs/Mars_5000532c"
#folder="/data/Herschel/Calibration/Inputs/Mars_5000532c/1arcsec"
folder="/data/Herschel/Calibration/Inputs/Mars_50004e4b/1arcsec"

# select final size of image (arcsec)
sizeImage = 1400

###################################################################################

# choose bands
bands = ["PSW", "PMW", "PLW"]

# choose cut level to fit gaussian
cutLevel = {"PSW":240, "PMW":170, "PLW":100}

# initial beam width guess
FWHM = {"PSW":18.1, "PMW":25.2, "PLW":36.9}

###################################################################################

# define function
def gaussFitter(p, *args):
    # function to fit a gaussian model to data
    
    # the model is f = p[0] * e^(-((x-p[1])^2/(2*p[3]^2) + (y-p[2])^2/(2*p[3]^2)) + p[4]
    
    # extract arguments
    image = args[0]
    xCo = args[1]
    yCo = args[2]
        
    # create blank map
    modelImage = numpy.zeros(image.shape)
    
#     # loop over every image pixel
#     for i in range(0, image.shape[0]):
#         for j in range(0, image.shape[1]):     
#             if numpy.isnan(image[i,j]) == True:
#                 continue
#             flux = p[0] * numpy.exp(-((i-p[1])**2.0/(2.0*p[3]**2.0) + (j-p[2])**2.0/(2.0*p[3]**2.0))) + p[4]
# 
#             # average the pixel value
#             modelImage[i,j] = flux
    modelFlux = p[0] * numpy.exp(-((yCo-p[1])**2.0/(2.0*p[3]**2.0) + (xCo-p[2])**2.0/(2.0*p[3]**2.0))) 
    
    # create difference image
    chisq = ((image - modelFlux)**2.0).sum()
    print chisq, p
    # check source is on image 
    #print chisq, p
    return chisq
    
###################################################################################

# list files in folder
files = os.listdir(folder)

# select just fits files
fitsFiles = [fitsFile for fitsFile in files if fitsFile[-5:] == ".fits"]

# identify the file in each band
fitsNames = {}
for band in bands:
    if len([True for fitsFile in fitsFiles if fitsFile.count(band)]) != 1:
        raise Exception("Not a single appropiate FITS file present")
    bandFile = [fitsFile for fitsFile in fitsFiles if fitsFile.count(band) == 1][0]
    fitsNames[band] = bandFile

# loop over each band
for band in bands:
    # load in data
    fits = pyfits.open(pj(folder,fitsNames[band]))
    data = fits[1].data
    header = fits[1].header
    print numpy.shape(data)
    # calculate the pixel size of the image
    try:
        cdelt1 = abs(header["CDELT1"] * 3600.0)
        cdelt2 = header["CDELT2"] * 3600.0
    except:
        cdelt1 = abs(header["CD1_1"] / math.cos((math.atan2(header["CD2_1"],header["CD1_1"]))) * 3600.0)
        cdelt2 = header["CD2_2"] / math.cos((math.atan2(-header["CD1_2"],header["CD2_2"]))) * 3600.0
    
    # Make input X and Y arrays for every pixel on the map
    xpix = numpy.zeros((header["NAXIS1"]*header["NAXIS2"]),dtype=int)
    for i in range(0,header["NAXIS2"]):
        xpix[i*header["NAXIS1"]:(i+1)*header["NAXIS1"]] = numpy.arange(0,header["NAXIS1"],1)
    ypix = numpy.zeros((header["NAXIS1"]*header["NAXIS2"]),dtype=int)
    for i in range(1,header["NAXIS2"]):
        ypix[(i)*header["NAXIS1"]:(i+1)*header["NAXIS1"]] = i
    xpix = xpix.reshape(data.shape)
    ypix = ypix.reshape(data.shape)
    sel1=numpy.where(numpy.isnan(data)==False)
    dataMin=data[sel1].min()
    data=data-dataMin
    print dataMin,data[sel1].min()
    
    # find the centre of the beam by fitting gaussian
    guessSig = FWHM[band] / (2.35482 * cdelt1)
    selection = numpy.where((numpy.isnan(data) == False) & (data > cutLevel[band]))
    dataMax = data[selection].max()
    guessCentre = numpy.where((data == dataMax) & (numpy.isnan(data) == False))
    p0 = [dataMax, guessCentre[0][0], guessCentre[1][0], guessSig]
    plsq = fmin(gaussFitter,p0,args=(data[selection],xpix[selection],ypix[selection]), full_output=True, disp=0, maxfun=1.0e9, maxiter=1.0e9)
    data=data/dataMax
    #data=numpy.where(numpy.isnan(data) == False,data/dataMax,data)
    # centre of image
    centre = [int(numpy.round(plsq[0][1])), int(numpy.round(plsq[0][2]))]
    #centre = [plsq[0][1], plsq[0][2]]
    
    # crop image to size
    if int(numpy.round(sizeImage/cdelt1)) %2 == 0:
        newSize = int(numpy.round(sizeImage/cdelt1)) + 1
    else:
        newSize = int(numpy.round(sizeImage/cdelt1))
        
    # select region of image to copy
    boxlim = [centre[0]-(newSize-1.0)/2.0,centre[0]+(newSize-1.0)/2.0+1,centre[1]-(newSize-1.0)/2.0,centre[1]+(newSize-1.0)/2.0+1]
    if boxlim[0] < 0:
        boxlim[0] = 0
    if boxlim[1] > data.shape[0]:
        boxlim[1] = data.shape[0]
    if boxlim[2] < 0:
        boxlim[2] = 0
    if boxlim[3] > data.shape[1]:
        boxlim[3] = data.shape[1]
    cutImage = data[boxlim[0]:boxlim[1],boxlim[2]:boxlim[3]]
    
    # create a new fits file for cut image
    hdu = pyfits.PrimaryHDU(cutImage)
    hdulist = pyfits.HDUList([hdu])
    hdu.header.update('EQUINOX', 2000.0)
    hdu.header.update('CTYPE1', header["CTYPE1"])
    hdu.header.update('CTYPE2', header["CTYPE2"])
    hdu.header.update('CRPIX1', (newSize-1.0)/2.0+1.0)
    hdu.header.update('CRPIX2', (newSize-1.0)/2.0+1.0)
    hdu.header.update('CRVAL1', 0.0)
    hdu.header.update('CRVAL2', 0.0)
    hdu.header.update('CDELT1', -cdelt1/3600.0)
    hdu.header.update('CDELT2', cdelt2/3600.0)
    try:
        hdu.header.update('LONPOLE', header["LONPOLE"])
        hdu.header.update('LATPOLE', header["LATPOLE"])
    except:
        pass
    outfile = fitsNames[band][:-5] + "-cutTrim.fits"
    hdulist.writeto(pj(folder,outfile))
    