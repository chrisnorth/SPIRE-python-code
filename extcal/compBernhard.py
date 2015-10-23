##Compare Bernhard's maps

import pyfits
from numpy import arange,zeros_like,sqrt,floor,isfinite
from scipy import where
from beam import beam2prof,beam2prof2,beamarea_az,modbeam_area
import aplpy as ap
import matplotlib.pyplot as plot
import sys
###############################################################
### Read in Bernhard's Beam Maps
### Produce beam profiles
###############################################################

######################################
### Set filenames                  ###
######################################

fitsIn=['../Bernhard/0x5000241aL_PSW_pmcorr_1arcsec.fits', \
   '../Bernhard/0x5000241aL_PMW_pmcorr_1arcsec.fits', \
   '../Bernhard/0x5000241aL_PLW_pmcorr_1arcsec.fits']

fitsClean=['../Bernhard/0x5000241aL_PSW_pmcorr_1arcsec_cln_bgsub.fits', \
    '../Bernhard/0x5000241aL_PMW_pmcorr_1arcsec_cln_bgsub.fits', \
    '../Bernhard/0x5000241aL_PLW_pmcorr_1arcsec_cln_bgsub.fits']

fitsNorm=['../Bernhard/0x5000241aL_PSW_pmcorr_1arcsec_norm.fits', \
   '../Bernhard/0x5000241aL_PMW_pmcorr_1arcsec_norm.fits', \
   '../Bernhard/0x5000241aL_PLW_pmcorr_1arcsec_norm.fits']

fitsDiff=['../Bernhard/0x5000241aL_PSW_pmcorr_1arcsec_diff.fits', \
   '../Bernhard/0x5000241aL_PMW_pmcorr_1arcsec_diff.fits', \
   '../Bernhard/0x5000241aL_PLW_pmcorr_1arcsec_diff.fits']
   
bands=['PSW','PMW','PLW']

##################################################
### Set background level and central positions ###
##################################################

bgLev=[-8.e-6,-2.e-5,-2.e-5]
peakNorm=[150.49599008,96.70399888,58.608041759]
cPxs=[[965,1119],[966,1139],[977,1138]]

#############################################
######## Read in maps & calculate radii   ###
#############################################
#####    
#####mapsIn=[]
#####mapsClean=[]
#####mapsNorm=[]
#####mapsDiff=[]
#####masksIn=[]
#####masksClean=[]
#####masksDiff=[]
#####xArrs=[]
#####yArrs=[]
#####radArrs=[]
#####radIntArrs=[]
#####for b in range(3):
#####    print 'Reading maps for band %d'%(b+1)
#####    mapsIn.append(pyfits.getdata(fitsIn[b]))
#####    masksIn.append(where(isfinite(mapsIn[b]),1.,0.))
#####
#####    mapsClean.append(pyfits.getdata(fitsClean[b]))
#####    masksClean.append(where(isfinite(mapsClean[b]),1.,0.))
#####    #mapsClean[b]=where(isfinite(mapsClean[b]),mapsClean[b],0.)
#####    
#####    mapsNorm.append(mapsIn[b]/peakNorm[b] - bgLev[b])
#####    
#####    #mapsIn[b]=where(isfinite(mapsIn[b]),mapsIn[b],0.)
#####
#####    mapsDiff.append(mapsClean[b] - mapsNorm[b])
#####
#####    masksDiff.append(where(masksIn[b] == 1.,1.,0.))
#####    masksDiff[b]=where(masksClean[b] == 1.,masksDiff[b],0.)
#####    
#####    beamx=zeros_like(mapsIn[b])
#####    beamy=zeros_like(mapsIn[b])
#####
#####    hdr0=pyfits.getheader(fitsIn[b],0)
#####    hdr=pyfits.getheader(fitsIn[b],1)
#####    nx0=hdr.get('NAXIS1')
#####    ny0=hdr.get('NAXIS2')
#####    dx0=abs(hdr.get('CDELT1'))*3600. #in arcsec
#####
#####    print 'Writing to FITS files...'
#####    hdrOut=hdr0
#####    priHduNorm=pyfits.PrimaryHDU(header=hdr0)
#####    hduNorm=pyfits.ImageHDU(mapsNorm[b],header=hdrOut)
#####    hduListNorm=pyfits.HDUList([priHduNorm,hduNorm])
#####    hduListNorm.writeto(fitsNorm[b],clobber=True)
#####    
#####    priHduDiff=pyfits.PrimaryHDU(header=hdr0)
#####    hduDiff=pyfits.ImageHDU(mapsDiff[b],header=hdrOut)
#####    hduListDiff=pyfits.HDUList([priHduDiff,hduDiff])
#####    hduListDiff.writeto(fitsDiff[b],clobber=True)

######################
### Plot figures   ###
######################

for b in range(3):
##    print 'Plotting band %d: Input'%(b)
##    fig=ap.FITSFigure(fitsIn[b],hdu=1)
##    hdr=pyfits.getheader(fitsIn[b])
##    raC=hdr.get('CRVAL1')
##    decC=hdr.get('CRVAL2')
##    fig.show_grayscale(stretch='log',vmin=1e-5,invert=True)
##    fig.recenter(raC,decC,width=1800/3600.,height=1800/3600.)
##    #fig.hide_tick_labels()
##    #fig.hide_axis_labels()
##    fig.save(fitsIn[b][0:-5]+'.png')
##    fig.save(fitsIn[b][0:-5]+'.eps')
    
##    print 'Plotting band %d: Clean'%(b)
##    fig=ap.FITSFigure(fitsClean[b])
##    hdr=pyfits.getheader(fitsClean[b])
##    raC=hdr.get('CRVAL1')
##    decC=hdr.get('CRVAL2')
##    fig.show_grayscale(stretch='log',vmin=1e-5,vmax=1,invert=True)
##    fig.recenter(raC,decC,width=1800/3600.,height=1800/3600.)
##    #fig.hide_tick_labels()
##    #fig.hide_axis_labels()
##    fig.save(fitsClean[b][0:-5]+'.png')
##    #fig.save(fitsClean[b][0:-5]+'.eps')

##    print 'Plotting band %d: Norm'%(b)
##    fig=ap.FITSFigure(fitsNorm[b],hdu=1)
##    hdr=pyfits.getheader(fitsNorm[b])
##    raC=hdr.get('CRVAL1')
##    decC=hdr.get('CRVAL2')
##    fig.show_grayscale(stretch='log',vmin=1e-5,vmax=1,invert=True)
##    fig.recenter(raC,decC,width=1800/3600.,height=1800/3600.)
##    #fig.hide_tick_labels()
##    #fig.hide_axis_labels()
##    fig.save(fitsNorm[b][0:-5]+'.png')
##    #fig.save(fitsNorm[b][0:-5]+'.eps')
    
    print 'Plotting band %d: Diff'%(b)
    fig=ap.FITSFigure(fitsDiff[b],hdu=1)
    hdr=pyfits.getheader(fitsDiff[b])
    raC=hdr.get('CRVAL1')
    decC=hdr.get('CRVAL2')
    fig.show_grayscale(stretch='linear',vmin=-1e-4,vmax=1e-4,invert=True)
    fig.recenter(raC,decC,width=1800/3600.,height=1800/3600.)
    #fig.hide_tick_labels()
    #fig.hide_axis_labels()
    fig.save(fitsDiff[b][0:-5]+'.png')
    #fig.save(fitsDiff[b][0:-5]+'.eps')

plot.show()
