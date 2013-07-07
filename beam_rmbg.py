# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 11:15:40 2012

@author: chris
"""

from scipy import where
from scipy.optimize import fmin
from numpy import isfinite,zeros,zeros_like,ones_like,arange,array
from numpy import sqrt,min,max,sum,average,arctan2,floor,pi,log10,mean
import pyfits
import sys
import matplotlib.pyplot as plot
import matplotlib.cm as cm
import string as str

from beam import getnewbeamprofile, beam_azsm, comb_beam


def maincode():

    #####################
    ###  Set Options  ###
    #####################
    
    plotallprofs=False
    plotsegs=False
    
    outdir='/data/Herschel/Calibration/Outputs'
    
    
    #beamdir='../Inputs/spire_beams_measured/'
    beamdir='/data/Herschel/Calibration/'
    bands=['psw','pmw','plw']
    
    filesIn=[]
    extnum=[]
    for b in range(3):
        filesIn.append('%sNepBeam_%s_1arcsec_combined.fits'%(beamdir,bands[b]))
        extnum.append(1)
        #files.append('%s%s_beam_1arcsec.fits'%(beamdir,bands[b]))
        #extnum.append(0)
    
    maps=[]
    hdrsIn=[]
    masks=[]
    #beamrads=[]
    #beamradints=[]
    #beamazs=[]
    beamxs=[]
    beamys=[]
    bgmasks=[]
    
    ###################################################
    ###  Read in beams
    ###################################################
    
    for b in range(3):
        print 'Reading band %d'%(b+1)    
        ##get fits data
        map0=pyfits.getdata(filesIn[b],extnum[b]) #read in map
        mapones=ones_like(map0) ## array of ones
        mapnonan=where(isfinite(map0),map0,0.) #replace NaNs and inf with 0
        maskval0=where(isfinite(map0),1,0) #replace NaNs and inf with 0
        hdr0=pyfits.getheader(filesIn[b],extnum[b])
        hdrsIn.append(hdr0)
        maps.append(map0)
        masks.append(maskval0)    
        
        #get map size
        nx0=hdr0.get('NAXIS1')
        ny0=hdr0.get('NAXIS2')
        #get pixel size    
        dx0=abs(hdr0.get('CDELT1'))*3600. #in arcsec
        dy0=abs(hdr0.get('CDELT2'))*3600. #in arcsec
        bpix=dx0
    
        #get centre pixels
        (cpx0,cpy0)=where(mapnonan == max(mapnonan))
    
        #make radius array
        #beamrad=zeros_like(map0)
        beamx=zeros_like(map0)
        beamy=zeros_like(map0)
        
        ##################################################
        ###   Calculate x position
        ##################################################
        
        for x in range(nx0):
            beamx[x,:]=(x - cpx0) * dx0
            beamy[x,:]=(arange(ny0) - cpy0) * dx0
            #beamrad[x,:]=sqrt((x - cpx0)**2 + (arange(ny0) - cpy0)**2) * dx0
        
        ####################################################
        ###   Compute radius and azimuth
        ####################################################
        
        print 'Calculating radius...'
        beamrad = sqrt(beamx**2 + beamy**2)
        ##print 'beamradint'    
        #beamradint=floor(beamrad)
        ##print 'beamaz'    
        #beamaz = arctan2(beamy,beamx) + pi
    
        beamxs.append(beamx)
        beamys.append(beamy)
        #beamrads.append(beamrad)
        #beamradints.append(beamradint)
        #beamazs.append(beamaz)
    
        ##refine to smaller box    
        #bwid = 2000
        #inbox = zeros_like(map0)
        #inbox[cpx0-bwid/2:cpx0+bwid/2,cpy0-bwid/2,cpy0+bwid/2] = 1
        
        #sys.exit()
        
        #map0=map0[cpx0-bwid/2:cpx0+bwid/2,cpy0-bwid/2,cpy0+bwid/2]
        #beamrad=beamrad[cpx0-bwid/2:cpx0+bwid/2,cpy0-bwid/2,cpy0+bwid/2]
        #beamrad=beamrad[cpx0-bwid/2:cpx0+bwid/2,cpy0-bwid/2,cpy0+bwid/2]
        
        ##make beam profile
        radarr=arange(1001.)
        nrad=len(radarr)    
    
        #print 'beamprof0,2,3'
        #beamprof0=beam_azsm(beamx,beamy,map0,radarr)
        #beamprof2=beam_azsm(beamx,beamy,map0,radarr,nsamp=720)
        #beamprof3=beam_azsm(beamx,beamy,map0,radarr,nsamp=1440)
        
        ##read in my beam
        #print 'beamprof1'    
        #(beam_scl,beam_fix)=getnewbeamprofile(b,radarr,bzlim=None)
        #beamprof1=comb_beam(beam_scl,beam_fix)
    
    
        bgmaskrad=240. ##radius of mask in arcsec
        bgmask=where(beamrad < bgmaskrad, 0, maskval0)
        bgmasks.append(bgmask)
    
    ##############################################
    ####  Bin maps to large pixels for BG removal
    ##############################################

    binxs=[]
    binys=[]
    binmaps=[]
    binmasks=[]
    binbgmasks=[]
        
    for b in range(3):
        ##Bin to 1arcmin
        bin1=60.
        nxb1=int(floor(nx0/bin1))
        nyb1=int(floor(ny0/bin1))
        print 'Binning band %d to %.1f" pixels (%d x %d)'%(b+1,bin1,nxb1,nyb1)
        
        binx=zeros((nxb1,nyb1))
        biny=zeros_like(binx)
        binmap=zeros_like(binx)
        binmask=zeros_like(binx)
        binbgmask=zeros_like(binx)
        #ylist=bin1*
        for x in range(nxb1):
            for y in range(nyb1):
                binx[x,y]=beamx[bin1*(x+0.5),bin1*(y+0.5)]
                biny[x,y]=beamy[bin1*(x+0.5),bin1*(y+0.5)]
                binmap[x,y]=mean(maps[b][bin1*x:bin1*(x+1),bin1*y:bin1*(y+1)])
                binmask[x,y]=mean(masks[b][bin1*x:bin1*(x+1),bin1*y:bin1*(y+1)])
                binbgmask[x,y]=mean(bgmasks[b][bin1*x:bin1*(x+1),bin1*y:bin1*(y+1)])

        binxs.append(binx)
        binys.append(biny)
        binmaps.append(binmap)
        binmasks.append(binmask)
        binbgmasks.append(binbgmask)
        
#    for b in range(3):
#        
#        extent=(min(binxs[b]),max(binxs[b]),min(binys[b]),max(binys[b]))
#        print extent
#        plot.figure()
#        plot.clf()
#        plot.subplot(2,2,1)
#        plot.imshow(maps[b],vmin=-1.e-3,vmax=1.e-3,extent=extent)
#        plot.colorbar()
#        plot.title('orig map (%s)'%bands[b])
#        
#        plot.subplot(2,2,2)
#        plot.imshow(masks[b],extent=extent)
#        plot.colorbar()
#        plot.title('orig mask')
#        
#        plot.subplot(2,2,3)
#        plot.imshow(binmaps[b],vmin=-1.e-3,vmax=1.e-3,extent=extent)
#        plot.colorbar()
#        plot.title('binned map')
#    
#        plot.subplot(2,2,4)
#        plot.imshow(binbgmasks[b],extent=extent)
#        plot.colorbar()
#        plot.title('binned mask')
#        
    #sys.exit()

    ##############################################
    ####  Calculate background fit
    ##############################################

    fixmaps=[]    
    
    for b in range(3):
        print ' '
        print 'Fitting band %d...'%(b+1)
        #x0init=8.521e-5
        #gradXinit=1.809e-6
        #gradYinit=3.502e-6
        x0init=0.
        gradXinit=0.
        gradYinit=0.
        
        paramInit=array([x0init,gradXinit,gradYinit])
    
        ##bgFit = fmin(fitbg,paramInit,args=(beamxs[b],beamys[b],maps[b],bgmasks[b]),disp=True,maxiter=25,xtol=0.0001)
        bgFit = fmin(fitbg,paramInit,args=(binxs[b],binys[b],binmaps[b],binbgmasks[b]),disp=True,maxiter=25,xtol=0.0001)
        print bgFit

        chisq,fitMap,chisqMap = fitbg_test(bgFit,binxs[b],binys[b],binmaps[b],binbgmasks[b])
        maskedMap=where(binbgmasks[b] == 1,binmaps[b],0.)    
        fixMap=maskedMap - fitMap
        reldifMap=(maskedMap - fixMap)/maskedMap
        reldifMap=where(binbgmasks[b] == 1,reldifMap,0.)
        meanDif=mean(reldifMap)/mean(binbgmasks[b])
        print 'Mean relative difference:',meanDif
        
        extentbin=(min(binxs[b]),max(binxs[b]),min(binys[b]),max(binys[b]))
        extent=(min(beamxs[b]),max(beamxs[b]),min(beamys[b]),max(beamys[b]))
        plot.figure()
        plot.subplot(2,2,1)
        plot.imshow(maskedMap,vmin=-1e-2,vmax=1e-2,extent=extentbin,interpolation='nearest',cmap=cm.hot)
        plot.title('Input Map (%s)'%bands[b])
        plot.colorbar()
        
        plot.subplot(2,2,2)
        plot.imshow(fitMap,extent=extentbin,interpolation='nearest',cmap=cm.hot)
        plot.title('Fit map')
        plot.colorbar()
        
        plot.subplot(2,2,3)
        plot.imshow(reldifMap,extent=extentbin,interpolation='nearest',cmap=cm.hot,vmin=-10.,vmax=10)
        plot.title('(Input-Fix)/Input')
        plot.colorbar()
        
        plot.subplot(2,2,4)
        plot.imshow(fixMap,extent=extentbin,vmin=-1e-2,vmax=1e-2,interpolation='nearest',cmap=cm.hot)
        plot.title('Fixed Map')
        plot.colorbar()    

        plot.savefig('../Outputs/beam_rmbg_%s.png'%(bands[b]),dpi=300)
        
        fixmapfull=maps[b] - (bgFit[0] + bgFit[1]*beamxs[b] + bgFit[2]*beamys[b])

        plot.figure()
        plot.subplot(1,2,1)
        plot.imshow(log10(maps[b]),extent=extent)
        plot.title('Input map')
        plot.colorbar()

        plot.subplot(1,2,2)
        plot.imshow(log10(fixmapfull),extent=extent)        
        plot.title('BG-removed map')
        plot.colorbar()

        plot.savefig('../Outputs/beam_rmbg_full_%s.png'%(bands[b]),dpi=300)

        ###################################
        ####   Write to FITS file
        ###################################

        print 'Writing to FITS file...'
        hdrOut=hdrsIn[b]
        fileOut='%sNepBeam_%s_1arcsec_combined_BGremoved.fits'%(beamdir,bands[b])
        hdrOut.add_comment('Background removed of form A + Bx + Cy')
        hdrOut.update('BG_A',bgFit[0],comment='background parameter A (zero level_')
        hdrOut.update('BG_B',bgFit[1],comment='background parameter B (x gradient)')
        hdrOut.update('BG_C',bgFit[2],comment='background parameter C (y gradient)')
        hdrOut.update('BG_SM',bin1,comment='Pixel bin size (") used for BG calc')
        hdrOut.update('BGSRCRAD',bgmaskrad,comment='Source radius excluded for BG calc')
        #hdrOut.update('NAXIS1',nx0)
        #hdrOut.update('NAXIS2',ny0)

        prihdu=pyfits.PrimaryHDU(header=hdrOut)
        hdunew=pyfits.ImageHDU(fixmapfull,header=hdrOut)
        hdulist=pyfits.HDUList([prihdu,hdunew])
        hdulist.writeto(fileOut,clobber=True)
        
    
    plot.show()
     

def fitbg(params,xMap,yMap,realMap,maskMap):

    bg0=params[0]
    gradX=params[1]
    gradY=params[2]
    
    #xMap=maps_in[0]
    #yMap=maps_in[1]
    #realMap=maps_in[2]
    #maskMap=maps_in[3]

    nMask=sum(maskMap)
    
    realMap=realMap + 1.
    fitMap=1. + bg0 + gradX*xMap + gradY*yMap
    
    chisqMap=where(maskMap == 1.,((realMap - fitMap)**2)/fitMap,0.)
    
    #chisq = sum(chisqMap)/(nMask - 3.)
    chisq = sum(chisqMap)

    #print params,chisq
    realMap=realMap-1.
    return(chisq)

def fitbg_test(params,xMap,yMap,realMap,maskMap):

    bg0=params[0]
    gradX=params[1]
    gradY=params[2]
    
    #xMap=maps_in[0]
    #yMap=maps_in[1]
    #realMap=maps_in[2]
    #maskMap=maps_in[3]

    nMask=sum(maskMap)
    
    fitMap=where(maskMap==1.,bg0 + gradX*xMap + gradY*yMap,0.)
    
    chisqMap=where(maskMap == 1.,((realMap - fitMap)**2)/fitMap,0.)
    
    chisq = sum(chisqMap)/(nMask - 3.)

    return(chisq,fitMap,chisqMap)
    
if __name__ == "__main__":
    maincode()