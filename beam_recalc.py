# -*- coding: utf-8 -*-
"""
Created on Mon Apr  2 13:58:25 2012

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
    
    ## Set map type (orig | bgrem)
    maptype='bgrem'
    print 'Reading map of type "%s"'%(maptype)
    plotallprofs=True
    plotsegs=False
    
    outdir='/data/Herschel/Calibration/Outputs'
    
    
    #beamdir='../Inputs/spire_beams_measured/'
    beamdir='/data/Herschel/Calibration/'
    bands=['psw','pmw','plw']
    
    files=[]
    extnum=[]
    for b in range(3):
        if maptype=='orig':
            files.append('%sNepBeam_%s_1arcsec_combined.fits'%(beamdir,bands[b]))
        elif maptype=='bgrem':
            files.append('%sNepBeam_%s_1arcsec_combined_BGremoved.fits'%(beamdir,bands[b]))
        else:
            print 'Unknown maptype. Must be "orig" or "bgrem"'
            sys.exit()
        extnum.append(1)
        #files.append('%s%s_beam_1arcsec.fits'%(beamdir,bands[b]))
        #extnum.append(0)
    
    maps=[]
    hdrsIn=[]
    masks=[]
    beamrads=[]
    beamradints=[]
    beamazs=[]
    beamxs=[]
    beamys=[]
    bgmasks=[]
    
    ###################################################
    ###  Read in beams
    ###################################################
    
    for b in range(3):
        print 'Reading band %d'%(b+1)    
        ##get fits data
        map0=pyfits.getdata(files[b],extnum[b]) #read in map
        mapones=ones_like(map0) ## array of ones
        mapnonan=where(isfinite(map0),map0,0.) #replace NaNs and inf with 0
        maskval0=where(isfinite(map0),1,0) #replace NaNs and inf with 0
        hdr0=pyfits.getheader(files[b],extnum[b])
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
        
        print 'Calculating radius and azimuth'
        beamrad = sqrt(beamx**2 + beamy**2)
        #print 'beamradint'    
        beamradint=floor(beamrad)
        #print 'beamaz'    
        beamaz = arctan2(beamy,beamx) + pi
    
        beamxs.append(beamx)
        beamys.append(beamy)
        beamrads.append(beamrad)
        beamradints.append(beamradint)
        beamazs.append(beamaz)
    
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

    
    #########################################################
    ####  Plot beam profiles
    #########################################################
    
    
    if plotallprofs:
        beamprofs=[]
        weightss=[]
    
        ##set up figure
        plot.figure(1,figsize=(6,10))
        cols=['b','g','r']
    
        for b in range(3):
            print 'beamprof4'
            beamprof4=zeros_like(radarr)
            beamprof4c=zeros_like(radarr)
            weights=zeros_like(radarr)
            mask=zeros_like(map0)
            for r in range(nrad):
                beamprof4[r]=average(mapnonan[where(beamradint == radarr[r])])
                ##compute weight to compensate for rectangular map and nans        
                weights[r]=sum(mapones[where(beamradint == radarr[r])]) / sum(maskval0[where(beamradint == radarr[r])])
                print r,radarr[r],beamprof4[r],weights[r]
                beamprof4c[r]=beamprof4[r]*weights[r]
        
            beamprofs.append(beamprof4)
            weightss.append(weights)
    
            ##plot to chart
            plot.subplot(3,1,b+1)
            plot.axhline(0.,color='k',linestyle=':')
            plot.ylim(-6e-2,6e-2)
            plot.plot(radarr,beamprof4,color=cols[b])
            plot.ylabel('Beam profile [%s]'%(str.upper(bands[b])))
            plot.ticklabel_format(axis='y',style='sci',scilimits=(12,12))
    
            if b == 2:
                plot.xlabel('Radius [arcsec]')
    
        plot.savefig(outdir+'beamprof_zoom.png',dpi=300,transparent=True)
        plot.savefig(outdir+'beamprof_zoom.eps',transparent=True)

        if maptype=='orig':
            fileout='%s/Nepbeam_Orig_profile.dat'%(outdir)
        elif maptype=='bgrem':
            fileout='%s/Nepbeam_BGremoved_profile.dat'%(outdir)
        fprof=open(fileout,'w')
        fprof.write('# radius ["], PSW, PMW, PLW\n')
        for r in range(nrad):
            fprof.write('%.1f , %.7g, %.7g, %.7g\n'%(radarr[r],beamprofs[0][r],beamprofs[1][r],beamprofs[2][r]))
        fprof.close()
    
    ###############################################################
    ### Plot segments
    ###############################################################
    
    if plotsegs:
        print 'Plotting segments...'
        ##set segments
        nseg=3
        segmid=(arange(nseg)+0.5)*2.*pi/nseg
        segwid=2.*pi/nseg
        segrad=radarr
        segprof=zeros([nrad,nseg])
        segprofs=[]
    
        plot.figure(figsize=(10,6))
        
        for b in range(1):
            print 'band %d'%(b+1)
            for s in range(nseg):
                azrange=[2.*pi*s/nseg,2.*pi*(s+1)/nseg]
                print 'segment %d: %.1f-%.1f (%.1f +/- %.1f)'%(s+1,azrange[0]*180/pi,azrange[1]*180/pi,segmid[s]*180/pi,segwid/2.*180/pi)
                segmap=where(abs(beamazs[b] - segmid[b]) <= segwid/2.,maps[b],0.)
                segmap=where(isfinite(segmap),segmap,0.)
                segones=where(abs(beamazs[b] - segmid[b]) <= segwid/2.,1.,0.)
                segval=where(abs(beamazs[b] - segmid[b]) <= segwid/2.,masks[b],0.)
                
                for r in range(nrad):
                    segprof[r,s]=average(segmap[where(beamradint == radarr[r])])
    
                if b == 0:
                    plot.subplot(1,nseg,s+1)
                    plot.imshow(log10(sqrt(segmap**2)))
                    plot.title('Segment %d'%(s+1))
                    
            segprofs.append(segprof)
    
        plot.savefig(outdir+'beamprof_segs_maps.png',dpi=300,transparent=True)
    
        plot.figure(figsize=(6,10))
        for b in range(1):
            for s in range(nseg):
                plot.plot(radarr,segprofs[b][:,s],label='Segment %d'%(s+1))
            plot.ylim(-1e-2,1e-2)
            plot.legend(loc='lower left')
            
        
        plot.savefig(outdir+'beamprof_segs_zoom.png',dpi=300,transparent=True)



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