# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 08:30:25 2012

@author: chris
"""

##measure beam areas

import pyfits
import sys
from numpy import array,zeros,zeros_like,min,max,mean,floor,remainder,NaN
from scipy import where
from beam import beamarea,arcsec2sr

fitsfiles=['../CalProducts/0x5000241aL_PSW_pmcorr_1arcsec_cln_bgsub.fits', \
      '../CalProducts/0x5000241aL_PMW_pmcorr_1arcsec_cln_bgsub.fits', \
      '../CalProducts/0x5000241aL_PLW_pmcorr_1arcsec_cln_bgsub.fits']

#fitsfiles=['../Inputs/spire_beams_measured/psw_beam_1arcsec.fits', \
#      '../Inputs/spire_beams_measured/pmw_beam_1arcsec.fits', \
#      '../Inputs/spire_beams_measured/plw_beam_1arcsec.fits']

crpix=[[965,1118],\
    [965,1138],\
    [976,1139]]
    
regrid=[6,10,14]

fitsfiles_regrid=['../CalProducts/0x5000241aL_PSW_pmcorr_6arcsec_cln_bgsub.fits', \
            '../CalProducts/0x5000241aL_PMW_pmcorr_10arcsec_cln_bgsub.fits', \
            '../CalProducts/0x5000241aL_PLW_pmcorr_14arcsec_cln_bgsub.fits']

for b in range(3):
    print '---\nBand %d'%(b+1)
    file=fitsfiles[b]
    beamin=pyfits.getdata(file,0)
    beamin=where(beamin==beamin,beamin,0.)
    hdr=pyfits.getheader(file,0)
    nxin=hdr.get('NAXIS1')
    nyin=hdr.get('NAXIS2')
    dx=hdr.get('CDELT1')*3600. #deg->arcsec
    dy=hdr.get('CDELT2')*3600. #deg->arcsec
    cpx=crpix[b][0]
    cpy=crpix[b][1]
    xin=abs(dx)*(array(range(0,nxin)) - cpx)
    yin=abs(dy)*(array(range(0,nyin)) - cpy)
    xarrin=zeros((nyin,nxin))
    for x in range(0,nxin):
        xarrin[:,x]=float(x) - cpx
    yarrin=zeros((nyin,nxin))
    for y in range(0,nyin):
        yarrin[y,:]=float(y) - cpy
    print '(%d,%d):%.4f'%(xarrin[where(beamin==max(beamin))],yarrin[where(beamin==max(beamin))],beamin[where(beamin==max(beamin))])
    xz=where(xin==0)[0]
    yz=where(yin==0)[0]
    #print 'beam[xmax,ymax]:',beamin[yz,xz]
    #print 'max(beam)',max(beamin)
    areain=beamarea(xarrin,yarrin,beamin,brad=700.,pix=1.)
    print 'Area (1" pixels): %.2f [%.4g]'%(areain,arcsec2sr(areain))
    
    print '...regridding...'
    xr=720
    yr=xr
    ixlim=[cpx-xr-regrid[b]/2,cpx+xr+regrid[b]/2]
    iylim=[cpy-yr-regrid[b]/2,cpy+yr+regrid[b]/2]
    ixout=range(ixlim[0],ixlim[1],regrid[b])
    iyout=range(iylim[0],iylim[1],regrid[b])    
    nxout=len(ixout)
    nyout=len(iyout)
    beamout=zeros((nxout,nyout))
    xarrout=zeros_like(beamout)
    yarrout=zeros_like(beamout)
    bmax=0.
    for x in range(len(ixout)):
        for y in range(len(iyout)):
            x0=ixout[x]
            x1=x0+regrid[b]
            y0=iyout[y]
            y1=y0+regrid[b]

            beampatch=beamin[x0:x1,y0:y1]
            beamout[x,y]=mean(beamin[x0:x1,y0:y1])
            xarrout[x,y]=mean(xarrin[x0:x1,y0:y1])
            yarrout[x,y]=mean(yarrin[x0:x1,y0:y1])
            if max(beampatch) > 0.9:
                print 'max(%d:%d,%d:%d)=%.4g ; (%.1f,%.1f)=%.4g)'%(x0,x1,y0,y1,max(beampatch),xarrout[x,y],yarrout[x,y],beamout[x,y])
            if beamout[x,y]>bmax:
                bmax=beamout[x,y]
                xmax=x
                ymax=y
    print 'max(beamin): %.5g'%max(beamin)
    print 'max(beamout): %.5g'%max(beamout)
    beamout=beamout/max(beamout)
    areaout=beamarea(xarrout,yarrout,beamout,brad=700.,pix=regrid[b])
    print 'Area (%d" pixels): %.2f (%.2f x 1" beam)'%(regrid[b],areaout,areaout/areain)
    
    ###making full regrid
    nxr=int(floor(nxin/regrid[b]))
    nyr=int(floor(nyin/regrid[b]))
    xoff=remainder(cpx,regrid[b])
    yoff=remainder(cpy,regrid[b])
    print 'offset:',xoff,yoff
    beamregrid=zeros((nxr-1,nyr-1))
    brmax=0
    cpxr=floor((cpx-xoff)/regrid[b])
    cpyr=floor((cpy-yoff)/regrid[b])
    for x in range(nxr-1):
        for y in range(nyr-1):
            beamregrid[x,y]=mean(beamin[x*regrid[b]+xoff:(x+1)*regrid[b]+xoff,y*regrid[b]+yoff:(y+1)*regrid[b]+yoff])
            if beamregrid[x,y]==0:
                beamregrid[x,y]=NaN
            if beamregrid[x,y] >= brmax:
                brmax=beamregrid[x,y]
                #cpxr=x
                #cpyr=y
    print brmax,max(beamregrid),cpxr,cpyr
    beamregrid=beamregrid/brmax
    ##write regridded beam maps to file
    hdrout=hdr.copy()
    hdrout.update('CDELT1',dx*regrid[b]/3600.)
    hdrout.update('CDELT2',dy*regrid[b]/3600.)
    hdrout.update('CRPIX1',cpxr)
    hdrout.update('CRPIX2',cpyr)
    pyfits.writeto(fitsfiles_regrid[b],data=beamregrid,header=hdrout,clobber=True)
    #sys.exit()