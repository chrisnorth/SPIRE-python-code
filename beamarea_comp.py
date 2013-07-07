# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:38:07 2012

@author: chris
"""

from numpy import zeros_like,zeros,arange,isfinite,sqrt,max
from scipy import where
import pyfits
import matplotlib.pyplot as plot
import sys

beamdir='../Inputs/spire_beams_varpix/'
beamdir0='../Inputs/spire_beams_measured/'
obsid='5000241a'

pix1=[1,1,1]
pix2=[6,10,14]
bands=['Psw','Pmw','Plw']
bands0=['psw','pmw','plw']
cpx0arr=[1000,1005,1013]
cpy0arr=[999,1004,1013]

files0=[]
files1=[]
files2=[]
for b in range(3):
    files0.append('%s%s_beam_1arcsec.fits'%(beamdir0,bands0[b]))
    files1.append('%s%s_Beams_multiPixelSize/%s_%sBeam_pxSize_%.1f.fits'%(beamdir,obsid,obsid,bands[b],pix1[b]))
    files2.append('%s%s_Beams_multiPixelSize/%s_%sBeam_pxSize_%.1f.fits'%(beamdir,obsid,obsid,bands[b],pix2[b]))

maps0=[]
maps1=[]
maps2=[]
rads0=[]
rads1=[]
rads2=[]

for b in range(3):
    print 'Band %d'%(b)
    print 'Map 0:'
    map0=pyfits.getdata(files0[b],0) #read in map
    map0=where(isfinite(map0),map0,0.) #replace NaNs and inf with 0
    hdr0=pyfits.getheader(files0[b],0)
    #get map size
    nx0=hdr0.get('NAXIS1')
    ny0=hdr0.get('NAXIS2')
    #get pixel size    
    dx0=hdr0.get('CDELT1')*3600. #in arcsec
    dy0=hdr0.get('CDELT2')*3600. #in arcsec
    #get centre pixels
    (cpx0,cpy0)=where(map0 == max(map0))
    print 'Centre:',cpx0,cpy0
    #make radius array
    rad0=zeros_like(map0)
    for x in range(nx0):
        rad0[x,:]=sqrt((x - cpx0)**2 + (arange(ny0) - cpy0)**2) * dx0

    maps0.append(map0)
    rads0.append(rad0)

    plot.figure(1) 
    plot.imshow(map0)
    plot.show()
    
    print 'Map 1:'
    map1=pyfits.getdata(files1[b],1)
    map1=where(isfinite(map1),map1,0.) #replace NaNs and inf with 0
    hdr1=pyfits.getheader(files1[b],1)
    #get map size
    nx1=hdr1.get('NAXIS1')
    ny1=hdr1.get('NAXIS2')
    #get pixel size    
    dx1=hdr1.get('CDELT1')*3600. #in arcsec
    dy1=hdr1.get('CDELT2')*3600. #in arcsec
    #get centre pixels
    (cpx1,cpy1)=where(map1 == max(map1))
    print 'Centre:',cpx1,cpy1
    #make radius array (in arcsec)
    rad1=zeros_like(map1)
    for x in range(nx1):
        rad1[x,:]=sqrt((x - cpx1)**2 + (arange(ny1) - cpy1)**2) * dx1
    
    maps1.append(map1)
    rads1.append(rad1)
    
    print 'Map 2:'
    map2=pyfits.getdata(files2[b],1)
    map2=where(isfinite(map2),map2,0.) #replace NaNs and inf with 0
    hdr2=pyfits.getheader(files2[b],1)
    #get map size
    nx2=hdr2.get('NAXIS1')
    ny2=hdr2.get('NAXIS2')
    #get pixel size    
    dx2=hdr2.get('CDELT1')*3600. #in arcsec
    dy2=hdr2.get('CDELT2')*3600. #in arcsec
    #get centre pixels
    (cpx2,cpy2)=where(map2 == max(map2))
    print 'Centre:',cpx2,cpy2
    #make radius array (in arcsec)
    rad2=zeros_like(map2)
    for x in range(nx2):
        rad2[x,:]=sqrt((x - cpx2)**2 + (arange(ny2) - cpy2)**2) * dx2
    
    maps2.append(map2)
    rads2.append(rad2)

maxr=1000.
radarr=arange(10,maxr+10.,10.)
nrad=len(radarr)
arear=zeros((nrad,3,3))

for r in range(nrad):
    for b in range(3):
        map0r=where(rads0[b] <= radarr[r],maps0[b],0.)
        arear[r,0,b]=sum(map0r)
        map1r=where(rads1[b] <= radarr[r],maps1[b],0.)
        arear[r,1,b]=sum(map1r)
        map2r=where(rads2[b] <= radarr[r],maps2[b],0.)
        arear[r,2,b]=sum(map2r)
        