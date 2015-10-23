# -*- coding: utf-8 -*-
"""
Created on Tue Dec 11 17:41:10 2012

@author: chris
"""
from numpy import arange,zeros,array,min,max,zeros_like,mean
from scipy import where
import pyfits
from rsrf import getrsrf,getapf,bandedge
from beam import getbeammaps,getnewbeamprofile,comb_beam,beamarea_az,beamarea
from beam import get_effnu_az_newbeam,arcsec2sr,beamarea_az_nu_new,measbeam_new
import sys
#################################
##Set parameters
#################################

bzlim=None
bzval=0.
brad=600.
ind=0.85
apip=-1.

c=299792458. #speed of light
wlc=(250.e-6,350.e-6,500.e-6)
wlc=array(wlc)
nuc=c/wlc
band=arange(3)
anep=array((1.26,1.39,1.45))

############################
##  Set frequency limits  ##
############################
numin=300.e9 #300GHz
numax=1800.e9 #1790 GHz
dnu=0.3e9 #0.3GHz step
nnu=1+(numax - numin)/dnu
nul=range(int(nnu))
nuarr=array(nul)*dnu + numin

wlmin=c/numax
wlmax=c/numin

ilim=zeros((3,2))
nulim=zeros((3,2))
ilim_2=zeros((3,2))
nulim_2=zeros((3,2))
############################
##  Read in RSRF and APF  ##
############################
print 'Reading in RSRF and Aperture efficiency...'
rsrfarr=zeros((nnu,3))
apfarr=zeros((nnu,3))
for b in band:
    rsrfarr[:,b]=getrsrf('M',b,nuarr)
    apfarr[:,b]=getapf('R',b,nuarr)
    (ilim[b,:],nulim[b,:])=bandedge(nuarr,rsrfarr[:,b])
    (ilim_2[b,:],nulim_2[b,:])=bandedge(nuarr,rsrfarr[:,b],fact=2.)
    
#############################
##  Read in beam profiles  ##
#############################
print 'Reading in beam profiles...'
radarr=arange(brad)
nrad=radarr.size
beam_scl=zeros((nrad,3))
beam_fix=zeros((nrad,3))
beam_cmb=zeros((nrad,3))
areas_prof=zeros(3)
areas_prof_sr=zeros(3)
bsw=array([[70,76,78,87],[95,103,110,130],[135,145,169,180]])
fracarea_prof=zeros((nrad,3))
for b in band:
    ##read in beam profiles
    beam_scl[:,b],beam_fix[:,b]=getnewbeamprofile(b,radarr,bzlim=bzlim,bzval=bzval)
    beam_cmb[:,b]=comb_beam(beam_scl[:,b],beam_fix[:,b])
    areas_prof[b]=beamarea_az(radarr,beam_cmb[:,b],brad=brad)
    for r in arange(nrad):
        fracarea_prof[r,b]=beamarea_az(radarr,beam_cmb[:,b],brad=radarr[r])/areas_prof[b]
    areas_prof_sr[b]=arcsec2sr(areas_prof[b])

#########################
##  Read in beam maps  ##
#########################

print 'Reading in beam maps...'
#fitsfiles=['../CalProducts/0x5000241aL_PSW_pmcorr_1arcsec_cln_bgsub.fits', \
#      '../CalProducts/0x5000241aL_PMW_pmcorr_1arcsec_cln_bgsub.fits', \
#      '../CalProducts/0x5000241aL_PLW_pmcorr_1arcsec_cln_bgsub.fits']
#
#crpix=[[965,1118],\
#    [965,1138],\
#    [976,1139]]
#    
#areas_map=zeros(3)
#areas_map_sr=zeros(3)
#beam_map=[]
#xarr_map=[]
#yarr_map=[]
#for b in range(3):
#    #print '---\nReading in beam map for band %d'%(b+1)
#    file=fitsfiles[b]
#    beamin=pyfits.getdata(file,0)
#    beamin=where(beamin==beamin,beamin,0.)
#    hdr=pyfits.getheader(file,0)
#    nxin=hdr.get('NAXIS1')
#    nyin=hdr.get('NAXIS2')
#    dx=hdr.get('CDELT1')*3600. #deg->arcsec
#    dy=hdr.get('CDELT1')*3600. #deg->arcsec
#    cpx=crpix[b][0]
#    cpy=crpix[b][1]
#    xin=abs(dx)*(array(range(0,nxin)) - cpx)
#    yin=abs(dy)*(array(range(0,nyin)) - cpy)
#    xarrin=zeros((nyin,nxin))
#    for x in range(0,nxin):
#        xarrin[:,x]=float(x) - cpx
#    yarrin=zeros((nyin,nxin))
#    for y in range(0,nyin):
#        yarrin[y,:]=float(y) - cpy
#    #print '(%d,%d):%.4f'%(xarrin[where(beamin==max(beamin))],yarrin[where(beamin==max(beamin))],beamin[where(beamin==max(beamin))])
#    xz=where(xin==0)[0]
#    yz=where(yin==0)[0]
#    #print 'beam[xmax,ymax]:',beamin[yz,xz]
#    #print 'max(beam)',max(beamin)
#    
#    beam_map.append(beamin)
#    xarr_map.append(xarrin)
#    yarr_map.append(yarrin)
#
xarr_map,yarr_map,beam_map=getbeammaps()
areas_map=zeros(3)
areas_map_sr=zeros(3)
for b in band:
    areas_map[b]=beamarea(xarr_map[b],yarr_map[b],beam_map[b],brad=brad,pix=1.)
    areas_map_sr[b]=arcsec2sr(areas_map[b])

########################
##  Regrid beam maps  ##
########################

print 'Regridding beam maps...'
regrid=[6,10,14]
areas_map_pix=zeros(3)
areas_map_pix_sr=zeros(3)
#for b in band:
#    xr=720
#    yr=xr
#    beamin=beam_map[b]
#    xarrin=xarr_map[b]
#    yarrin=yarr_map[b]
#    cpx=crpix[b][0]
#    cpy=crpix[b][1]
#    ixlim=[cpx-xr-regrid[b]/2,cpx+xr+regrid[b]/2]
#    iylim=[cpy-yr-regrid[b]/2,cpy+yr+regrid[b]/2]
#    ixout=range(ixlim[0],ixlim[1],regrid[b])
#    iyout=range(iylim[0],iylim[1],regrid[b])    
#    nxout=len(ixout)
#    nyout=len(iyout)
#    beamout=zeros((nxout,nyout))
#    xarrout=zeros_like(beamout)
#    yarrout=zeros_like(beamout)
#    for x in range(len(ixout)):
#        for y in range(len(iyout)):
#            x0=ixout[x]
#            x1=x0+regrid[b]
#            y0=iyout[y]
#            y1=y0+regrid[b]
#    
#            beampatch=beamin[x0:x1,y0:y1]
#            beamout[x,y]=mean(beamin[x0:x1,y0:y1])
#            xarrout[x,y]=mean(xarrin[x0:x1,y0:y1])
#            yarrout[x,y]=mean(yarrin[x0:x1,y0:y1])
#            #if max(beampatch) > 0.9:
#                #print 'max(%d:%d,%d:%d)=%.4g ; (%.1f,%.1f)=%.4g)'%(x0,x1,y0,y1,max(beampatch),xarrout[x,y],yarrout[x,y],beamout[x,y])
#    #print 'max(beamin): %.5g'%max(beamin)
#    #print 'max(beamout): %.5g'%max(beamout)
#    beamout=beamout/max(beamout)

xarr_map_pix,yarr_map_pix,beam_map_pix=getbeammaps(regrid=[6,10,14])
for b in band:
    areas_map_pix[b]=beamarea(xarr_map_pix[b],yarr_map_pix[b],beam_map_pix[b],brad=brad,pix=regrid[b])
    #print 'Area (%d" pixels): %.2f (%.2f x 1" beam)'%(regrid[b],areas_map_pix[b],areas_map_pix[b]/areas_map[b])
    areas_map_pix_sr[b]=arcsec2sr(areas_map_pix[b])

#sys.exit()
##########################################
##  Calculate effective beam frequency  ##
##########################################

print 'Calculating Effective frequencies (radial beam)...'
aprec=0.0001
if ind==0:
    nueff=nuc
else:
    nueff=get_effnu_az_newbeam(radarr,beam_scl,beam_fix,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=False,ind=ind)
dnu_a=10
#area_bm=beamarea_az_nu_new(nuarr,radarr,beam_scl,beam_fix,nueff,ilim_2,dnu_a=dnu_a,verbose=False,ind=ind)

#####################################
##  Calculate pipeline beam areas  ##
#####################################

print 'Calculating pipeline areas...'
areas_pip=zeros(3)
areas_pip_sr=zeros(3)
pipfact=zeros(3)
areas_map_pip=zeros(3)
areas_map_pip_sr=zeros(3)
areas_map_pix_pip=zeros(3)
areas_map_pix_pip_sr=zeros(3)
pipe_conv=zeros(3)
for b in band:
    areas_pip[b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],apip,ind=ind)
    areas_pip_sr[b]=arcsec2sr(areas_pip)[b]

    pipfact[b]=areas_pip[b]/areas_prof[b]
    areas_map_pip[b]=areas_map[b] * pipfact[b]
    areas_map_pip_sr[b]=arcsec2sr(areas_map_pip[b])

    pipe_conv[b] = 1.e-6 / areas_map_pip_sr[b]

    areas_map_pix_pip[b]=areas_map_pix[b] * pipfact[b]
    areas_map_pix_pip_sr[b]=arcsec2sr(areas_map_pix_pip[b])

#fracarea_pip=zeros((nrad,3))
#for b in band:
#    for r in arange(nrad):
#        fracarea_pip[r,b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],apip,ind=ind,brad=radarr[r])/areas_pip[b]

######################################
##  Calculate aperture corrections  ##
######################################

print 'Calculating aperture efficiency...'

apRad=[22.,30.,42.]
anRad=[60.,90.]
fracAp=zeros(3)
fracAn=zeros(3)
fracAp_pip=zeros(3)
fracAn_pip=zeros(3)

for b in band:
    fracAp_pip[b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],apip,ind=ind,brad=apRad[b])/areas_pip[b]
    anIn_pip=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],apip,ind=ind,brad=anRad[0])
    anOut_pip=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],apip,ind=ind,brad=anRad[1])
    fracAn_pip[b]=(anOut_pip-anIn_pip)/areas_pip[b]

    fracAp[b]=beamarea(xarr_map[b],yarr_map[b],beam_map[b],brad=apRad[b],pix=1.)/areas_map[b]
    anIn=beamarea(xarr_map[b],yarr_map[b],beam_map[b],brad=anRad[0],pix=1.)
    anOut=beamarea(xarr_map[b],yarr_map[b],beam_map[b],brad=anRad[1],pix=1.)
    fracAn[b]=(anOut-anIn)/areas_map[b]

apCorr=1./(fracAp-fracAn)
apCorr_pip=1./(fracAp_pip-fracAn_pip)

#######################
##  Writing to file  ##
#######################

print 'Writing values to file...'

fileout='../CalProducts/Beams.dat'
fout=open(fileout,'w')
fout.write('Parameters\n----------\n')
fout.write('  Maximum beam radius: %.1f arcsec\n'%(brad))
fout.write('  Central wavelengths: (%.0f, %.0f, %.0f) um\n'%(wlc[0]*1.e6,wlc[1]*1.e6,wlc[2]*1.e6))
fout.write('  Central Frequencies: (%.0f, %.0f, %.0f) GHz\n'%(nuc[0]/1.e9,nuc[1]/1.e9,nuc[2]/1.e9))
fout.write('  Power law index of core beam profile with frequency: %.2f\n'%(ind))
fout.write('  Neptune spectral index: (%.2f, %.2f, %.2f)\n'%(anep[0],anep[1],anep[2]))
fout.write('  Pipeline source spectral index: %.1f\n'%(apip))

fout.write('\nBeam area as measured on Neptune\n----------\n')
fout.write('  Beam map integrated to %d": (%.0f, %.0f, %.0f) arcsec^2\n'%(int(brad),areas_map[0],areas_map[1],areas_map[2]))
fout.write('                             : (%.3f, %.3f, %.3f) x 10^-8 sr\n'%(areas_map_sr[0]*1.e8,areas_map_sr[1]*1.e8,areas_map_sr[2]*1.e8))
fout.write('Areas from model beam profiles\n----------\n')
fout.write('  Beam profile integrated to %d": (%.0f, %.0f, %.0f) arcsec^2\n'%(int(brad),areas_prof[0],areas_prof[1],areas_prof[2]))
fout.write('                                 : (%.3f, %.3f, %.3f) x 10^-8 sr\n'%(areas_prof_sr[0]*1.e8,areas_prof_sr[1]*1.e8,areas_prof_sr[2]*1.e8))

fout.write('\nNeptune Beam areas with (6,10,14)" pixels\n----------\n'%())
fout.write('  Beam area: (%.0f, %.0f, %.0f) arcsec^2\n'%(areas_map_pix[0],areas_map_pix[1],areas_map_pix[2]))
fout.write('           : (%.3f, %.3f, %.3f) x 10^-8 sr\n'%(areas_map_pix_sr[0]*1.e8,areas_map_pix_sr[1]*1.e8,areas_map_pix_sr[2]*1.e8))

fout.write('\nEffective monochromatic frequency of Neptune maps:\n----------\n')
fout.write('  Effective frequencies:  (%.1f, %.1f, %.1f) GHz\n' % (nueff[0]/1.e9,nueff[1]/1.e9,nueff[2]/1.e9))
fout.write('  Effective wavelengths:  (%.1f, %.1f, %.1f) um\n' % (c/nueff[0]*1.e6,c/nueff[1]*1.e6,c/nueff[2]*1.e6))

fout.write('\nPipeline beam areas (alpha=-1)\n----------\n')
fout.write('  Beam area with 1" pixels: (%.0f, %.0f, %.0f) arcsec^2\n'%(areas_map_pip[0],areas_map_pip[1],areas_map_pip[2]))
fout.write('                          : (%.3f, %.3f, %.3f) x 10^-8 sr\n'%(areas_map_pip_sr[0]*1.e8,areas_map_pip_sr[1]*1.e8,areas_map_pip_sr[2]*1.e8))
fout.write('  Beam area with (6,10,14)" pixels: (%.0f, %.0f, %.0f) arcsec^2\n'%(areas_map_pix_pip[0],areas_map_pix_pip[1],areas_map_pix_pip[2]))
fout.write('                                  : (%.3f, %.3f, %.3f) x 10^-8 sr\n'%(areas_map_pix_pip_sr[0]*1.e8,areas_map_pix_pip_sr[1]*1.e8,areas_map_pix_pip_sr[2]*1.e8))
fout.write('\nPoint source -> extended source pipeline conversion factor\n----------\n')
fout.write('  Conversion factor: (%.2f, %.2f, %.2f) MJy/sr per Jy/beam\n'%(pipe_conv[0],pipe_conv[1],pipe_conv[2]))

#fout.write('\nAperture photometry parameters\n---------\n')
#fout.write('Aperture radius: (%.1f, %.1f, %.1f) arcsec\n'%(apRad[0],apRad[1],apRad[2]))
#fout.write('Background inner-outer radius: %.1f-%.1f arcsec (all bands)\n'%(anRad[0],anRad[1]))
#fout.write('\nAperture correction for Neptune source (1" pixels)\n----------\n')
#fout.write('Fraction of beam in aperture: (%.5g, %.5g, %.5g)\n'%(fracAp[0],fracAp[1],fracAp[2]))
#fout.write('Fraction of beam in annulus: (%.5g, %.5g, %.5g)\n'%(fracAn[0],fracAn[1],fracAn[2]))
#fout.write('Aperture correction factor: (%.5g, %.5g, %.5g)\n'%(apCorr[0],apCorr[1],apCorr[2]))

#fout.write('Aperture photometry for alpha=-1 source\n---------\n')
#fout.write('Fraction of beam in aperture: (%.5g, %.5g, %.5g)\n'%(fracAp_pip[0],fracAp_pip[1],fracAp_pip[2]))
#fout.write('Fraction of beam in annulus: (%.5g, %.5g, %.5g)\n'%(fracAn_pip[0],fracAn_pip[1],fracAn_pip[2]))
#fout.write('Aperture correction factor: (%.5g, %.5g, %.5g)\n'%(apCorr_pip[0],apCorr_pip[1],apCorr_pip[2]))

fout.write('\nAperture photometry parameters\n---------\n')
fout.write('  Aperture radius: (%.0f, %.0f, %.0f) arcsec\n'%(apRad[0],apRad[1],apRad[2]))
fout.write('  Background inner-outer radius: %.0f-%.0f arcsec (all bands)\n'%(anRad[0],anRad[1]))
fout.write('  Aperture correction factor for Neptune source: (%.3f, %.3f, %.3f)\n'%(apCorr[0],apCorr[1],apCorr[2]))
fout.write('  Aperture correction factor for alpha=-1 source: (%.3f, %.3f, %.3f)\n'%(apCorr_pip[0],apCorr_pip[1],apCorr_pip[2]))

fout.close()