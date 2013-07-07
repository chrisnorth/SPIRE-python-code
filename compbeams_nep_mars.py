# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 09:28:28 2013

@author: chris
"""
from numpy import zeros,arange,shape,array,zeros_like
from beam import getnewbeamprofile,comb_beam,getbeammaps
from csv import reader
import string
import sys
from scipy.interpolate import interp1d
from scipy.optimize import fmin
import matplotlib.pyplot as plot
#
# Compare beam profiles measured on Neptune and Mars
#

def maincode():
    
    #read in Neptune beam profile
    print 'Read in measured beam profiles...'    
    brad=700
    bzlim=None
    bzval=0.
    nep_rad=[arange(brad)]*3
    nrad=len(nep_rad[0])
    beam_scl=zeros((nrad,3))
    beam_fix=zeros((nrad,3))
    nep_prof=[]
    for b in range(3):
        beam_scl[:,b],beam_fix[:,b]=getnewbeamprofile(b,nep_rad[0],bzlim=bzlim,bzval=bzval)
        nep_prof.append(comb_beam(beam_scl[:,b],beam_fix[:,b]))
    
        #areas_prof[b]=beamarea_az(radarr,beam_cmb[:,b],brad=brad)
    #print 'Measured Beam areas (from profile): [%.2f,%.2f,%.2f]'%(areas_prof[0],areas_prof[1],areas_prof[2])

##    print 'Reading in Mars (1) beam maps...'
##    
##    xarr_mars1,yarr_mars1,beam_mars1=getbeammaps(src="Mars1",regrid=None)
##    xarr_mars1=array(xarr_mars1)
##    yarr_mars1=array(yarr_mars1)
##    beam_mars1=array(beam_mars1)
##    print shape(xarr_mars1)
##
##    print 'Making azimuthal beam profile...'
##    mars1_profs=[]
##    mars1_pmins=[]
##    mars1_pmaxs=[]
##    mars1_psds=[]
##    for b in range(3):
##        (mars1_prof,mars1_pmin,mars1_pmax,mars1_psd)=beam_azsm2(xarr_mars1[b,:,:],yarr_mars1[b,:,:],beam_mars1[b,:,:],nep_rad[b],retall=True,pess=None,nsamp=360)
##        mars1_profs.append(mars1_prof)
##        mars1_pmins.append(mars1_pmin)
##        mars1_pmaxs.append(mars1_pmax)
##        mars1_psds.append(mars1_psd)

    ###print 'Reading in Mars (2) beam maps...'
    ###xarr_mars2,yarr_mars2,beam_mars2=getbeammaps(src="Mars2",regrid=None)
    ###
    ###print 'Making azimuthal beam profile...'
    ###mars2_profs=[]
    ###mars2_pmins=[]
    ###mars2_pmaxs=[]
    ###mars2_psds=[]
    ###for b in range(3):
    ###    (mars2_prof,mars2_pmin,mars2_pmax,mars2_psd)=beam_azsm2(xarr_mars2[b,:,:],yarr_mars2[b,:,:],beam_mars2[b,:,:],radarr,retall=True,pess=None,nsamp=360)
    ###    mars2_profs.append(mars2_prof)
    ###    mars2_pmins.append(mars2_pmin)
    ###    mars2_pmaxs.append(mars2_pmax)
    ###    mars2_psds.append(mars2_psd)
    ###
    
    print 'Reading in real profiles (Mars 1)...'
    files=['../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_1arcsec_PSW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_1arcsec_PMW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_1arcsec_PLW.dat']
    mars1_rad=[]
    mars1_prof=[]
    for b in range(3):
        rows=reader(open(files[b],'r'))
        rad=[]
        prof=[]
        for row in rows:
            if string.find(row[0],'#')<0:
                rad.append(float(row[0]))
                prof.append(float(row[1]))
        mars1_rad.append(array(rad))
        mars1_prof.append(array(prof))
                
    print 'Reading in real profiles (Mars 2)...'
    files=['../Inputs/Mars_50004e4b/beamProfs_Mars_50004e4b_1arcsec_PSW.dat', \
        '../Inputs/Mars_50004e4b/beamProfs_Mars_50004e4b_1arcsec_PMW.dat', \
        '../Inputs/Mars_50004e4b/beamProfs_Mars_50004e4b_1arcsec_PLW.dat']
    mars2_rad=[]
    mars2_prof=[]
    for b in range(3):
        rows=reader(open(files[b],'r'))
        rad=[]
        prof=[]
        for row in rows:
            if string.find(row[0],'#')<0:
                rad.append(float(row[0]))
                prof.append(float(row[1]))
        mars2_rad.append(array(rad))
        mars2_prof.append(array(prof))

    print shape(mars1_rad),shape(mars1_rad[0]),shape(mars1_rad[1]),shape(mars1_rad[2])
    print shape(mars1_prof),shape(mars1_prof[0]),shape(mars1_rad[1]),shape(mars1_rad[2])
    #sys.exit()
    
    mfit1_prof=[]
    mfit2_prof=[]
    radint=[arange(40,200,1),arange(50,200,1),arange(75,200,1)]
    nri=zeros(3)
    fitpar1=[]
    fitpar2=[]
    print 'Matching Mars & Neptune'
    for b in range(3):
        nri[b]=len(radint[b])
        intnep=interp1d(nep_rad[b],nep_prof[b])
        intmars1=interp1d(mars1_rad[b],mars1_prof[b])
        intmars2=interp1d(mars2_rad[b],mars2_prof[b])
        pnep=intnep(radint[b])
        pmars1=intmars1(radint[b])
        pmars2=intmars2(radint[b])
        guessOff1=pnep[nri[b]/2] - pmars1[nri[b]/2]
        guessOff2=pnep[nri[b]/2] - pmars2[nri[b]/2]
        guessFac1=(pnep[0]-pnep[nri[b]-1])/(pmars1[0]-pmars1[nri[b]-1])
        guessFac2=(pnep[0]-pnep[nri[b]-1])/(pmars2[0]-pmars2[nri[b]-1])
        guess1=[guessOff1,guessFac1]
        guess2=[guessOff2,guessFac2]
        fitpar1.append(fmin(matchProfs,guess1,args=(pnep,pmars1),full_output=False,disp=0,maxfun=1.e9,maxiter=1.e9))
        fitpar2.append(fmin(matchProfs,guess2,args=(pnep,pmars2),full_output=False,disp=0,maxfun=1.e9,maxiter=1.e9))
        print fitpar1[b]
        print fitpar2[b]
        fit1=zeros_like(mars1_prof[b])
        fit2=zeros_like(mars2_prof[b])
        print shape(fit1),shape(fitpar1),shape(mars1_prof[b]),shape(fitpar1)
        #fit1[:]=fitpar1[0] + mars1_prof[b]*fitpar1[1]
        #fit2[:]=fitpar2[0] + mars2_prof[b]*fitpar2[1]
        mfit1_prof.append(fitpar1[b][0] + mars1_prof[b]*fitpar1[b][1])
        mfit2_prof.append(fitpar2[b][0] + mars2_prof[b]*fitpar2[b][1])
        #mfit1_prof.append(fit1)
        #mfit2_prof.append(fit2)
    print shape(nep_prof),shape(nep_prof[2])    
    print shape(mfit1_prof),shape(mfit1_prof[2])

    plot.figure(1)
    plot.clf()
    cols=['b','g','r']
    for b in range(3):
        #plot.axvline(radint[b][0],c='0.7',lw=1.5,ls='-')
        #plot.axvline(radint[b][0],c=cols[b],lw=1.5,ls='--')
        #plot.axvline(radint[b][nri[b]-1],c='0.7',lw=1.5,ls='-')
        #plot.axvline(radint[b][nri[b]-1],c=cols[b],lw=1.5,ls='--')
        plot.plot(nep_rad[b],nep_prof[b],c=cols[b],lw=2,ls='-')
        plot.plot(mars1_rad[b],mfit1_prof[b],c=cols[b],lw=2,ls='--')
        plot.plot(mars2_rad[b],mfit2_prof[b],c=cols[b],lw=2,ls=':')
#    plot.plot(nep_rad[1],nep_prof[1],'g-')
#    plot.plot(nep_rad[2],nep_prof[2],'r-')
#    plot.plot(mars1_rad[0],mfit1_prof[0],'b--')
#    plot.plot(mars1_rad[1],mfit1_prof[1],'g--')
#    plot.plot(mars1_rad[2],mfit1_prof[2],'r--')
#    plot.plot(mars2_rad[0],mfit2_prof[0],'b:')
#    plot.plot(mars2_rad[1],mfit2_prof[1],'g:')
#    plot.plot(mars2_rad[2],mfit2_prof[2],'r:')
    plot.yscale('log')
    plot.ylim(1e-6,1)
    plot.xlim(0,300)
    plot.title('Neptune vs Mars')
    plot.xlabel('Angular distance [arcsec]')
    plot.ylabel('Beam response')
    plot.savefig('../Beam_plots/Profiles_Mars_Nep.png')
    plot.draw()

    plot.figure(2)
    plot.clf()
    plot.title('Neptune vs Mars')
    for b in range(3):
        plot.subplot(3,1,b+1)
        plot.axvline(radint[b][0],c='0.7',lw=1.5,ls='-')
        plot.axvline(radint[b][0],c=cols[b],lw=1.5,ls='--')
        plot.axvline(radint[b][nri[b]-1],c='0.7',lw=1.5,ls='-')
        plot.axvline(radint[b][nri[b]-1],c=cols[b],lw=1.5,ls='--')
        plot.plot(nep_rad[b],nep_prof[b],c=cols[b],lw=2,ls='-')
        plot.plot(mars1_rad[b],mfit1_prof[b],c=cols[b],lw=2,ls='--')
        plot.plot(mars2_rad[b],mfit2_prof[b],c=cols[b],lw=2,ls=':')
        plot.yscale('log')
        plot.ylim(1e-8,1)
        plot.xlim(0,500)
        if b==1:plot.ylabel('Beam response')
        if b==2:plot.xlabel('Angular distance [arcsec]')
        str=r'Mars = $%.2f\times$Nep + %.3g'%(fitpar1[b][1],fitpar1[b][0])
        plot.annotate(str,(0.99,0.85),ha='right',xycoords='axes fraction')
    plot.savefig('../Beam_plots/Profiles_Mars_Nep_bands.png')
    plot.draw()        

    plot.show()


# define function
def matchProfs(p, *args):
    # function to fit a gaussian model to data
    
    # the model is f = p[0] * e^(-((x-p[1])^2/(2*p[3]^2) + (y-p[2])^2/(2*p[3]^2)) + p[4]
    
    # extract arguments
    prof0= args[0]
    prof1= args[1]
        
    # create blank map
    modelProfile = zeros(prof1.shape)
    
    modelProfile = p[0] + prof1*p[1]
    
    # create difference image
    chisq = ((prof0 - modelProfile)**2.0).sum()
    print chisq, p
    # check source is on image 
    #print chisq, p
    return chisq

def beam_azsm2(beamx,beamy,beamin,radarr,retall=None,pess=None,nsamp=360):

    ## Azimuthally-smooth the beam (alternative method)
    
    from numpy import zeros,arange,pi,cos,sin,average,min,max,std,reshape,array
    import scipy.interpolate as interpolate
    
    nrad=radarr.size
    npt=beamx.size
    x1=reshape(beamx,npt)
    y1=reshape(beamy,npt)
    beamin1=reshape(beamin,npt)
    
    pos1=zeros((npt,2))
    pos1[:,0]=x1
    pos1[:,1]=y1

    print 'interpolating map...'    
    #interp=interpolate.interp2d(x1,y1,beamin1)
    interp=interpolate.LinearNDInterpolator(pos1,beamin1)
    #interp=interpolate.RectBivariateSpline(beamx[:,0],beamy[0,:],beamin)
    
    print 'making profile'
    beam_sm=zeros(nrad)
    beam_max=zeros(nrad)
    beam_min=zeros(nrad)
    beam_sd=zeros(nrad)
    for r in arange(nrad):
        thetalist=arange(nsamp)*pi/180. * (360/nsamp) # full circle in radians
        xlist = radarr[r] * cos(thetalist)
        ylist = -radarr[r] * sin(thetalist)
        beam_r=interp.ev(xlist,ylist)
        beam_sm[r]=average(beam_r)
        beam_min[r]=min(beam_r)
        beam_max[r]=max(beam_r)
        beam_sd[r]=std(beam_r)
    
    if retall == True:
        return(beam_sm,beam_min,beam_max,beam_sd)
    else:
        return(beam_sm)

def beam_azsm(beamx,beamy,beamin,radarr,retall=None,pess=None,nsamp=360):

    ## Azimuthally-smooth a 3d beam map to produce a radial beam profile
    
    from numpy import zeros,arange,pi,cos,sin,average,min,max,std
    import scipy.interpolate as interpolate
    from scipy import where
    
    nrad=radarr.size
    interp=interpolate.RectBivariateSpline(beamx[:,0],beamy[0,:],beamin)
    
    beam_sm=zeros(nrad)
    beam_max=zeros(nrad)
    beam_min=zeros(nrad)
    beam_sd=zeros(nrad)
    for r in arange(nrad):
        thetalist=arange(nsamp)*(pi/180.) * (360./nsamp)# full circle in radians
        xlist = radarr[r] * cos(thetalist)
        ylist = -radarr[r] * sin(thetalist)
        beam_r=interp.ev(xlist,ylist)
        beam_sm[r]=average(beam_r)
        beam_min[r]=min(beam_r)
        beam_max[r]=max(beam_r)
        beam_sd[r]=std(beam_r)

    if pess:
        inmap=where(beamin > 0)
        maxx=max(abs(beamx[inmap]))
        maxy=max(abs(beamy[inmap]))
        maxrad=min((maxx,maxy))
        print 'Limiting beam profile to %d arcsec'%(maxrad)
        beam_sm=where(radarr > maxrad,0.,beam_sm)
    
    if retall == True:
        return(beam_sm,beam_min,beam_max,beam_sd)
    else:
        return(beam_sm)
        
if __name__=="__main__":
    maincode()