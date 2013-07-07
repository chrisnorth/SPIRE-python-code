from __future__ import division
from numpy import array,arange,ceil,floor,min,max,zeros,log10,pi,sin,cos
import argparse
import matplotlib.pyplot as plot
from scipy.interpolate import RectBivariateSpline
from scipy import where

from beam import getbeam,rotbeam,beam_azsm,modbeam,get_effnu
from rsrf import getrsrf,bandedge

def maincode():
    
    parser=argparse.ArgumentParser(description='Herschel-SPIRE Beam comparison')
    parser.add_argument('--bwid',action='store',default=750.,type=float,dest='bwid',help='Full width of beam read from file (arcsec). Default=2000')
    
    args=parser.parse_args()
    bwid=args.bwid
    band = arange(3)
    
    c=299792458. #speed of light
    #set band centres
    wlc=(250.e-6,350.e-6,500.e-6)
    wlc=array(wlc)
    nuc=c/wlc

#################
##  Get beams  ##
#################

    bgrid=1. #beam grid size in arcsec
    npx=bwid/bgrid
    beams_m=zeros((npx,npx,3))
    beams_t=zeros((npx,npx,3))
    print 'Getting beams...'
    for b in band:
        beamx_m,beamy_m,beams_m[:,:,b]=getbeam('M',b,bgrid,bwid)
        beamx_t,beamy_t,beams_t[:,:,b]=getbeam('T',b,bgrid,bwid)
    
    angs=arange(0,25,5)    
    ncut=angs.size

######################
##  Make beam cuts  ##
######################
    
    print 'Making beam cuts...'
    beam_mxcut=zeros((npx,3,ncut))
    beam_mycut=zeros((npx,3,ncut))
    beam_txcut=zeros((npx,3))
    beam_tycut=zeros((npx,3))

    xcut0_x=beamx_m[:,0]
    xcut0_y=zeros(npx)
    ycut0_x=zeros(npx)
    ycut0_y=beamy_m[0,:]
    for b in band:        
        beam_mspl=RectBivariateSpline(beamx_m[:,0],beamy_m[0,:],beams_m[:,:,b])
        for c in arange(ncut):
            ang=angs[c]*pi/180.
            xcut_x=xcut0_x*cos(ang)-xcut0_y*sin(ang)
            xcut_y=xcut0_x*sin(ang)+xcut0_y*cos(ang)
            ycut_x=ycut0_x*cos(ang)-ycut0_y*sin(ang)
            ycut_y=ycut0_x*sin(ang)+ycut0_y*cos(ang)
            beam_mxcut[:,b,c]=beam_mspl.ev(xcut_x,xcut_y)
            beam_mycut[:,b,c]=beam_mspl.ev(ycut_x,ycut_y)
        beam_tspl=RectBivariateSpline(beamx_t[:,0],beamy_t[0,:],beams_t[:,:,b])
        beam_txcut[:,b]=beam_tspl.ev(xcut0_x,xcut0_y)
        beam_tycut[:,b]=beam_tspl.ev(ycut0_x,ycut0_y)

    plotbcuts(xcut0_x,ycut0_y,beam_mxcut,beam_mycut,beam_txcut,beam_tycut)

##############################
##  Make az-smoothed beams  ##
##############################

    print 'Making az-smoothed beams...'
    radarr=arange(0.,max(beamx_m[:,0]),1.)
    nrad=radarr.size
    beams_maz=zeros((nrad,3))
    beams_taz=zeros((nrad,3))
    for b in band:    
        beams_maz[:,b]=beam_azsm(beamx_m,beamy_m,beams_m[:,:,b],radarr)
        beams_taz[:,b]=beam_azsm(beamx_t,beamy_t,beams_t[:,:,b],radarr)
        
    #plot.figure(4)
    #plot.subplot

##########################
##  Define frequencies  ##
##########################

    numin=300.e9 #300GHz
    numax=1800.e9 #1790 GHz
    dnu=0.3e9 #0.3GHz step
    nnu=1+(numax - numin)/dnu
    nul=range(int(nnu))
    nuarr=array(nul)*dnu + numin
    
################
##  Get RSRF  ##
################

    print 'Getting RSRF...'
    rsrfarr=zeros((nnu,3))
    ilim=zeros((3,2))
    nulim=zeros((3,2))
    ilim_2=zeros((3,2))
    nulim_2=zeros((3,2))
    for b in band:
        rsrfarr[:,b]=getrsrf('M',b,nuarr)
        (ilim[b,:],nulim[b,:])=bandedge(nuarr,rsrfarr[:,b])
        (ilim_2[b,:],nulim_2[b,:])=bandedge(nuarr,rsrfarr[:,b],fact=2.)

##########################################
##  Calculate effective beam frequency  ##
##########################################
    #brad=350.
    print 'Calculating Effective frequencies...'
    aprec=0.0001
    brad=350.
    nueff_m=get_effnu(beamx_m,beamy_m,beams_m,rsrfarr,nuarr,nuc,brad,aprec=aprec)
    print 'Effective wavelengths (M): [%.2f, %.2f, %.2f] um' % (c/nueff_m[0]*1.e6,c/nueff_m[1]*1.e6,c/nueff_m[2]*1.e6)
    nueff_t=get_effnu(beamx_t,beamy_t,beams_t,rsrfarr,nuarr,nuc,brad,aprec=aprec)
    print 'Effective wavelengths (T): [%.2f, %.2f, %.2f] um' % (c/nueff_t[0]*1.e6,c/nueff_t[1]*1.e6,c/nueff_t[2]*1.e6)

##############################
##  Model theoretical beam  ##
##############################

    alpha_nep=array((1.26,1.39,1.45))
    modbeams=zeros((nrad,3))
    for b in band:
        nuarr_b=nuarr[ilim_2[b,0]:ilim_2[b,1]]
        rsrfarr_b=rsrfarr[ilim_2[b,0]:ilim_2[b,1],b]
        modbeams[:,b]=modbeam(radarr,beams_taz[:,b],nuarr,rsrfarr[:,b],nuc[b],alpha_nep[b])

    plotbazsm(radarr,beams_maz,beams_taz)

    #print 'Rotating beams...'
    #beamx_mrot,beamy_mrot=rotbeam(beamx_m,beamy_m,20)

    #plotbcomp(beamx_m,beamy_m,beams_m,beamx_mrot,beamy_mrot,beams_m)
    #plotbcomp(beamx_m,beamy_m,beams_m,beamx_t,beamy_t,beams_t)
    
    plot.show()

def plotbcomp(beamx_m,beamy_m,beams_m,beamx_t,beamy_t,beams_t):
    print 'Plotting beams...'
    plot.figure(1,figsize=(7,10))
    plot.clf()
    plot.hot()
    tit=('PSW','PMW','PLW')

    for b in arange(3):        
        pnum_m=2*b+1
        pnum_t=2*b+2
        
        #margins
        ym=1.1
        xm=1.0

        #color levels
        ncol=16
        ntick=8
        pl_m=array([-8.,0])
        pl_t=array([-8.,0])        
        pr_m=pl_m[1]-pl_m[0]
        pr_t=pl_t[1]-pl_t[0]
        cl_m=arange(pl_m[0],pl_m[1]+pr_m/ncol,pr_m/ncol)
        cl_t=arange(pl_t[0],pl_t[1]+pr_t/ncol,pr_t/ncol)
        tl_m=arange(pl_m[0],pl_m[1]+pr_m/ntick,pr_m/ntick)
        tl_t=arange(pl_t[0],pl_t[1]+pr_t/ntick,pr_t/ntick)
        
        plot.subplot(3,2,pnum_m)
        plot.axis('equal')
        plot.contourf(beamx_m,beamy_m,log10(beams_m[:,:,b]),levels=cl_m)
        plot.xlim((xm*min(beamx_m),xm*max(beamx_m)))
        if b==2:
            plot.xticks((floor(min(beamx_m)),0,ceil((max(beamx_m)))))
            plot.xlabel('x ["]')
        else:plot.xticks((floor(min(beamx_m)),0,ceil((max(beamx_m)))),('','',''))
        plot.ylim((ym*min(beamy_m),ym*max(beamy_m)))
        plot.yticks((floor(min(beamy_m)),0,ceil((max(beamy_m)))),rotation=90)
        plot.ylabel(tit[b])
        if b==0: plot.title('Measured')
        plot.colorbar(format='%.1f',ticks=tl_m)

        plot.subplot(3,2,pnum_t)
        plot.axis('equal')
        plot.contourf(beamx_t,beamy_t,log10(beams_t[:,:,b]),levels=cl_t)
        plot.xlim((xm*min(beamx_t),xm*max(beamx_t)))
        if b==2:
            plot.xticks((floor(min(beamx_t)),0,ceil((max(beamx_t)))))
            plot.xlabel('x ["]')
        else: plot.xticks((floor(min(beamx_t)),0,ceil((max(beamx_t)))),('','',''))
        plot.ylim((ym*min(beamy_t),ym*max(beamy_t)))
        plot.yticks((floor(min(beamy_t)),0,ceil((max(beamy_t)))),('','',''),rotation=90)
        if b==0: plot.title('Theoretical')
        plot.colorbar(format='%.1f',ticks=tl_t)

def plotbcuts(xcut0_x,ycut0_y,beam_mxcut,beam_mycut,beam_txcut,beam_tycut):
    
    plot.figure(2,figsize=(10,8))
    plot.clf()

    ncut=beam_mxcut[0,0,:].size
    clist=arange(ncut)/ncut
    #colc=[]
    #for n in ncut:colc.append((0,0,clist[n]))
    colc=('k','r','g','b','c','m','y')
    tit=('PSW','PMW','PLW')
    
    for b in arange(3):
        pnum_x=2*b+1
        pnum_y=2*b+2
        plot.subplot(3,2,pnum_x)
        plot.plot(xcut0_x,beam_txcut[:,b],color='r')
        for c in arange(ncut):
            plot.plot(xcut0_x,beam_mxcut[:,b,c],color=colc[c],linestyle='--')
        plot.yscale('log')
        plot.ylim(1e-5,1)
        plot.ylabel(tit[b])
        plot.xscale('log')
        plot.xlim(10,max(xcut0_x))
        if b==2: plot.xlabel('x ["]')
        
        plot.subplot(3,2,pnum_y)
        plot.plot(ycut0_y,beam_tycut[:,b],color='r')
        for c in arange(ncut):
            plot.plot(ycut0_y,beam_mycut[:,b,c],color=colc[c],linestyle='--')
        plot.yscale('log')
        plot.ylim(1e-5,1)
        plot.xscale('log')
        plot.xlim(10,max(ycut0_y))
        if b==2: plot.xlabel('y ["]')
        

def plotbazsm(radarr,beam_meas,beam_mod):
    plot.figure(3,figsize=(10,8))
    plot.clf()
    
    cols=('b','g','r')    
    tit=('PSW','PMW','PLW')
    for b in arange(3):
        plot.subplot(3,1,b+1)
        plot.plot(radarr,beam_meas[:,b],color=cols[b],linestyle='-',label='Measured beam')
        plot.plot(radarr,beam_mod[:,b],color=cols[b],linestyle='--',label='Modelled beam')
        plot.yscale('log')
        plot.ylabel(tit[b])
        plot.legend()


if __name__ == "__main__":
    maincode()
