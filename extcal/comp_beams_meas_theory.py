from __future__ import division
from numpy import array,arange,ceil,floor,min,max,zeros,log10,pi,sin,cos,shape
import argparse
import matplotlib.pyplot as plot
from scipy.interpolate import RectBivariateSpline
from scipy import where

from beam import getbeammaps,getbeammaps_theory,beam_azsm,modbeam_new,modbeam_area,measbeam_new
from beam import getnewbeamprofile,comb_beam,beamarea_az,beamarea,get_effnu_az_newbeam
from rsrf import getrsrf,bandedge


#
# Compare measured and theoretical beam profiles
#

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
    print 'Read in measured beam profiles...'    
    brad=350
    bzlim=None
    bzval=0.
    radarr=arange(brad)
    nrad=radarr.size
    beam_scl=zeros((nrad,3))
    beam_fix=zeros((nrad,3))
    beam_cmb=zeros((nrad,3))
    areas_prof=zeros(3)
    for b in band:
        beam_scl[:,b],beam_fix[:,b]=getnewbeamprofile(b,radarr,bzlim=bzlim,bzval=bzval)
        beam_cmb[:,b]=comb_beam(beam_scl[:,b],beam_fix[:,b])
        areas_prof[b]=beamarea_az(radarr,beam_cmb[:,b],brad=brad)
    print 'Measured Beam areas (from profile): [%.2f,%.2f,%.2f]'%(areas_prof[0],areas_prof[1],areas_prof[2])

    print 'Read in theoretical beam profiles...'    
    beam_scl_t=zeros((nrad,3))
    beam_fix_t=zeros((nrad,3))
    beam_cmb_t=zeros((nrad,3))
    areas_prof_t=zeros(3)
    for b in band:
        beam_scl_t[:,b],beam_fix_t[:,b]=getnewbeamprofile(b,radarr,bzlim=bzlim,bzval=bzval,type='T')
        beam_cmb_t[:,b]=comb_beam(beam_scl_t[:,b],beam_fix_t[:,b])
        areas_prof_t[b]=beamarea_az(radarr,beam_cmb_t[:,b],brad=brad)
    print 'Theoretical Beam areas (from profile): [%.2f,%.2f,%.2f]'%(areas_prof_t[0],areas_prof_t[1],areas_prof_t[2])

    print 'Reading in beam maps...'
    #xarr_map=[]
    #yarr_map=[]
    #beam_map=[]
    #bwid=brad*1.1
    #bpix=1.
    areas_map_t=zeros(3)
    xarr_map_t,yarr_map_t,beam_map_t=getbeammaps(src='Theory',regrid=None)
    pix=[6.,10.,14.]
    for b in band:
        #xarrin,yarrin,beamin=getbeam('M',b,bpix,bwid)
        #xarr_map.append(xarrin)
        #yarr_map.append(yarrin)
        #beam_map.append(beamin)
        areas_map_t[b]=beamarea(xarr_map_t[b],yarr_map_t[b],beam_map_t[b],brad=brad)*pix[b]**2
    print 'Beam Areas (from map) [sq. arcsec]: %.2f, %.2f, %.2f'%(areas_map_t[0], areas_map_t[1], areas_map_t[2])
    prof2map_t=areas_map_t/areas_prof_t
    print 'Prof -> Map factor: %.4f, %.4f, %.4f'%(prof2map_t[0],prof2map_t[1],prof2map_t[2])

#############################
#####  Define frequencies  ##
#############################
###
###    numin=300.e9 #300GHz
###    numax=1800.e9 #1790 GHz
###    dnu=0.3e9 #0.3GHz step
###    nnu=1+(numax - numin)/dnu
###    nul=range(int(nnu))
###    nuarr=array(nul)*dnu + numin
###    
###################
#####  Get RSRF  ##
###################
###
###    print 'Getting RSRF...'
###    rsrfarr=zeros((nnu,3))
###    ilim=zeros((3,2))
###    nulim=zeros((3,2))
###    ilim_2=zeros((3,2))
###    nulim_2=zeros((3,2))
###    for b in band:
###        rsrfarr[:,b]=getrsrf('M',b,nuarr)
###        (ilim[b,:],nulim[b,:])=bandedge(nuarr,rsrfarr[:,b])
###        (ilim_2[b,:],nulim_2[b,:])=bandedge(nuarr,rsrfarr[:,b],fact=2.)
###
#############################################
#####  Calculate effective beam frequency  ##
#############################################
###    #brad=350.
###    print 'Calculating Effective frequencies (radial beam)...'
###    alpha_nep=array((1.29,1.42,1.47))
###    aprec=0.001
###    ind=0.85
###    if ind==0:
###        nueff=nuc
###    else:
###        nueff=get_effnu_az_newbeam(radarr,beam_scl,beam_fix,rsrfarr,nuarr,nuc,brad,alpha_nep,aprec=aprec,verbose=False,ind=ind)
###    print 'Effective wavelengths:  [%.2f, %.2f, %.2f] um' % (c/nueff[0]*1.e6,c/nueff[1]*1.e6,c/nueff[2]*1.e6)
###    
###    print 'Calculating effective beam area (alpha=-1)...'
###    alphapip=-1    
###    area_pip=zeros(3)
###    for b in band:
###        print 'Band %d...'%(b)
###        area_pip[b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alphapip,ind=ind)
###
###    print 'Comparing modelled and measured areas...'
###    beam_mod=zeros((nrad,3))
###    beam_comp=zeros((nrad,3))
###    #areas_mod=zeros((nrad_a,3))
###    area_rarr2=arange(1,101,1)
###    nrad_a2=area_rarr2.size
###    areas_nep=zeros((nrad_a2,3))
###    areas_mod=zeros((nrad_a2,3))
###    areas_rel=zeros((nrad_a2,3))
###    indmod=0.85
###    for b in band:
###        print 'Band %d...' % b
###        beam_mod[:,b]=modbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],0,ind=indmod)
###        beam_comp[:,b]=(beam_mod[:,b]-beam_cmb[:,b])/areas_prof[b]
###        areas_nep[:,b]=modbeam_area(radarr,beam_cmb[:,b],area_rarr2)
###        areas_mod[:,b]=modbeam_area(radarr,beam_mod[:,b],area_rarr2)
###        ##correct for profile->map effect
###        areas_nep[:,b]=areas_nep[:,b]*prof2map[b]
###        areas_mod[:,b]=areas_mod[:,b]*prof2map[b]
###        areas_rel[:,b]=(areas_mod[:,b]-areas_nep[:,b])/areas_nep[:,b]
###        #areas_rel[:,b]=(modbeam_area[:,b]-areas_r[:,b])/areas_r[:,b]
###        #print 'areas (Nep,Mod,Rel):',areas_nep[:,b],areas_mod[:,b],areas_rel[:,b]
###        print 'min/max',b,min(areas_rel[:,b]),max(areas_rel[:,b])

    plot.figure(1)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    radarrp=where(radarr>0,radarr,1e-1)
    pmeas,=plot.plot(radarrp,beam_cmb[:,0],c='k',ls='-',label='Measured')
    ptheo,=plot.plot(radarrp,beam_cmb_t[:,0],c='k',ls='--',label='Theoretical')
    #pmod,=plot.plot(radarrp,beam_mod[:,0],c='k',ls=':',label=r'Modelled ($\gamma=1$)')
    for b in range(3):
        plot.plot(radarrp,beam_cmb[:,b],c=cols[b],ls='-',label=tit[b])
        plot.plot(radarrp,beam_cmb_t[:,b],c=cols[b],ls='--')
        #plot.plot(radarrp,beam_mod[:,b],c=cols[b],ls=':')
    plot.title(r'Measured vs Theoretical')
    plot.ylabel('Beam response')
    plot.ylabel('Angular distance [arcsec]')
    plot.xscale('log')
    plot.yscale('log')
    plot.legend(loc='lower left')
    plot.gca().lines.remove(pmeas)
    plot.gca().lines.remove(ptheo)
    #plot.gca().lines.remove(pmod)
    plot.savefig('../Plots/Beam_meas-vs-theory.png')
    plot.draw()

    plot.figure(2)
    plot.clf()
    for b in range(3):
        plot.plot(radarrp,beam_cmb[:,b]/beam_cmb_t[:,b],c=cols[b],ls='-',label=tit[b])
    plot.xscale('log')
    plot.yscale('log')
    plot.title(r'Measured $/$ Theoretical')
    plot.ylabel('Ratio')
    plot.xlabel('Angular distance [arcsec]')
    plot.legend(loc='lower left')
    plot.savefig('../Plots/Beam_meas-div-theory.png')
    plot.draw()

    plot.figure(3)
    plot.clf()
    pos,=plot.plot(radarrp,1.-beam_cmb_t[:,0]/beam_cmb[:,0],c='k',ls='-',label='positive')
    neg,=plot.plot(radarrp,beam_cmb_t[:,0]/beam_cmb[:,0]-1.,c='k',ls='--',label='negative')
    for b in range(3):
        plot.plot(radarrp,1.-beam_cmb_t[:,b]/beam_cmb[:,b],c=cols[b],ls='-')
        plot.plot(radarrp,beam_cmb_t[:,b]/beam_cmb[:,b]-1.,c=cols[b],ls='--')
    plot.xscale('log')
    plot.yscale('log')
    plot.title(r'(Meas $-$ Theory)/Meas')
    plot.ylabel('Difference')
    plot.xlabel('Angular distance [arcsec]')
    plot.legend(loc='lower left')
    plot.gca().lines.remove(pos)
    plot.gca().lines.remove(neg)
    plot.savefig('../Plots/Beam_meas-sub-theory.png')
    plot.draw()
    
    plot.figure(4)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    radarrp=where(radarr>0,radarr,1e-1)
    pmeas,=plot.plot(radarrp,beam_cmb[:,0],c='k',ls='-',label='Measured')
    ptheo,=plot.plot(radarrp,beam_cmb_t[:,0],c='k',ls='--',label='Theoretical')
    #pmod,=plot.plot(radarrp,beam_mod[:,0],c='k',ls=':',label=r'Modelled ($\gamma=1$)')
    for b in range(3):
        plot.plot(radarrp,beam_cmb[:,b],c=cols[b],ls='-',label=tit[b])
        plot.plot(radarrp,beam_cmb_t[:,b],c=cols[b],ls='--')
        #plot.plot(radarrp,beam_mod[:,b],c=cols[b],ls=':')
    plot.title(r'Measured vs Theoretical')
    plot.ylabel('Beam response')
    plot.xlabel('Angular distance [arcsec]')
    plot.xscale('linear')
    plot.yscale('linear')
    plot.ylim(0.5,1)
    plot.xlim(0,30)
    plot.legend(loc='upper right')
    #plot.legend(loc='lower left')
    plot.gca().lines.remove(pmeas)
    plot.gca().lines.remove(ptheo)
    #plot.gca().lines.remove(pmod)
    plot.savefig('../Plots/Beam_meas-vs-theory_mainbeam.png')
    plot.draw()
    
    plot.show()

if __name__ == "__main__":
    maincode()
