#### Test beam areas and aperture correction

from numpy import zeros,ones,arange,array,floor,ceil,log,log10,sum,max,min,pi,sqrt
from numpy import concatenate,inf,ones_like,zeros_like
import sys
#import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plot
from scipy import where
from scipy.interpolate import interp1d
#from scipy.interpolate import RectBivariateSpline
import argparse
import datetime

from beam import getbeammaps,getnewbeamprofile,comb_beam
from beam import get_effnu_az_newbeam_map,effbeam
from beam import arcsec2sr,measbeam_new,beam_fwhm
from beam import beamarea,beamarea_az

from rsrf import getrsrf,getapf,bandedge

def maincode():
    parser=argparse.ArgumentParser(description='Herschel-SPIRE Extended Emission Calibration')
    parser.add_argument('--brad',action='store',default=600.,type=float,dest='brad',help='Radius to integrate beam out to (arcsec). Default=350')
    parser.add_argument('--ind',action='store',default=0.85,type=float,dest='ind',help='Power law index with which beam FWHM changes with frequency (FWHM \\propto \\nu^{ind}). Default=0.65')
    parser.add_argument('--rsrf', action='store',default='M',dest='rsrftype',help='RSRF type to use [M=Measured | T=Top-hat]')
    parser.add_argument('--apeff', action='store',default='R',dest='apftype',help='Aperture Efficiency type to use [R=Real | U=Uniform]')
    parser.add_argument('--nepmod', action='store',default='ESA4',dest='nepmod',help='Model to use for Neptune spectrum ["ESA2"|"ESA4"]')
    parser.add_argument('--setmap', action='store_true',default=False,dest='setmap',help='Set to match effective frequencies to map areas')
    parser.add_argument('--donueff', action='store_true',default=False,dest='donueff',help='Set to recalculate effective frequency')

    args=parser.parse_args()
    brad=args.brad
    ind=args.ind    
    rsrftype=args.rsrftype
    apftype=args.apftype
    nepmod=args.nepmod
    setmap=args.setmap
    donueff=args.donueff
    
    bzlim=None
    bzval=0.
    aunit='arcsec'
    

    if nepmod=="ESA2":
        alpha_nep=array((1.26,1.39,1.45))
    elif nepmod=="ESA4":
        alpha_nep=array((1.29,1.42,1.47))
    else:
        print 'Invalid value for NEPMOD: "%s". Must be "ESA2" or "ESA4"'
        sys.exit()
    
    fsuff='_test-beam_Br%d_Ind%.2f_%s'%(int(brad),ind,nepmod)
    if setmap:
        fsuff=fsuff+'_MapAreas'
    else:
        fsuff=fsuff+'_ProfAreas'
    if rsrftype == 'T':
        fsuff=fsuff+'_noRSRF'
    if apftype == 'U':
        fsuff=fsuff+'_noApEff'

    mpl.rcParams['legend.fontsize']='medium'

    
    c=299792458. #speed of light
    #set band centres
    wlc=(250.e-6,350.e-6,500.e-6)
    wlc=array(wlc)
    nuc=c/wlc
    print nuc/1.e12
    band=arange(3)

    #########################
    ##  Source properties  ##
    #########################
    
    #Set source brightness and spectral index
    
    area0=arcsec2sr(array([426.,771.,1626.])) #area in steradians
    areapip=arcsec2sr(array([465.,822.,1768.]))
    fwhm0=array([17.6,23.9, 35.2])
    
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
    print 'Wavelength range: %.2f-%.2f microns' % (wlmin*1.e6,wlmax*1.e6)
    

    ####################
    ##  Read in RSRF  ##
    ####################
    
    #RSRFtype options:
    #  T=Top Hat
    #  M=Measured
    #rsrftype='M'
    print 'Getting RSRF of type "%s"...' % (rsrftype)
    rsrfarr=zeros((nnu,3))
    ilim=zeros((3,2))
    nulim=zeros((3,2))
    ilim_2=zeros((3,2))
    nulim_2=zeros((3,2))
    if rsrftype=="M":
        for b in band:
            rsrfarr[:,b]=getrsrf(rsrftype,b,nuarr)
            (ilim[b,:],nulim[b,:])=bandedge(nuarr,rsrfarr[:,b])
            (ilim_2[b,:],nulim_2[b,:])=bandedge(nuarr,rsrfarr[:,b],fact=2.)
            print 'Band %d limits: [%.2f:%.2f] GHz ([%.2f:%.2f] um)' % (b,nulim[b,0]/1.e9,nulim[b,1]/1.e9,c/nulim[b,1]*1.e6,c/nulim[b,0]*1.e6)
            #print 'Band %d limits: [%.2f:%.2f] um' % (b,c/nulim[b,1]*1.e6,c/nulim[b,0]*1.e6)
    elif rsrftype == 'T':
        nulim=array([[1.03e12,1.418e12],\
                   [0.741e12,1.008e12],\
                   [0.4906e12,0.7323e12]])
        #print nulim
        for b in band:
            rsrfarr_x = zeros(nnu)
            lim0=nulim[b,0]
            lim1=nulim[b,1]
            rsrfarr_x = where(nuarr > lim0,1.,0.)
            rsrfarr_x = where(nuarr > lim1,0.,rsrfarr_x)
            rsrfarr[:,b] = rsrfarr_x
            (ilim[b,:],nulim[b,:])=bandedge(nuarr,rsrfarr[:,b])
            (ilim_2[b,:],nulim_2[b,:])=bandedge(nuarr,rsrfarr[:,b],fact=2.)
    else:
        print 'Unknown value for RSRF. Must be one of [M (Measured)|T (Top Hat)'
        sys.exit()
    print nulim   
    
    #Aperture Efficiency options:
    # R=Real
    # U=Uniform
    #apftype='R'
    print 'Getting Aperture Efficiency of type "%s"...' % (apftype)
    
    if apftype == 'R':
        apfarr=zeros((nnu,3))
        for b in band:
            apfarr[:,b]=getapf(apftype,b,nuarr)
    elif apftype == 'U':
        apfarr=ones((nnu,3))
    else:
        print 'Unknown value for APEFF. Must be one of [R (Real)|U (Uniform)'
        sys.exit()
    
    #Compute 'extended source' rsrf
    rsrfearr=zeros((nnu,3))
    for b in band:
        rsrfearr[:,b]=rsrfarr[:,b]*(nuc[b]/nuarr)**(2.*ind)

    
    alpharr=arange(-4,5.5,0.5)
    nalph=alpharr.size
    alphapip=-1
    print 'Reading in beams...'
    
    radarr=arange(brad)
    nrad=radarr.size
    beam_scl=zeros((nrad,3))
    beam_fix=zeros((nrad,3))
    beam_cmb=zeros((nrad,3))
    areas_prof=zeros(3)
    bsw=array([[70,76,78,87],[95,103,110,130],[135,145,169,180]])
    for b in band:
        beam_scl[:,b],beam_fix[:,b]=getnewbeamprofile(b,radarr,bzlim=bzlim,bzval=bzval)
        beam_cmb[:,b]=comb_beam(beam_scl[:,b],beam_fix[:,b])
        areas_prof[b]=beamarea_az(radarr,beam_cmb[:,b],brad=brad)
    print 'Beam areas (from profile): [%.2f,%.2f,%.2f]'%(areas_prof[0],areas_prof[1],areas_prof[2])
    
    print 'Reading in beam maps...'
    #xarr_map=[]
    #yarr_map=[]
    #beam_map=[]
    #bwid=brad*1.1
    #bpix=1.
    areas_map=zeros(3)
    xarr_map,yarr_map,beam_map=getbeammaps(regrid=None)
    for b in band:
        #xarrin,yarrin,beamin=getbeam('M',b,bpix,bwid)
        #xarr_map.append(xarrin)
        #yarr_map.append(yarrin)
        #beam_map.append(beamin)
        areas_map[b]=beamarea(xarr_map[b],yarr_map[b],beam_map[b],brad=brad)
    print 'Beam Areas (from map) [sq. arcsec]: %.2f, %.2f, %.2f'%(areas_map[0], areas_map[1], areas_map[2])
    prof2map=areas_map/areas_prof
    print 'Prof -> Map factor: %.4f, %.4f, %.4f'%(prof2map[0],prof2map[1],prof2map[2])

    if setmap:
        print "Matching effective frequencies to map areas"
        area_effnu=areas_map
        prof2map[:]=1.
        pre_nueff=[1216.3e9,867.5e9,611.01e9]
    else:
        print "Matching effective frequencies to profile areas"
        area_effnu=None
        pre_nueff=[1222.0e9,871.56e9,613.17e9]
    
    
    print 'Calculating Effective frequencies (radial beam)...'
    aprec=0.0001
    if ind==0:
        nueff=nuc
    else:
        if donueff:
            nueff=get_effnu_az_newbeam_map(radarr,beam_scl,beam_fix,rsrfarr,nuarr,nuc,brad,alpha_nep,area_map=area_effnu,aprec=aprec,verbose=True,ind=ind)
        else:
            nueff=pre_nueff
    print 'Effective wavelengths:  [%.2f, %.2f, %.2f] um' % (c/nueff[0]*1.e6,c/nueff[1]*1.e6,c/nueff[2]*1.e6)
    #print 'Calculating Effective frequencies (beam map)...'
    #nueff2=get_effnu(beamx,beamy,beams,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=True)
    #print 'Effective wavelengths: [%.2f, %.2f, %.2f] um' % (c/nueff2[0]*1.e6,c/nueff2[1]*1.e6,c/nueff2[2]*1.e6)
    
    print 'Calculating effective beam area...'
    area_eff=zeros((nalph,3))
    fwhm_eff=zeros((nalph,3))
    area_pip=zeros(3)
    beam_eff=zeros((nrad,nalph,3))
    for b in band:
        print 'Band %d...'%(b)
        area_pip[b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alphapip,ind=ind)
        area_pip[b]=area_pip[b]*prof2map[b]
        for a in arange(nalph):
            beam_eff[:,a,b]=effbeam(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpharr[a],ind=ind)
            #area_eff[a,b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpharr[a],ind=ind)
            area_eff[a,b]=beamarea_az(radarr,beam_eff[:,a,b],brad=brad)
            fwhm_eff[a,b]=beam_fwhm(fwhm0[b],nuarr,rsrfarr[:,b],nueff[b],alpharr[a],ind)
        #area_eff[:,b]=area_eff[:,b]*prof2map[b]
        
    #plotareaa(alpharr,area_eff,areas_prof,fsuff,alpha_nep,unit=aunit)

    for a in range(nalph):
        print alpharr[a],area_eff[a,0],area_eff[a,1],area_eff[a,2]

    print 'Computing beam area profiles and aperture corrections...'
    #aperture and annulus radii in arcsec
    rad_ap=[22,30,42]
    rad_bg=[60,90]    
    omega_ap=zeros((nalph,3))
    omega_bg=zeros((nalph,3))
    size_ap=zeros(3)
    size_ap_est=zeros(3)
    size_bg=zeros(3)
    size_bg_est=zeros(3)
    apcorr_ap=zeros((nalph,3))
    apcorr_all=zeros((nalph,3))
    b_ones=ones_like(beam_cmb[:,0])
    area_rarr=arange(0,brad+10,10)
    nrad_a=area_rarr.size
    beam_eff_rad=zeros((nrad_a,nalph,3))
    for b in band:
        print 'Band %d...' % b
        #size within aperture and annulus [ int(2.pi.theta dtheta) ]
        size_ap[b]=beamarea_az(radarr,b_ones,brad=rad_ap[b])
        size_ap_est[b]=pi*rad_ap[b]**2
        size_bg[b]=beamarea_az(radarr,b_ones,brad=rad_bg[1])-beamarea_az(radarr,b_ones,brad=rad_bg[0])
        size_bg_est[b]=pi*(rad_bg[1]**2-rad_bg[0]**2)
        print 'Aperture:',[0,rad_ap[b]],size_ap[b],size_ap_est[b],size_ap_est[b]/size_ap[b]
        print 'Annulus: ',rad_bg,size_bg[b],size_bg_est[b],size_bg_est[b]/size_bg[b]
        for a in range(nalph):
            #calculate curve of growth (limited values of r)
            for r in range(nrad_a):
                beam_eff_rad[r,a,b]=beamarea_az(radarr,beam_eff[:,a,b],brad=area_rarr[r])
            
            #beam within aperture and annulus [ int(b(theta).2.pi.theta dtheta) ]
            omega_ap[a,b]=beamarea_az(radarr,beam_eff[:,a,b],brad=rad_ap[b])
            omega_bg[a,b]=beamarea_az(radarr,beam_eff[:,a,b],brad=rad_bg[1])-beamarea_az(radarr,beam_eff[:,a,b],brad=rad_bg[0])
        
            apcorr_ap[a,b]=area_eff[a,b]/omega_ap[a,b]
            apcorr_all[a,b]=area_eff[a,b]/(omega_ap[a,b] - omega_bg[a,b]*size_ap[b]/size_bg[b])
        #        for r in arange(nrad_a):
#            #areas_r[r,b]=beamarea_az(radarr,beam_sm[:,b],brad=area_rarr[r])
#            radx=float(area_rarr[r])
#            areas_r[r,b]=beamarea_az(radarr,beam_cmb[:,b],brad=radx)
#            #print '  Radius %.1f: Area %.2f'%(radx,areas_r[r,b])
        ##correct for profile->map effect

    print 'Areas (Aperture | Annulus)'
    for a in range(nalph):
        print '%.1f, %.4f, %.4f, %.4f | %.4f, %.4f, %.4f'%(alpharr[a],omega_ap[a,0],omega_ap[a,1],omega_ap[a,2],omega_bg[a,0],omega_bg[a,1],omega_bg[a,2])

    print 'Aperture corrections (Aperture only | Inc background)'
    for a in range(nalph):
        print '%.1f, %.4f, %.4f, %.4f | %.4f, %.4f, %.4f'%(alpharr[a],apcorr_ap[a,0],apcorr_ap[a,1],apcorr_ap[a,2],apcorr_all[a,0],apcorr_all[a,1],apcorr_all[a,2])
    
    plot.show()

def plotareaa(alpharr,areaeff,areas,fsuff,alpha_nep,unit='sr'):
    
    from numpy import array
    from beam import arcsec2sr    
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    area0=array([426.,771.,1626.])
    #alpha_nep=array((1.26,1.39,1.45))
    
    if unit=='sr':
        #convert to steradians
        area0=arcsec2sr(area0)*1.e8
        areas=arcsec2sr(areas)*1.e8
        areaeff=arcsec2sr(areaeff)*1.e8
        ytitu=r'[sr] $\times 10^8$'
    else: ytitu=r'[arcsec${^2}]$'
    
    xtickv=arange(-4,5,1)
    nt=xtickv.size
    xtickl=[]
    xtickbl=[]
    for n in arange(nt): 
        xtickl.append(str(xtickv[n]))
        xtickbl.append('')
    
    plot.figure(13)
    plot.clf()
    plot.subplot(3,1,1)
    plot.plot(alpharr,areaeff[:,0],'b-',label=r'$\Omega_\mathrm{eff}$')
    plot.axhline(y=areas[0],color='b',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    plot.axhline(y=area0[0],color='b',linestyle=':',label=r'$\Omega_0$')
    plot.axvline(x=alpha_nep[0],color='b',linestyle=':')
    #plot.ylabel('Beam Area '+ytitu)
    plot.xticks(xtickv,xtickbl)
    #plot.title('PSW')
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    plot.text(x0+0.05*xr,y0+0.05*yr,'PSW')
    plot.legend(loc='upper right',ncol=3)

    plot.subplot(3,1,2)
    plot.plot(alpharr,areaeff[:,1],'g-',label=r'$\Omega_\mathrm{eff}$')
    plot.axhline(y=areas[1],color='g',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    plot.axhline(y=area0[1],color='g',linestyle=':',label=r'$\Omega_0$')
    plot.axvline(x=alpha_nep[1],color='g',linestyle=':')
    plot.ylabel('Beam Area '+ytitu)
    plot.xticks(xtickv,xtickbl)
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    plot.text(x0+0.05*xr,y0+0.05*yr,'PMW')
    #plot.title('PMW')
    
    plot.subplot(3,1,3)
    plot.plot(alpharr,areaeff[:,2],'r-',label=r'$\Omega_\mathrm{eff}$')
    plot.axhline(y=areas[2],color='r',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    plot.axhline(y=area0[2],color='r',linestyle=':',label=r'$\Omega_0$')
    plot.axvline(x=alpha_nep[2],color='r',linestyle=':')
    #plot.ylabel('Beam Area '+ytitu)
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.xticks(xtickv,xtickl)
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    plot.text(x0+0.05*xr,y0+0.05*yr,'PLW')
    #plot.title('PLW')

    plot.savefig('../Outputs/beamarea_alpha'+fsuff+'_'+unit+'.png',transparent=False,dpi=300.)
    #plot.savefig('../Outputs/beamarea_alpha'+fsuff+'_'+unit+'.eps',transparent=False)
    
if __name__=="__main__":
    maincode()