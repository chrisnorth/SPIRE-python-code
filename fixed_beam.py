##########################################
## Test effect of fixed portion of beam ##
##########################################

import matplotlib.pyplot as plot
from numpy import arange,array,zeros,ones,shape
from scipy import where
import sys

from beam import getnewbeamprofile,comb_beam,beamarea_az,get_effnu_az_newbeam
from beam import beamarea_az_nu_new,avgbeam,measbeam_new,arcsec2sr
from rsrf import getrsrf,getapf,bandedge

def maincode():
   
    band=range(3)
    brad=500
    bzlim=None
    bzval=0.
    ind=0.8
    rsrftype='M'
    apftype='R'
    aunit='arcsec'
    
    c=299792458. #speed of light
    #set band centres
    wlc=(250.e-6,350.e-6,500.e-6)
    wlc=array(wlc)
    nuc=c/wlc
    print nuc/1.e12
    band=arange(3)
    
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
            print 'Frequency: %.1f [%.1f,%.1f]'%(nuc[b],nulim[b,0],nulim[b,1])
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

    print 'Reading in beams...'
    
    radarr=arange(brad)
    nrad=radarr.size
    beam_scl1=zeros((nrad,3))
    beam_fix1=zeros((nrad,3))
    beam_cmb1=zeros((nrad,3))
    beam_scl2=zeros((nrad,3))
    beam_fix2=zeros((nrad,3))
    beam_cmb2=zeros((nrad,3))
    areas1=zeros(3)
    areas2=zeros(3)
    bsw=array([[70,76,78,87],[95,103,110,130],[135,145,169,180]])
    for b in band:
        beam_scl1[:,b],beam_fix1[:,b]=getnewbeamprofile(b,radarr,bzlim=bzlim,bzval=bzval)
        beam_cmb1[:,b]=comb_beam(beam_scl1[:,b],beam_fix1[:,b])
        areas1[b]=beamarea_az(radarr,beam_cmb1[:,b],brad=brad)
        beam_scl2[:,b]=beam_cmb1[:,b]
        beam_fix2[:,b]=1.e-8
        beam_cmb2[:,b]=comb_beam(beam_scl2[:,b],beam_fix2[:,b])
        areas2[b]=beamarea_az(radarr,beam_cmb2[:,b],brad=brad)
        
    
    print 'Beam areas (1): [%.2f,%.2f,%.2f]'%(areas1[0],areas1[1],areas1[2])
    print 'Beam areas (2): [%.2f,%.2f,%.2f]'%(areas2[0],areas2[1],areas2[2])
    
    
    ##########################################
    ##  Calculate effective beam frequency  ##
    ##########################################
    #brad=350.
    print 'Calculating Effective frequencies (radial beam)...'
    aprec=0.0001
    if ind==0:
        nueff1=nuc
    else:
        nueff1=get_effnu_az_newbeam(radarr,beam_scl1,beam_fix1,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=True,ind=ind)
        nueff2=get_effnu_az_newbeam(radarr,beam_scl2,beam_fix2,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=True,ind=ind)
    print 'Effective wavelengths (1):  [%.2f, %.2f, %.2f] um' % (c/nueff1[0]*1.e6,c/nueff1[1]*1.e6,c/nueff1[2]*1.e6)
    print 'Effective wavelengths (2):  [%.2f, %.2f, %.2f] um' % (c/nueff2[0]*1.e6,c/nueff2[1]*1.e6,c/nueff2[2]*1.e6)
    #print 'Calculating Effective frequencies (beam map)...'
    #nueff2=get_effnu(beamx,beamy,beams,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=True)
    #print 'Effective wavelengths: [%.2f, %.2f, %.2f] um' % (c/nueff2[0]*1.e6,c/nueff2[1]*1.e6,c/nueff2[2]*1.e6)
    
    alpharr=arange(-4,5.5,0.5)
    nalph=alpharr.size
    
    alpha_nep=array((1.26,1.39,1.45))
    
    print 'Calculating beam area over frequency...'
    #Set frequency arrays for beam area calculations
    dnu_a=10
    area_bm1=beamarea_az_nu_new(nuarr,radarr,beam_scl1,beam_fix1,nueff1,ilim_2,dnu_a=dnu_a,verbose=False,ind=ind)
    area_bm2=beamarea_az_nu_new(nuarr,radarr,beam_scl2,beam_fix2,nueff2,ilim_2,dnu_a=dnu_a,verbose=False,ind=ind)

    area_bmd=zeros((nnu,3))
    for b in band:
        area_bmd[:,b]=where(area_bm1[:,b] != 0,(area_bm2[:,b]-area_bm1[:,b])/area_bm1[:,b],0.)

    avgnu_bm1=zeros(3)
    avgarea_bm1=zeros(3)
    avgnu_bm2=zeros(3)
    avgarea_bm2=zeros(3)
    print 'area_bm1:',area_bm1.shape
    print 'area_bm2:',area_bm2.shape
    print 'area_bmd:',area_bmd.shape
    for b in band:
        (avgnu_bm1[b],avgarea_bm1[b])=avgbeam(nuarr,area_bm1[:,b])
        (avgnu_bm2[b],avgarea_bm2[b])=avgbeam(nuarr,area_bm2[:,b])
    print 'Plotting Area vs. freq'
    plotareaf(nuarr,area_bm1,area_bm2,area_bmd,ilim,avgnu_bm1,avgarea_bm1,avgnu_bm2,avgarea_bm2,unit=aunit)

    print 'Calculating effective beam area...'
    area_eff1=zeros((nalph,3))
    area_eff2=zeros((nalph,3))
    area_effd=zeros((nalph,3))
    for b in band:
        print 'Band %d...'%(b)
        for a in arange(nalph):
            area_eff1[a,b]=measbeam_new(radarr,beam_scl1[:,b],beam_fix1[:,b],nuarr,rsrfarr[:,b],nueff1[b],alpharr[a],ind=ind)
            area_eff2[a,b]=measbeam_new(radarr,beam_scl2[:,b],beam_fix2[:,b],nuarr,rsrfarr[:,b],nueff2[b],alpharr[a],ind=ind)
        area_effd[:,b]=(area_eff1[:,b]-area_eff2[:,b])/area_eff1[:,b]
            

    plotareaa(alpharr,area_eff1,areas1,area_eff2,areas2,area_effd,unit=aunit)

    plot.show()

def plotareaa(alpharr,area_eff1,areas1,area_eff2,areas2,area_effd,unit='arcsec'):
    
    from numpy import array
    from beam import arcsec2sr    
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    area0=array([426.,771.,1626.])
    alpha_nep=array((1.26,1.39,1.45))
    
    if unit=='sr':
        #convert to steradians
        area0=arcsec2sr(area0)*1.e8
        areas1=arcsec2sr(areas1)*1.e8
        area_eff1=arcsec2sr(area_eff1)*1.e8
        areas2=arcsec2sr(areas2)*1.e8
        area_eff2=arcsec2sr(area_eff2)*1.e8
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
    plot.subplot(2,1,1)
    lpart,=plot.plot(alpharr,area_eff1[:,0],'k-',label=r'Part-scaled')
    lfull,=plot.plot(alpharr,area_eff1[:,0],'k--',label=r'Full-scaled')
    plot.plot(alpharr,area_eff1[:,0],'b-')
    plot.plot(alpharr,area_eff2[:,0],'b--')
    #plot.axhline(y=areas1[0],color='k',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    #plot.axhline(y=area0[0],color='k',linestyle=':',label=r'$\Omega_0$')
    plot.axvline(x=alpha_nep[0],color='b',linestyle=':')
    #plot.ylabel('Beam Area '+ytitu)
    plot.xticks(xtickv,xtickl)
    plot.ylim(0,2000)
    #plot.title('PSW')

    plot.plot(alpharr,area_eff1[:,1],'g-')
    plot.plot(alpharr,area_eff2[:,1],'g--')
    plot.axvline(x=alpha_nep[1],color='g',linestyle=':')
        
    plot.plot(alpharr,area_eff1[:,2],'r-')
    plot.plot(alpharr,area_eff2[:,2],'r--')
    plot.axvline(x=alpha_nep[2],color='r',linestyle=':')
    
    plot.ylabel('Effective beam Area '+ytitu)
    #plot.xlabel(r'Spectral index ($\alpha$)')
    plot.legend(loc='lower left',ncol=3)

    plot.gca().lines.remove(lpart)
    plot.gca().lines.remove(lfull)

    plot.subplot(2,1,2)
    
    plot.plot(alpharr,area_effd[:,0],'b-',label='PSW')
    plot.plot(alpharr,area_effd[:,1],'g-',label='PMW')
    plot.plot(alpharr,area_effd[:,2],'r-',label='PLW')
    plot.xticks(xtickv,xtickl)
    
    plot.ylabel('([Full] - [Part]) / [Part]')
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.legend(loc='lower right',ncol=3)

#    plot.plot(alpharr,area_eff1[:,1],'g-',label=r'$\Omega_\mathrm{eff}$ (part-scaled)')
#    plot.plot(alpharr,area_eff2[:,1],'g--',label=r'$\Omega_\mathrm{eff}$ (full-scaled)')
#    plot.axhline(y=areas1[1],color='k',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
#    plot.axhline(y=area0[1],color='k',linestyle=':',label=r'$\Omega_0$')
#    
#    plot.ylabel('Beam Area '+ytitu)
#    plot.xticks(xtickv,xtickbl)
#    (x0,x1)=plot.xlim()
#    (y0,y1)=plot.ylim()
#    xr=x1-x0
#    yr=y1-y0
#    plot.text(x0+0.05*xr,y0+0.05*yr,'PMW')
#    #plot.title('PMW')
#    
#    plot.subplot(3,1,3)
#    plot.plot(alpharr,area_eff1[:,2],'r-',label=r'$\Omega_\mathrm{eff}$ (part-scaled)')
#    plot.plot(alpharr,area_eff2[:,2],'r--',label=r'$\Omega_\mathrm{eff}$ (full-scaled)')
#    plot.axhline(y=areas1[2],color='k',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
#    plot.axhline(y=area0[2],color='k',linestyle=':',label=r'$\Omega_0$')
#    plot.axvline(x=alpha_nep[2],color='r',linestyle=':')
#    #plot.ylabel('Beam Area '+ytitu)
#    plot.xlabel(r'Spectral index ($\alpha$)')
#    plot.xticks(xtickv,xtickl)
#    (x0,x1)=plot.xlim()
#    (y0,y1)=plot.ylim()
#    xr=x1-x0
#    yr=y1-y0
#    plot.text(x0+0.05*xr,y0+0.05*yr,'PLW')
#    #plot.title('PLW')
    
    plot.savefig('../Outputs/beamarea_alpha_comp_'+unit+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamarea_alpha_comp_'+unit+'.eps',transparent=False)

def plotareaf(nuarr,area_bm1,area_bm2,area_bmd,ilim,avgnu1,avgarea1,avgnu2,avgarea2,unit='arcsec'):
    
    from numpy import array
    from beam import arcsec2sr    
    ############################################################
    #Plot beam area vs frequency

    c=299792458. #speed of light
    wl0=array([250.,350.,500.])*1.e-6
    nu0=(c/wl0)/1.e9 #in GHz
    area0=array([426.,771.,1626.])

    avgnu1=avgnu1/1.e9
    avgnu2=avgnu2/1.e9

    if unit=='sr':
        #convert to steradians
        area0=arcsec2sr(area0)*1.e8
        area_bm1=arcsec2sr(area_bm1)*1.e8
        avgarea1=arcsec2sr(avgarea1)*1.e8
        ytit='Beam Area [sr]'
    else: ytit=r'Beam Area [arcsec$^2$]'
    
    plot.figure(5)
    plot.clf()
    plot.subplot(2,1,1)
    
    ##plot area vs freq
    l1,=plot.plot(nuarr[ilim[0,0]:ilim[0,1]]/1.e9,area_bm1[ilim[0,0]:ilim[0,1],0],'k-',label='Part-scaled')
    plot.plot(nuarr[ilim[0,0]:ilim[0,1]]/1.e9,area_bm1[ilim[0,0]:ilim[0,1],0],'b-')
    plot.plot(nuarr[ilim[1,0]:ilim[1,1]]/1.e9,area_bm1[ilim[1,0]:ilim[1,1],1],'g-')
    plot.plot(nuarr[ilim[2,0]:ilim[2,1]]/1.e9,area_bm1[ilim[2,0]:ilim[2,1],2],'r-')
    
    l2,=plot.plot(nuarr[ilim[0,0]:ilim[0,1]]/1.e9,area_bm2[ilim[0,0]:ilim[0,1],0],'k--',label='Full-scaled')
    plot.plot(nuarr[ilim[0,0]:ilim[0,1]]/1.e9,area_bm2[ilim[0,0]:ilim[0,1],0],'b--')
    plot.plot(nuarr[ilim[1,0]:ilim[1,1]]/1.e9,area_bm2[ilim[1,0]:ilim[1,1],1],'g--')
    plot.plot(nuarr[ilim[2,0]:ilim[2,1]]/1.e9,area_bm2[ilim[2,0]:ilim[2,1],2],'r--')
    
    #plot nu-averaged beam areas
    plot.plot(nu0,area0,mec='k',marker='^',ls='None',label='Current beam areas')
    lavg1,=plot.plot(avgnu1[0],avgarea1[0],mec='k',marker='o',ls='None',label=r'$\nu$-avg area')
    #lavg2,=plot.plot(avgnu2[0],avgarea2[0],mec='k',mfc='k',marker='o',ls='None',label=r'$\nu$-avg area (full-scaled)')

    plot.plot(avgnu1[0],avgarea1[0],'bo',mec='b',mfc='b')
    plot.plot(avgnu2[0],avgarea2[0],'bo',mec='b',mfc='None')
    #plot.plot(array([0,avgnu1[0]]),array([avgarea1[0],avgarea1[0]]),'b:')
    #plot.plot(array([avgnu1[0],avgnu1[0]]),array([0,avgarea1[0]]),'b:')

    plot.plot(avgnu1[1],avgarea1[1],'go',mec='g',mfc='g')
    plot.plot(avgnu2[1],avgarea2[1],'go',mec='g',mfc='None')
    #plot.plot(array([0,avgnu1[1]]),array([avgarea1[1],avgarea1[1]]),'g:')
    #plot.plot(array([avgnu1[1],avgnu1[1]]),array([1,avgarea1[1]]),'g:')

    plot.plot(avgnu1[2],avgarea1[2],'ro',mec='r',mfc='r')
    plot.plot(avgnu2[2],avgarea2[2],'ro',mec='r',mfc='None')
    #plot.plot(array([0,avgnu1[2]]),array([avgarea1[2],avgarea1[2]]),'r:')
    #plot.plot(array([avgnu1[2],avgnu1[2]]),array([2,avgarea1[2]]),'r:')
    
    plot.ylabel(ytit)
    plot.xlabel('Frequency [GHz]')
    plot.legend(loc='upper right',numpoints=1,ncol=2)
    plot.gca().lines.remove(lavg1)
    #plot.gca().lines.remove(lavg2)
    plot.gca().lines.remove(l1)
    plot.gca().lines.remove(l2)
    
    plot.subplot(2,1,2)
    plot.plot(nuarr[ilim[0,0]:ilim[0,1]]/1.e9,area_bmd[ilim[0,0]:ilim[0,1],0],'b-',label='PSW')
    plot.plot(nuarr[ilim[1,0]:ilim[1,1]]/1.e9,area_bmd[ilim[1,0]:ilim[1,1],1],'g-',label='PMW')
    plot.plot(nuarr[ilim[2,0]:ilim[2,1]]/1.e9,area_bmd[ilim[2,0]:ilim[2,1],2],'r-',label='PLW')

    plot.ylabel('([Full] - [Part]) / [Part]')
    plot.xlabel('Frequency [GHz]')
    plot.legend(loc='upper right',numpoints=1,ncol=2)
    
    plot.savefig('../Outputs/beamarea_nu_comp_'+unit+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamarea_nu_com_'+unit+'.eps',transparent=False)

if __name__=="__main__":
    maincode()
