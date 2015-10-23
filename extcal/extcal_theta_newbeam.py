#from __future__ import division
#General python modules
#import pylab
from numpy import zeros,ones,arange,array,floor,ceil,log,log10,sum,max,min,pi,sqrt
from numpy import concatenate,inf
import sys
#import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plot
from scipy import where
#from scipy.interpolate import RectBivariateSpline
import argparse
import datetime

#Specific modules
from beam import getbeam,getnewbeamprofile,comb_beam,get_effnu_az_newbeam,beamprof_nu
from beam import arcsec2sr,avgbeam,modbeam_new,modbeam_area,measbeam_new
from beam import beamarea_az,beamarea_az_nu_new,beamarea_az_th_nu_new,beamarea_az_th
from rsrf import getrsrf,getapf,bandedge
from calc_conv import calc_k4,calc_kc,calc_k5

def maincode():
    parser=argparse.ArgumentParser(description='Herschel-SPIRE Extended Emission Calibration')
    parser.add_argument('--brad',action='store',default=350.,type=float,dest='brad',help='Radius to integrate beam out to (arcsec). Default=350')
    parser.add_argument('--ind',action='store',default=0.65,type=float,dest='ind',help='Power law index with which beam FWHM changes with frequency (FWHM \\propto \\nu^{ind}). Default=0.65')
    parser.add_argument('--bzlim',action='store',default=None,type=float,dest='bzlim',help='Zero-limit for beamfiles, below which the values are replaced by BZVAL. Default=None (i.e. don\'t set limit)')
    parser.add_argument('--bzval',action='store',default=0.,type=float,dest='bzval',help='Value with which to replace beam values below the zero-limit (BZLIM). Default=0')
    parser.add_argument('--aunit',action='store',default="sr",dest='aunit',help='Unit of area to use for area plots [arcsec|sr]. Default=sr')
    parser.add_argument('--newa', action='store_true',default=False,dest='newa',help='Set to use new area in colour-correction calculations')
    parser.add_argument('--rsrf', action='store',default='M',dest='rsrftype',help='RSRF type to use [M=Measured | T=Top-hat]')
    parser.add_argument('--apeff', action='store',default='R',dest='apftype',help='Aperture Efficiency type to use [R=Real | U=Uniform]')

    parser.add_argument('--prsrf', action='store_true',default=False,dest='prsrf',help='Set to plot RSRF and Ap. Eff.')
    parser.add_argument('--pblim', action='store_true',default=False,dest='pblim',help='Set to plot beam profiles at band limits')
    parser.add_argument('--pbcmb', action='store_true',default=False,dest='pbcmb',help='Set to plot combination of scaled and fixed beams')
    parser.add_argument('--pbcomp', action='store_true',default=False,dest='pbcomp',help='Set to plot beam comparison profiles (measured vs. theoretical)')
    parser.add_argument('--nopbeam', action='store_true',default=False,dest='nopbeam',help='Set to NOT plot beam maps (overrides PBEAM)')
    parser.add_argument('--pcorr', action='store_true',default=False,dest='pcorr',help='Set to plot colour correction factors')
    parser.add_argument('--pareaf', action='store_true',default=False,dest='pareaf',help='Set to plot beam area against frequency')
    parser.add_argument('--parear', action='store_true',default=False,dest='parear',help='Set to plot beam area against radius')
    parser.add_argument('--pareaa', action='store_true',default=False,dest='pareaa',help='Set to plot beam area against spectral index')
    parser.add_argument('--pareath', action='store_true',default=False,dest='pareath',help='Set to plot beam area against source FWHM')
    parser.add_argument('--pall', action='store_true',default=False,dest='pall',help='Plot all plots (set NOPBEAM to exclude beam maps')

    args=parser.parse_args()
    bzlim=args.bzlim
    bzval=args.bzval
    brad=args.brad
    pall=args.pall
    aunit=args.aunit
    newa=args.newa
    ind=args.ind    
    rsrftype=args.rsrftype
    apftype=args.apftype
    
    if pall:
        prsrf=True
        pblim=True
        pbcmb=True
        pbcomp=True
        pcorr=True
        pareaf=True
        parear=True
        pareaa=True
        pareath=True
    else:
        prsrf=args.prsrf
        pblim=args.pblim
        pbcmb=args.pbcmb
        pcorr=args.pcorr
        pareaf=args.pareaf
        parear=args.parear
        pareaa=args.pareaa
        pbcomp=args.pbcomp
        pareath=args.pareath
    
    fsuff='_theta_newBprofBS_Br%d_Ind%.2f'%(int(brad),ind)
    if bzlim != None:
        fsuff=fsuff+'_Bl%1g_Bv%1g'%(bzlim,bzval)
    if rsrftype == 'T':
        fsuff=fsuff+'_noRSRF'
    if apftype == 'U':
        fsuff=fsuff+'_noApEff'
    if newa == True:
        fsuff=fsuff+'_newArea'
    print 'File_suffix= '+fsuff

    file_summ='../Outputs/Summ'+fsuff+'.dat'
    fsumm=open(file_summ,'w')
    
    now = datetime.datetime.now()
    lsum='Date/Time: %s'%(now)    
    fsumm.write(lsum)
    
    lsum='\nInputs parameters:\n'
    fsumm.write(lsum)
    lsum='Beam radius: %d arcsec \n'%(brad)
    fsumm.write(lsum)
    lsum='Beam freq-dependent power-law index: %.2g \n'%(ind)
    fsumm.write(lsum)
    if bzlim:
        lsum='Beam zero-limit: %.2g \n'%(bzlim)
        fsumm.write(lsum)
        lsum='Beam zero-value: %.2g \n'%(bzval)
        fsumm.write(lsum)
    else:
        lsum='Beam zero-limit: None \n'
        fsumm.write(lsum)
        lsum='Beam zero-value: None \n'
        fsumm.write(lsum)
    if rsrftype=='T':
        lsum='RSRF: Tophat \n'
        fsumm.write(lsum)
    elif rsrftype=='M':
        lsum='RSRF: Measured \n'
        fsumm.write(lsum)
    if apftype=='R':
        lsum='Ap. Eff.: Real \n'
        fsumm.write(lsum)
    elif apftype=='U':
        lsum='Ap. Eff.: Uniform \n'
        fsumm.write(lsum)

    
    plotarr=array([prsrf,pblim,pcorr,pareaf,parear,pareaa,pbcomp,pareath])

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
        rsrfearr[:,b]=rsrfarr[:,b]*(nuc[b]/nuarr)**2
    
    if prsrf:
        plotrsrf(nuarr,rsrfarr,rsrfearr,apfarr,fsuff)

    lsum='\nPassbands:\n'
    fsumm.write(lsum)
    lsum='Bandedge (PSW) [THz]: %.4f, %.4f \n' % (nulim[0,0]*1.e-12,nulim[0,1]*1.e-12)
    fsumm.write(lsum)
    lsum='Bandedge (PMW) [THz]: %.4f, %.4f \n' % (nulim[1,0]*1.e-12,nulim[1,1]*1.e-12)
    fsumm.write(lsum)
    lsum='Bandedge (PLW) [THz]: %.4f, %.4f \n' % (nulim[2,0]*1.e-12,nulim[2,1]*1.e-12)
    fsumm.write(lsum)

    print "Writing RSRF, ApEff to files..."
    
    file_rsrf='../Outputs/RSRF'+fsuff+'.csv'
    file_apf='../Outputs/ApEff'+fsuff+'.csv'
    
    frsrf=open(file_rsrf,'w')
    fapf=open(file_apf,'w')

    line='#nu (THz), PSW, PMW, PLW \n'
    frsrf.write(line)
    fapf.write(line)

    for n in range(int(nnu)):
        line='%.5f , %.5g , %.5g , %.5g \n' %(nuarr[n]*1.e-12,rsrfarr[n,0],rsrfarr[n,1],rsrfarr[n,2])
        frsrf.write(line)
        line='%.5f , %.5g , %.5g , %.5g \n' %(nuarr[n]*1.e-12,apfarr[n,0],apfarr[n,1],apfarr[n,2])
        fapf.write(line)

    frsrf.close()
    fapf.close()
    
    ####################
    ##  Calculate K4  ##
    ####################
    
    print 'Calculating K4,Kc...'
    alpharr=arange(-4,5.5,0.5)
    nalph=alpharr.size
    K4=zeros((nalph,3))
    K4e=zeros((nalph,3))
    K4pip=zeros(3)
    #K4epip=zeros(3)
    alphapip=-1.
    
    Kc=zeros((nalph,3))
    Kce=zeros((nalph,3))
    
    for b in band:
        #Calculate K4 (pipeline)
        K4pip[b] = calc_k4(alphapip,rsrfarr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
        #Calculate K4 (extended,pipeline)
        #K4epip[b] = calc_k4(alphapip,rsrfearr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
        
        for a in arange(nalph):
            #Calculate K4 (point source)
            K4[a,b] = calc_k4(alpharr[a],rsrfarr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
            #Calculate K4 (fake-extended)    
            K4e[a,b] = calc_k4(alpharr[a],rsrfearr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
            
            #Calculate correction factor (point & fake-extended)        
            Kc[a,b] = calc_kc(alpharr[a],rsrfarr[:,b],apfarr[:,b],nuarr,dnu,nuc[b],alphapip)
            Kce[a,b] = calc_kc(alpharr[a],rsrfearr[:,b],apfarr[:,b],nuarr,dnu,nuc[b],alphapip)

            
    if pcorr:
        plotK4(alpharr,K4,K4e,fsuff)
        plotKc(alpharr,Kc,Kce,fsuff)

    #####################################################
    ## K4: Jy per Jy_meas for point source [f(alpha)]
    ## K4pip: Jy_pip per Jy_meas for point source (assumes alpha=-1)
    ## K4e: Jy per Jy_meas for extended source [f(alpha)]
    ## Kc: Jy per Jy_pip for point source [f(alpha)]
    ## Kce: Jy per Jy_pip for extended source [f(alpha)]
    ######################################################
    
    ####################
    ##  Read in beam  ##
    ####################
    
    #Beamtype options:
    #  G=Gaussian
    #  E=Elliptical Gaussian
    #  M=Measured beams
    #  T=Theoretical (modelled) beams)
    #  S=Spliced beams
    
    print 'Reading in beams...'
    
    radarr=arange(brad)
    nrad=radarr.size
    beam_scl=zeros((nrad,3))
    beam_fix=zeros((nrad,3))
    beam_cmb=zeros((nrad,3))
    areas=zeros(3)
    bsw=array([[70,76,78,87],[95,103,110,130],[135,145,169,180]])
    for b in band:
        beam_scl[:,b],beam_fix[:,b]=getnewbeamprofile(b,radarr,bzlim=bzlim,bzval=bzval)
        beam_cmb[:,b]=comb_beam(beam_scl[:,b],beam_fix[:,b])
        areas[b]=beamarea_az(radarr,beam_cmb[:,b],brad=brad)
    
    print 'Reading in beam maps...'
    xarr_map=[]
    yarr_map=[]
    beam_map=[]
    bwid=brad*1.1
    bpix=1.
    for b in band:
        xarr_map.append(getbeam('M',b,bpix,bwid,bzlim=bzlim,bzval=bzval,verbose=None))
        yarr_map.append(getbeam('M',b,bpix,bwid,bzlim=bzlim,bzval=bzval,verbose=None))
        beam_map.append(getbeam('M',b,bpix,bwid,bzlim=bzlim,bzval=bzval,verbose=None))
    print 'Beam areas: [%.2f,%.2f,%.2f]'%(areas[0],areas[1],areas[2])
    
    lsum='\nBeams:\n'
    fsumm.write(lsum)
    lsum='Beam Areas [sq. arcsec]: %.2f, %.2f, %.2f \n'%(areas[0], areas[1], areas[2])
    fsumm.write(lsum)
    
    if pblim:
        beam_scl_lim=zeros((nrad,3,3)) #rad,nu,band
        beam_cmb_lim=zeros((nrad,3,3)) #rad,nu,band
        for b in band:
            nulist=array([nulim[b,0],nuc[b],nulim[b,1]])
            beam_cmb_lim[:,:,b]=beamprof_nu(radarr,beam_scl[:,b],beam_fix[:,b],nulist,nuc[b],ind=ind)
            beam_scl_lim[:,:,b]=beamprof_nu(radarr,beam_scl[:,b],zeros(nrad),nulist,nuc[b],ind=ind)
        
        plotblim(radarr,beam_scl_lim,beam_cmb_lim,beam_fix,brad,fsuff)

    #sys.exit()

    if pbcmb:
        plotbcomb(radarr,beam_scl,beam_fix,beam_cmb,bsw,fsuff)
    
    #sys.exit()
    
    ##########################################
    ##  Calculate effective beam frequency  ##
    ##########################################
    #brad=350.
    print 'Calculating Effective frequencies (radial beam)...'
    aprec=0.0001
    if ind==0:
        nueff=nuc
    else:
        nueff=get_effnu_az_newbeam(radarr,beam_scl,beam_fix,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=True,ind=ind)
    print 'Effective wavelengths:  [%.2f, %.2f, %.2f] um' % (c/nueff[0]*1.e6,c/nueff[1]*1.e6,c/nueff[2]*1.e6)
    #print 'Calculating Effective frequencies (beam map)...'
    #nueff2=get_effnu(beamx,beamy,beams,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=True)
    #print 'Effective wavelengths: [%.2f, %.2f, %.2f] um' % (c/nueff2[0]*1.e6,c/nueff2[1]*1.e6,c/nueff2[2]*1.e6)
    
    lsum='Effective frequencies [THz]: %.5f, %.5f, %.5f \n'%(nueff[0]*1.e-12,nueff[1]*1.e-12,nueff[2]*1.e-12)    
    fsumm.write(lsum)
    lsum='Effective/nominal frequencies: %.5f, %.5f, %.5f \n'%(nueff[0]/nuc[0], nueff[1]/nuc[1],nueff[2]/nuc[2])
    fsumm.write(lsum)
    
    ### thetarr is in FWHM
    thetarr=concatenate(([1.e-8],arange(1.,10.,0.2)*1.e-3, \
                        arange(1.,10.,0.2)*1.e-2, \
                        arange(1.,10.,0.2)*1.e-1, \
                        arange(1.,10.,0.2), \
                        arange(1.,10.,0.2)*10.,\
                        arange(1.,11.,0.2)*100.,[10000]))
    ithi=array(where(thetarr == max(thetarr)))[0]
    print 'ithi= %d'%ithi
    #thetarr=array((1.e-3,1000.))
    #ithi=1.
    
    nth=thetarr.size

    #for t in range(nth):
    #    print t,thetarr[t]
    #sys.exit()

    srcarea1=arcsec2sr(pi*(thetarr/(2.*sqrt(log(2.))))**2)
    #srcarea2=zeros(nth)
    #for t in arange(nth):
    #    srcarea2[t]=srcarea(thetarr[t])
    
    nalph=alpharr.size
    nth=thetarr.size
    arr_a=zeros((nalph,nth))
    arr_t=zeros((nalph,nth))
    for a in arange(nalph):
        for t in arange(nth):
            arr_a[a,t]=alpharr[a]
            arr_t[a,t]=thetarr[t]
    

    alpha_nep=array((1.26,1.39,1.45))
    
    print 'Calculating beam area over frequency...'
    #Set frequency arrays for beam area calculations
    dnu_a=10
    area_bmth=beamarea_az_th_nu_new(nuarr,thetarr,radarr,beam_scl,beam_fix,nueff,ilim_2,dnu_a=dnu_a,verbose=True,ind=ind)
    area_bm=beamarea_az_nu_new(nuarr,radarr,beam_scl,beam_fix,nueff,ilim_2,dnu_a=dnu_a,verbose=False,ind=ind)

    if pareaf:
        avgnu_bm=zeros((nth,3))
        avgarea_bm=zeros((nth,3))
        print 'area_bmth:',area_bmth.shape
        for t in arange(nth):
            for b in band:
                (avgnu_bm[t,b],avgarea_bm[t,b])=avgbeam(nuarr,area_bmth[:,t,b])
        print 'Plotting Area vs. freq'
        plotareaf(nuarr,area_bmth,ilim_2,avgnu_bm,avgarea_bm,ithi,fsuff,unit=aunit)

    print 'Calculating effective beam area...'
    area_eff=zeros((nalph,3))
    for b in band:
        print 'Band %d...'%(b)
        for a in arange(nalph):
            area_eff[a,b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpharr[a],ind=ind)
    filebeam='../Outputs/beamarea'+fsuff+'.csv'
    fbm=open(filebeam,'w')
    line='#alpha , PSW , PMW , PLW \n'
    fbm.write(line)
    for a in range(nalph):
        line='%.2f, %.5f, %.5f, %.5f\n'%(alpharr[a],area_eff[a,0],area_eff[a,1],area_eff[a,2])
        fbm.write(line)
    fbm.close()
    
    if pareaa:
        print 'Plotting Area vs. alpha'
        plotareaa(alpharr,area_eff,areas,fsuff,unit=aunit)
    
    #convert area to sr    
    area_bmth_sr=arcsec2sr(area_bmth)
    area_bm_sr=arcsec2sr(area_bm)
    area_eff_sr=arcsec2sr(area_eff)
    
    if pareath:
        print 'Plotting Area vs. theta'
        areas_th=zeros((3,nth))
        for t in arange(nth):
            for b in band:
                areas_th[b,t]=beamarea_az_th(radarr,beam_cmb[:,b],thetarr[t])
            print '%.3g": %.2g , %.2g , %.2g'%(thetarr[t],areas_th[0,t],areas_th[1,t],areas_th[2,t])
        plotareath(thetarr,areas_th,fsuff,unit=aunit)                
    
    if parear or pbcomp:
        print 'Computing beam area profiles...'
        area_rarr=arange(0,brad+10,10)
        nrad_a=area_rarr.size
        areas_r=zeros((nrad_a,3))
        for b in band:
            print 'Band %d...' % b
            for r in arange(nrad_a):
                #areas_r[r,b]=beamarea_az(radarr,beam_sm[:,b],brad=area_rarr[r])
                radx=float(area_rarr[r])
                areas_r[r,b]=beamarea_az(radarr,beam_cmb[:,b],brad=radx)
                #print '  Radius %.1f: Area %.2f'%(radx,areas_r[r,b])
            print min(beam_cmb[:,b]),max(beam_cmb[:,b])
        print 'Areas: [%.2f , %.2f , %.2f] sq-arcsec' % (areas_r[nrad_a-1,0],areas_r[nrad_a-1,1],areas_r[nrad_a-1,2])
        
        file_brad='../Outputs/beamarea_rad'+fsuff+'.csv'
        fbrad=open(file_brad,'w')
        line='# Cumulative beam area (square arcmin)\n'
        fbrad.write(line)
        line='#Radius (arcmin), PSW, PMW, PLW\n'
        fbrad.write(line)
        for r in range(nrad_a):
            line='%.1f, %.2f, %.2f, %.2f\n'%(area_rarr[r],areas_r[r,0],areas_r[r,1],areas_r[r,2])
            fbrad.write(line)
        fbrad.close()
        
        if parear:
            plotarear(area_rarr,areas_r,fsuff,unit=aunit)

        plot.figure(30)
        plot.plot(radarr,beam_scl[:,0],'b-')
        plot.plot(radarr,beam_scl[:,1],'g-')
        plot.plot(radarr,beam_scl[:,2],'r-')
        plot.plot(radarr,beam_fix[:,0],'b--')
        plot.plot(radarr,beam_fix[:,1],'g--')
        plot.plot(radarr,beam_fix[:,2],'r--')
        plot.yscale('log')
        
    if pbcomp:
        print 'Comparing modelled and measured areas...'
        beam_mod=zeros((nrad,3))
        beam_comp=zeros((nrad,3))
        #areas_mod=zeros((nrad_a,3))
        area_rarr2=arange(1,101,1)
        nrad_a2=area_rarr2.size
        areas_nep=zeros((nrad_a2,3))
        areas_mod=zeros((nrad_a2,3))
        areas_rel=zeros((nrad_a2,3))
        for b in band:
            print 'Band %d...' % b
            beam_mod[:,b]=modbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpha_nep[b],ind=ind)
            beam_comp[:,b]=(beam_mod[:,b]-beam_cmb[:,b])/areas[b]
            areas_nep[:,b]=modbeam_area(radarr,beam_cmb[:,b],area_rarr2)
            areas_mod[:,b]=modbeam_area(radarr,beam_mod[:,b],area_rarr2)
            areas_rel[:,b]=(areas_mod[:,b]-areas_nep[:,b])/areas_nep[:,b]
            #areas_rel[:,b]=(modbeam_area[:,b]-areas_r[:,b])/areas_r[:,b]
            #print 'areas (Nep,Mod,Rel):',areas_nep[:,b],areas_mod[:,b],areas_rel[:,b]
            print 'min/max',b,min(areas_rel[:,b]),max(areas_rel[:,b])
        plotbeamcomp(radarr,beam_cmb,beam_mod,area_rarr2,areas_rel,fsuff)
    
    ############################
    ##  Calculate K5 and Kc2  ##
    ############################
    
    K5=zeros((nalph,3))
    Kc2=zeros((nalph,3))
    Kc2b=zeros((nalph,3))
    K6=zeros((nalph,nth,3))
    Kc3=zeros((nalph,nth,3))
    Kc3t=zeros((nalph,nth,3))
    
    for b in band:
        print 'band %d'%b
        for a in arange(nalph):
            ##calculate fully extended sources
            K5[a,b]=calc_k5(alpharr[a],nuarr,rsrfarr[:,b],apfarr[:,b],area_bm_sr[:,b],nuc[b])
            Kc2[a,b]=K5[a,b]/K4pip[b]
            Kc2b[a,b]=Kc2[a,b]*area_eff_sr[a,b]
            for t in arange(nth):
                K6[a,t,b]=calc_k5(alpharr[a],nuarr,rsrfarr[:,b],apfarr[:,b],area_bmth_sr[:,t,b],nuc[b])
                #print '(alpha %d,theta %.3f): %g , %f'%(alpharr[a],thetarr[t],K5[a,t,b],log10(K5[a,t,b]))
                Kc3[a,t,b]=K6[a,t,b]/(K4pip[b])#/arcsec2sr(areas[b]))
                Kc3t[a,t,b]=K6[a,t,b]/K4pip[b]*srcarea1[t]
     
    ##################################################
    ## K5: Jy/sr per Jy_meas for fully extended source [f(alpha)]
    ## Kc2: Jy/sr per Jy_pip for fully extended source [f(alpha)]
    ## Kc2b: Jy/beam per Jy_pip for fully extended source [f(alpha)]
    ## K6: Jy/sr per Jy_meas for partial source [f(alpha,theta)]
    ## Kc3: Jy/sr per Jy_pip for partial source [f(alpha,theta)]
    ## Kc3t: Jy per Jy_pip for partial source [f(alpha,theta)]
    ##################################################
    
    if pcorr:
        #plotK5(arr_a,arr_t,K5,fsuff)
        #plotKc2(arr_a,arr_t,Kc2,fsuff)

        aplot=14 #alpha=3
        iths=0
        thsmall=thetarr[iths]
        ithb=nth-1
        thbig=thetarr[ithb]
        
        plotK6_va(alpharr,thetarr,K6,K5,aplot,fsuff)
        plotKc3_va(alpharr,thetarr,Kc3,Kc2,aplot,fsuff)
        plotKc3t_va(alpharr,thetarr,Kc3t,Kc,aplot,fsuff)
        plotKc_all(alpharr,Kc,Kce,Kc3[:,ithb,:],Kc3t[:,0,:],fsuff)
        plotKcKc3t_thsmall(alpharr,Kc,Kc3t[:,iths,:],thsmall,fsuff)
        plotKc2Kc3_thbig(alpharr,Kc2,Kc3[:,ithb,:],thbig,fsuff)
        #plotsrcarea(thetarr,srcarea1,srcarea2,fsuff)
        #plotKc3(arr_a,arr_t,Kc3,fsuff)        
        plotKc2(alpharr,Kce,Kc2,area0,fsuff)
        #plotKc2rel(alpharr,Kce,Kc2,area0,fsuff)
    
    print Kc3[:,nth-1,0]
    print Kc2[:,0]
    print (Kc3[:,nth-1,0]-Kc2[:,0])/Kc2[:,0]

    Gfact=zeros((nalph,3))
    for b in band:
        Gfact[:,b]=(Kc[:,b]/arcsec2sr(areas[b]))/Kc2[:,b]
    
    lsum='\nOutputs:\n'
    fsumm.write(lsum)
    asum1=6 #alpha=-1
    asum2=14 #alpha=3
    asum3=12 #alpha=2
    
    oldDiff=100.*(Kce-Kc2b)/Kc2b
    #lsum='alpha 1: %.1f \n'%(alpharr[alpharr[asum1]])
    #fsumm.write(line)
    #lsum='alpha 2: %.1f \n'%(alpharr[asum2])
    #fsumm.write(line)
    lsum='K4(%d) [KMonP]: %.5f, %.5f, %.5f \n'%(alpharr[asum1],K4[asum1,0],K4[asum1,1],K4[asum1,2])
    fsumm.write(lsum)
    lsum='K4e(%d) [KMonE-old]: %.5f, %.5f, %.5f \n'%(alpharr[asum1],K4e[asum1,0],K4e[asum1,1],K4e[asum1,2])
    fsumm.write(lsum)
    lsum='Kc(%d) [KColP]: %.5f, %.5f, %.5f \n'%(alpharr[asum1],Kc[asum1,0],Kc[asum1,1],Kc[asum1,2])
    fsumm.write(lsum)
    lsum='Kce(%d) [KColE-old]: %.5f, %.5f, %.5f \n'%(alpharr[asum1],Kce[asum1,0],Kce[asum1,1],Kce[asum1,2])
    fsumm.write(lsum)
    lsum='K5(%d) [KMonE]: %.5g, %.5g, %.5g \n'% (alpharr[asum1],K5[asum1,0], K5[asum1,1], K5[asum1,2])
    fsumm.write(lsum)
    lsum='Kc2(%d) [KColE]: %.5g, %.5g, %.5g \n'% (alpharr[asum1],Kc2[asum1,0], Kc2[asum1,1], Kc2[asum1,2])
    fsumm.write(lsum)
    lsum='Kc2b(%d) [KColE*Beam]: %.5f, %.5f, %.5f \n'% (alpharr[asum1],Kc2b[asum1,0], Kc2b[asum1,1], Kc2b[asum1,2])
    fsumm.write(lsum)
    lsum='Gfact(%d) [(KColP/Beam)/KColE]: %.5f, %.5f, %.5f \n'% (alpharr[asum1],Gfact[asum1,0],Gfact[asum1,1],Gfact[asum1,2])
    fsumm.write(lsum)
    lsum='(Old-New)/New(%d) [KColE-old/KColE*Beam] (%%): %.2f, %.2f, %.2f \n'%(alpharr[asum1],oldDiff[asum1,0],oldDiff[asum1,1],oldDiff[asum1,2])
    fsumm.write(lsum)
    lsum='K4(%d) [KMonP]: %.5f, %.5f, %.5f \n'%(alpharr[asum2],K4[asum2,0],K4[asum2,1],K4[asum2,2])
    fsumm.write(lsum)
    lsum='K4e(%d) [KMonE-old]: %.5f, %.5f, %.5f \n'%(alpharr[asum2],K4e[asum2,0],K4e[asum2,1],K4e[asum2,2])
    fsumm.write(lsum)
    lsum='Kc(%d) [KColP]: %.5f, %.5f, %.5f \n'%(alpharr[asum2],Kc[asum2,0],Kc[asum2,1],Kc[asum2,2])
    fsumm.write(lsum)
    lsum='Kce(%d) [KColE-old]: %.5f, %.5f, %.5f \n'%(alpharr[asum2],Kce[asum2,0],Kce[asum2,1],Kce[asum2,2])
    fsumm.write(lsum)
    lsum='K5(%d) [KMonE]: %.5g, %.5g, %.5g \n'% (alpharr[asum2],K5[asum2,0], K5[asum2,1], K5[asum2,2])
    fsumm.write(lsum)
    lsum='Kc2(%d) [KColE]: %.5g, %.5g, %.5g \n'% (alpharr[asum2],Kc2[asum2,0], Kc2[asum2,1], Kc2[asum2,2])
    fsumm.write(lsum)
    lsum='Kc2b(%d) [KColE*beam]: %.5f, %.5f, %.5f \n'% (alpharr[asum2],Kc2b[asum2,0], Kc2b[asum2,1], Kc2b[asum2,2])
    fsumm.write(lsum)
    lsum='Gfact(%d) [(KColP/Beam)/KColE]: %.5f, %.5f, %.5f \n'% (alpharr[asum2],Gfact[asum2,0],Gfact[asum2,1],Gfact[asum2,2])
    fsumm.write(lsum)
    lsum='(Old-New)/New(%d) [KColE-old/KColE*Beam] (%%): %.2f, %.2f, %.2f \n'%(alpharr[asum2],oldDiff[asum2,0],oldDiff[asum2,1],oldDiff[asum2,2])
    fsumm.write(lsum)
    lsum='Gfact(%d) [(KColP/Beam)/KColE]: %.5f, %.5f, %.5f \n'% (alpharr[asum3],Gfact[asum3,0],Gfact[asum3,1],Gfact[asum3,2])
    fsumm.write(lsum)
    
    ######################
    ##  Write to files  ##
    ######################
    print 'Writing K4,K4e,Kc,Kce,K5,Kc2,Kc2b to files...'
    
    file_K4='../Outputs/K4'+fsuff+'.csv'
    file_K4e='../Outputs/K4e'+fsuff+'.csv'
    file_Kc='../Outputs/Kc'+fsuff+'.csv'
    file_Kce='../Outputs/Kce'+fsuff+'.csv'
    file_K5='../Outputs/K5'+fsuff+'.csv'
    file_Kc2='../Outputs/Kc2'+fsuff+'.csv'
    file_Kc2b='../Outputs/Kc2b'+fsuff+'.csv'
    file_tex='../Docs/Tables/Kc_Kce'+fsuff+'.tex'
    file_tex2='../Docs/Tables/Kc2_Kc2b'+fsuff+'.tex'
    file_tex3='../Docs/Tables/Kc_Kc2'+fsuff+'.tex'
    
    fK4=open(file_K4,'w')
    fK4e=open(file_K4e,'w')
    fKc=open(file_Kc,'w')
    fKce=open(file_Kce,'w')
    fK5=open(file_K5,'w')
    fKc2=open(file_Kc2,'w')
    fKc2b=open(file_Kc2b,'w')
    ftex=open(file_tex,'w')
    ftex2=open(file_tex2,'w')
    ftex3=open(file_tex3,'w')
    
    line='#alpha , PSW , PMW , PLW \n'
    fK4.write(line)
    fK4e.write(line)
    fKc.write(line)
    fKce.write(line)
    fK5.write(line)
    fKc2.write(line)
    fKc2b.write(line)

    line='\\begin{tabular}{c|ccc|ccc}\n'
    ftex.write(line)
    ftex2.write(line)
    line='\\begin{tabular}{r|l|l}\n'
    ftex3.write(line)
    line='& \\multicolumn{3}{c|}{$K_\\mathrm{c,pt}$} & \\multicolumn{3}{c}{$K_\\mathrm{c,ext}$} \\\\\n'
    ftex.write(line)
    line='& \\multicolumn{3}{c|}{$K_\\mathrm{c2} \\times 10^{-7}$} & \\multicolumn{3}{c}{$K_\\mathrm{c2}^\\mathrm{bm}$} \\\\\n'
    ftex2.write(line)
    line='& \\multicolumn{3}{c|}{Point source ($K_\\mathrm{CP}$)} & \\multicolumn{3}{c}{Extended source ($K_\\mathrm{CE}$)} \\\\\n'
    ftex3.write(line)
    line='$\\alpha$ & PSW & PMW & PLW & PSW & PMW & PLW \\\\\n'
    ftex.write(line)
    ftex2.write(line)
    ftex3.write(line)
    line='\\hline\n'
    ftex.write(line)
    ftex2.write(line)
    ftex3.write(line)
    
    Kc2_p=Kc2
    Kc2_exp=floor(log10(Kc2_p))
    Kc2_man=Kc2_p/10.**Kc2_exp
    Kc2b_p=Kc2b
    Kc2b_exp=floor(log10(Kc2b_p))
    Kc2b_man=Kc2b_p/10.**Kc2b_exp

        
    for a in arange(nalph):
        line='%.1f , %.9f , %.9f , %.9f \n' % (alpharr[a],K4[a,0],K4[a,1],K4[a,2])
        fK4.write(line)
        line='%.1f , %.9f , %.9f , %.9f \n' % (alpharr[a],K4e[a,0],K4e[a,1],K4e[a,2])
        fK4e.write(line)
        line='%.1f , %.9f , %.9f , %.9f \n' % (alpharr[a],Kc[a,0],Kc[a,1],Kc[a,2])
        fKc.write(line)
        line='%.1f , %.9f , %.9f , %.9f \n' % (alpharr[a],Kce[a,0],Kce[a,1],Kce[a,2])
        fKce.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (alpharr[a],K5[a,0],K5[a,1],K5[a,2])
        fK5.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (alpharr[a],Kc2[a,0],Kc2[a,1],Kc2[a,2])
        fKc2.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (alpharr[a],Kc2b[a,0],Kc2b[a,1],Kc2b[a,2])
        fKc2b.write(line)
        line='%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n' % (alpharr[a],Kc[a,0],Kc[a,1],Kc[a,2],Kce[a,0],Kce[a,1],Kce[a,2])
        ftex.write(line)
        line='%.1f & %.4g & %.4g & %.4g & %.4f & %.4f & %.4f \\\\ \n' % (alpharr[a],Kc2[a,0]/1.e7,Kc2[a,1]/1.e7,Kc2[a,2]/1.e7,Kc2b[a,0],Kc2b[a,1],Kc2b[a,2])
        ftex2.write(line)
        line='%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n' % (alpharr[a],Kc[a,0],Kc[a,1],Kc[a,2],Kc2[a,0]/1.e6,Kc2[a,1]/1.e6,Kc2[a,2]/1.e6)
        ftex3.write(line)
        
    fK4.close()
    fK4e.close()
    fKc.close()
    fKce.close()
    fK5.close()
    fKc2.close()
    fKc2b.close()
    
    line='\\end{tabular}\n'
    ftex.write(line)
    ftex2.write(line)
    ftex3.write(line)
    ftex.close()
    ftex2.close()
    ftex3.close()


    print 'Writing K6,Kc3,Kc3i to files...'
    
    ithp=array((0,28,32,37,38,39,41,46,55,56))
    nthp=ithp.size

    th_p=thetarr
    th_exp=floor(log10(th_p))
    th_man=th_p/10.**th_exp
    Kc3_p=Kc3
    Kc3_exp=floor(log10(Kc3_p))
    Kc3_man=Kc3_p/10.**Kc3_exp
    Kc3t_p=Kc3t
    Kc3t_exp=floor(log10(Kc3t_p))
    Kc3t_man=Kc3t_p/10.**Kc3t_exp


    bname=['PSW','PMW','PLW']
    for b in band:
        #file names
        file_K6_at='../Outputs/K6_'+bname[b]+fsuff+'.csv'
        file_Kc3_at='../Outputs/Kc3_'+bname[b]+fsuff+'.csv'
        file_Kc3t_at='../Outputs/Kc3t_'+bname[b]+fsuff+'.csv'
        file_Kc3_tex='../Docs/Tables/Kc3_'+bname[b]+fsuff+'.tex'
        file_Kc3t_tex='../Docs/Tables/Kc3t_'+bname[b]+fsuff+'.tex'
        #open files        
        fK6at=open(file_K6_at,'w')
        fKc3at=open(file_Kc3_at,'w')
        fKc3tat=open(file_Kc3t_at,'w')
        fKc3_tex=open(file_Kc3_tex,'w')
        fKc3t_tex=open(file_Kc3t_tex,'w')

        #column headings
        fK6at.write('#theta/alpha')
        fKc3at.write('#theta/alpha')
        fKc3tat.write('#theta/alpha')

        for a in arange(nalph):
            fK6at.write(',%.1f'%(alpharr[a]))
            fKc3at.write(',%.1f'%(alpharr[a]))
            fKc3tat.write(',%.1f'%(alpharr[a]))
        fK6at.write('\n') #end-of-line
        fKc3at.write('\n') #end-of-line
        fKc3tat.write('\n') #end-of-line

        #rows
        for t in arange(nth):
            fK6at.write('%.3f'%thetarr[t]) #row headings
            fKc3at.write('%.3f'%thetarr[t]) #row headings
            fKc3tat.write('%.3f'%thetarr[t]) #row headings
            for a in arange(nalph):
                fK6at.write(',%.9g'%K6[a,t,b]) #conv factors
                fKc3at.write(',%.9g'%Kc3[a,t,b]) #conv factors
                fKc3tat.write(',%.9g'%Kc3t[a,t,b]) #conv factors
            fK6at.write('\n')
            fKc3at.write('\n')
            fKc3tat.write('\n')


        line='\\begin{tabular}{c|ccccccccccc}\n'
        fKc3_tex.write(line)
        fKc3t_tex.write(line)
        line='& \\multicolumn{10}{c}{Source FWHM, $\\theta_0$ (arcsec)} \\\\\n'
        fKc3_tex.write(line)
        fKc3t_tex.write(line)
        fKc3_tex.write('$\\alpha$ ')
        fKc3t_tex.write('$\\alpha$ ')
        for p in range(nthp):
            if thetarr[ithp[p]] < 0.1:
                fKc3_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
                fKc3t_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
            elif thetarr[ithp[p]] >= 10000.:
                fKc3_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
                fKc3t_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
            else:
                fKc3_tex.write('& %.1f '%thetarr[ithp[p]])
                fKc3t_tex.write('& %.1f '%thetarr[ithp[p]])
        fKc3_tex.write('\\\\\n')
        fKc3t_tex.write('\\\\\n')
        fKc3_tex.write('\\hline\n')
        fKc3t_tex.write('\\hline\n')

        for a in range(nalph):
            fKc3_tex.write('%.1f '%alpharr[a])
            fKc3t_tex.write('%.1f '%alpharr[a])
            for p in range(nthp):
                if Kc3[a,ithp[p],b] >= 10:
                    fKc3_tex.write('& $%.2f \\times 10^{%d}$ '%(Kc3_man[a,ithp[p],b],Kc3_exp[a,ithp[p],b]))
                elif Kc3[a,ithp[p],b] <= 0.01:
                    fKc3_tex.write('& $%.2f \\times 10^{%d}$ '%(Kc3_man[a,ithp[p],b],Kc3_exp[a,ithp[p],b]))
                else:
                    fKc3_tex.write('& %.4f '%Kc3[a,ithp[p],b])
                if Kc3t[a,ithp[p],b] >= 10:
                    fKc3t_tex.write('& $%.2f \\times 10^{%d}$ '%(Kc3t_man[a,ithp[p],b],Kc3t_exp[a,ithp[p],b]))
                elif Kc3t[a,ithp[p],b] <= 0.01:
                    fKc3t_tex.write('& $%.2f \\times 10^{%d}$ '%(Kc3t_man[a,ithp[p],b],Kc3t_exp[a,ithp[p],b]))
                else:
                    fKc3t_tex.write('& %.4f '%Kc3t[a,ithp[p],b])
            fKc3_tex.write('\\\\\n')
            fKc3t_tex.write('\\\\\n')
        
        fKc3_tex.write('\\end{tabular}\n')
        fKc3t_tex.write('\\end{tabular}\n')

        #close files
        fK6at.close()
        fKc3at.close()
        fKc3tat.close()
        fKc3_tex.close()
        fKc3t_tex.close()

    fsumm.close()    
    
    #Show plots
    if plotarr.any():
        print 'Displaying plots...'
        plot.show()
    
    
#===============================================================================
#     Plotting routines
#===============================================================================
        
def plotrsrf(nuarr,rsrfarr,rsrfearr,apfarr,fsuff):
    #########################################################
    nuplt=nuarr/1.e9
    #Plot rsrf and aperture function
    plot.figure(1)
    plot.clf()
    #plot rsrf
    plot.plot(nuplt,rsrfarr[:,0],'b-',label='RSRF (PSW)')
    plot.plot(nuplt,rsrfarr[:,1],'g-',label='RSRF (PMW)')
    plot.plot(nuplt,rsrfarr[:,2],'r-',label='RSRF (PLW)')
    
    #plot nu^2.RSRF
    plot.plot(nuplt,rsrfearr[:,0],'b--',label=r'RSRF x $\nu^{-2}$')
    plot.plot(nuplt,rsrfearr[:,1],'g--',label='_nolegend_')
    plot.plot(nuplt,rsrfearr[:,2],'r--',label='_nolegend_')
    
    #plot aperture efficieny
    plot.plot(nuplt,apfarr[:,0],'b:',label='Ap. Eff.')
    plot.plot(nuplt,apfarr[:,1],'g:',label='_nolegend_')
    plot.plot(nuplt,apfarr[:,2],'r:',label='_nolegend_')
    
    plot.ylim(-0.1,1.1)
    plot.xlabel('Frequency [GHz]')
    plot.ylabel('RSRF ; Aperture Efficiency')
    #plot.title('RSRF and Aperture Efficiency')
    plot.legend()
    
    plot.savefig('../Outputs/rsrf'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/rsrf'+fsuff+'.eps',transparent=False)
    
    ########################################################

def plotbeam(beamx,beamy,beams,brad,fsuff,plimin=None):
    
    from numpy import arange,cos,sin,pi
    
    print 'Plotting beams...'
    print 'WARNING: this can be slow'
    #plot beams
    plot.figure(2)
    plot.clf()
    plot.hot()

    titles=array(['PSM Beam','PMW Beam','PLW Beam'])

    #make circle array
    tharr=arange(361.)/(2.*pi)
    xrad=brad*cos(tharr)
    yrad=brad*sin(tharr)
    
    for b in arange(3):
        print 'plot band %d'%b
        if b == 2:
            plot.subplot(2,1,2)
        else:
            plot.subplot(2,2,b+1)
        plot.axis('equal')
        bplot=log10(beams[:,:,b])
        #bplot=where(bplot == -pylab.inf,pylab.NaN,bplot) #replace -inf with NaN
        #bplot=where(bplot == pylab.inf,pylab.NaN,bplot) #replace inf with NaN
        #minb=min(bplot[where(bplot > -1.e30 )])
        #maxb=max(bplot[where(bplot < 1.e30 )])
        #print 'minb,maxb',minb,maxb
        print 'plimin',plimin
        if plimin == None:
            plim=array([-8,0])
        else: plim=array([plimin,0])
        print 'plim',plim
        (plim[0],plim[1])=(int(plim[0]),int(plim[1]))
        #bplot=where(bplot != bplot,plim[0],bplot) #replace NaNs with min
        #bplot=where(bplot < plim[0],plim[0],bplot) #replace min vals with min
        #bplot=where(bplot > plim[1],plim[1],bplot) #replace max vals with max
        prange=plim[1]-plim[0]
        #ntick=8
        ncol=32
        #tlevs=arange(plim[0],plim[1]+prange/ntick,prange/ntick)
        clevs=arange(plim[0],plim[1]+prange/ncol,prange/ncol)
        #print 'Ticks:',tlevs
        print 'Colours:',clevs
        plot.contourf(beamx,beamy,bplot,levels=clevs)
        plot.plot(xrad,yrad,'k:')
        #plot.plot(xrad,yrad,'w:')
        plot.clim(plim[0],plim[1])
        plot.xlim(min(beamx),max(beamx))
        plot.ylim(min(beamy),max(beamy))
        plot.title(titles[b])

        if b == 0:
            plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))),('','',''))
        else:
            plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))))
            plot.xlabel('x ["]')
        if b == 1:
            plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:]))),('','',''),rotation=90)
        else:
            plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:]))),rotation=90)
            plot.ylabel('y ["]')
        cb=plot.colorbar(format='%.1f')#,ticks=tlevs)
        cb.set_label('Log scale')
        
#    print 'plot PMW'
#    plot.subplot(2,2,2)
#    plot.axis('equal')
#    plot.contourf(beamx,beamy,bplot[:,:,1]-max(bplot[:,:,1]),levels=levs)
#    plot.clim(plim[0],plim[1])
#    plot.title('PMW Beam')
#    plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))))
#    #plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:])))),rotation=90)
#    plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:]))),('','',''),rotation=90)
#    cb=plot.colorbar(format='%.1f',ticks=levs)
#    cb.set_label('Log scale')
#    plot.xlabel('x ["]')
#        
#    print 'plot PLW'
#    plot.subplot(2,2,3)
#    plot.axis('equal')
#    plot.contourf(beamx,beamy,bplot[:,:,2]-max(bplot[:,:,2]),levels=levs)
#    plot.clim(plim[0],plim[1])
#    plot.title('PLW Beam')
#    plot.xticks((floor((min(beamx[:,0]))),0,ceil(max(beamx[:,0]))))
#    plot.yticks((floor((min(beamy[0,:]))),0,ceil(max(beamy[0,:]))),rotation=90)
#    #plot.colorbar()
#    cb=plot.colorbar(format='%.1f',ticks=levs)
#    cb.set_label('Log scale')
#    #colorbar.set_label('Log scale')
#    plot.xlabel('x ["]')
#    plot.ylabel('y ["]')
        
    plot.savefig('../Outputs/beams'+fsuff+'.png',transparent=False,dpi=300.)
    #plot.savefig('../Outputs/beams'+fsuff+'.eps',transparent=False)
    
def plotKc(alpharr,Kc,Kce,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(3)
    plot.clf()
    plot.plot(alpharr,Kc[:,0],'b-',label=r'$K_\mathrm{c,pt}$ (PSW)')
    plot.plot(alpharr,Kc[:,1],'g-',label=r'$K_\mathrm{c,pt}$ (PMW)')
    plot.plot(alpharr,Kc[:,2],'r-',label=r'$K_\mathrm{c,pt}$ (PLW)')
    plot.plot(alpharr,Kce[:,0],'b--',label=r'$K_\mathrm{c,ext}$ (PSW)')
    plot.plot(alpharr,Kce[:,1],'g--',label=r'$K_\mathrm{c,ext}$ (PMW)')
    plot.plot(alpharr,Kce[:,2],'r--',label=r'$K_\mathrm{c,ext}$ (PLW)')
        
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel('Colour Correction factor')
    #plot.title('Colour correction factors')
    plot.legend(ncol=2,loc='lower left')
    
    plot.savefig('../Outputs/colcorr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/colcorr'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotbeam_cut(cutx,cuty,beam_cutx,beam_cuty,brad,fsuff):
    ##plot beam cuts
    plot.figure(4)
    plot.clf()
    plot.subplot(3,1,1)
    plot.plot(cutx,beam_cutx[:,0,0],'g-',label='x (mid)')
    plot.plot(cutx,beam_cutx[:,1,0],'r-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,2,0],'b-',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,0,0],'g:',label='y (mid)')
    plot.plot(cuty,beam_cuty[:,1,0],'r:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,2,0],'b:',label='_nolegend_')
    plot.yscale('log')
    plot.ylabel('PSW Beam [Log]')
    plot.xlim(0,brad)
    plot.ylim(1.e-8,1.)
    plot.legend(loc='upper right',ncol=2)
    
    plot.subplot(3,1,2)
    plot.plot(cutx,beam_cutx[:,0,1],'g-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,1,1],'r-',label='x (low)')
    plot.plot(cutx,beam_cutx[:,2,1],'b-',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,0,1],'g:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,1,1],'r:',label='y (low)')
    plot.plot(cuty,beam_cuty[:,2,1],'b:',label='_nolegend_')
    plot.yscale('log')
    plot.ylabel('PMW Beam (Log)')
    plot.xlim(0,brad)
    plot.ylim(1.e-8,1.)
    plot.legend(loc='upper right',ncol=2)
    
    plot.subplot(3,1,3)
    plot.plot(cutx,beam_cutx[:,0,2],'g-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,1,2],'r-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,2,2],'b-',label='x (high)')
    plot.plot(cuty,beam_cuty[:,0,2],'g:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,1,2],'r:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,2,2],'b:',label='y (high)')
    plot.yscale('log')
    plot.ylabel('PLW Beam [Log]')
    plot.xlabel('Arcsec')
    plot.xlim(0,brad)
    plot.ylim(1.e-8,1.)
    plot.legend(loc='upper right',ncol=2)
    
    plot.savefig('../Outputs/beamcuts'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamcuts'+fsuff+'.eps',transparent=False)
    
def plotareaf(nuarr,area_bm,ilim,avgnu,avgarea,ith,fsuff,unit='sr'):
    
    from numpy import array
    from beam import arcsec2sr    
    ############################################################
    #Plot beam area vs frequency

    c=299792458. #speed of light
    wl0=array([250.,350.,500.])*1.e-6
    nu0=(c/wl0)/1.e9 #in GHz
    area0=array([426.,771.,1626.])

    avgnu=avgnu/1.e9

    if unit=='sr':
        #convert to steradians
        area0=arcsec2sr(area0)*1.e8
        area_bm=arcsec2sr(area_bm)*1.e8
        avgarea=arcsec2sr(avgarea)*1.e8
        ytit='Area [sr]'
    else: ytit=r'Area [arcsec$^2$]'
    
    plot.figure(5)
    plot.clf()
    
    plot.plot(nuarr[ilim[0,0]:ilim[0,1]]/1.e9,area_bm[ilim[0,0]:ilim[0,1],ith,0],'b-',label='PSW')
    plot.plot(nuarr[ilim[1,0]:ilim[1,1]]/1.e9,area_bm[ilim[1,0]:ilim[1,1],ith,1],'g-',label='PMW')
    plot.plot(nuarr[ilim[2,0]:ilim[2,1]]/1.e9,area_bm[ilim[2,0]:ilim[2,1],ith,2],'r-',label='PLW')
    plot.plot(nu0,area0,markeredgecolor='k',marker='^',linestyle='None',label='Current beam areas')
    plot.plot(avgnu[ith,0],avgarea[ith,0],markeredgecolor='k',marker='o',linestyle='None',label=r'$\nu$-averaged beam area')
    
    plot.plot(avgnu[ith,0],avgarea[ith,0],'bo')
    plot.plot(array([0,avgnu[ith,0]]),array([avgarea[ith,0],avgarea[ith,0]]),'b--')
    plot.plot(array([avgnu[ith,0],avgnu[ith,0]]),array([0,avgarea[ith,0]]),'b--')
    
    plot.plot(avgnu[ith,1],avgarea[ith,1],'go')
    plot.plot(array([0,avgnu[ith,1]]),array([avgarea[ith,1],avgarea[ith,1]]),'g--')
    plot.plot(array([avgnu[ith,1],avgnu[ith,1]]),array([0,avgarea[ith,1]]),'g--')
    
    plot.plot(avgnu[ith,2],avgarea[ith,2],'ro')
    plot.plot(array([0,avgnu[ith,2]]),array([avgarea[ith,2],avgarea[ith,2]]),'r--')
    plot.plot(array([avgnu[ith,2],avgnu[ith,2]]),array([0,avgarea[ith,2]]),'r--')

    #plot.plot(nu0[0],area0[0],'b^')
    #plot.plot(nu0[1],area0[1],'g^')
    #plot.plot(nu0[2],area0[2],'r^')
    #plot.title('Beam areas')
    plot.ylabel(ytit)
    plot.xlabel('Frequency [GHz]')
    plot.legend(loc='upper right',numpoints=1)
    plot.savefig('../Outputs/beamarea_nu'+fsuff+'_'+unit+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamarea_nu'+fsuff+'_'+unit+'.eps',transparent=False)

def plotarear(radarr,areas,fsuff,unit='sr'):
    
    from numpy import array
    from beam import arcsec2sr    
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    area0=array([426.,771.,1626.])

    if unit=='sr':
        #convert to steradians
        area0=arcsec2sr(area0)*1.e8
        areas=arcsec2sr(areas)*1.e8
        ytit='Area [sr]'
    else: ytit=r'Area [sq-arcsec$^2$]'
      
    plot.figure(6)
    plot.clf()
    plot.subplot(2,1,1)    
    plot.plot(radarr,areas[:,0],'b-',label='PSW')
    plot.plot(radarr,areas[:,1],'g-',label='PMW')
    plot.plot(radarr,areas[:,2],'r-',label='PLW')
    plot.axhline(y=area0[0],color='k',linestyle='--',label='Current beam area')    
    plot.axhline(y=area0[0],color='b',linestyle='--')
    plot.axhline(y=area0[1],color='g',linestyle='--')
    plot.axhline(y=area0[2],color='r',linestyle='--')
    #plot.title('Beam areas')
    plot.ylabel(ytit)
    plot.xlabel('Radius ["]')
    #plot.legend(loc='lower right',numpoints=1,ncol=2)
    
    na=areas[:,0].size
    plot.subplot(2,1,2)
    plot.plot(radarr,(areas[:,0]-areas[na-1,0])/areas[na-1,0],'b-',label='PSW')
    plot.plot(radarr,(areas[:,1]-areas[na-1,1])/areas[na-1,1],'g-',label='PMW')
    plot.plot(radarr,(areas[:,2]-areas[na-1,2])/areas[na-1,2],'r-',label='PLW')
    
    plot.axhline((area0[0]-areas[na-1,0])/areas[na-1,0],color='k',linestyle='--',label='Current beam area')
    plot.axhline((area0[0]-areas[na-1,0])/areas[na-1,0],color='r',linestyle='--')
    plot.axhline((area0[1]-areas[na-1,1])/areas[na-1,1],color='g',linestyle='--')
    plot.axhline((area0[2]-areas[na-1,2])/areas[na-1,2],color='b',linestyle='--')
    plot.axhline(0,color='k',linestyle='-')
    plot.xlabel('Radius ["]')
    plot.ylabel(r'$\frac{Area(r)}{Area(r_\mathrm{max})} - 1$')
    plot.legend(loc='upper right',ncol=2)
    plot.ylim(-0.05,0.05)

    plot.savefig('../Outputs/beamarea_rad'+fsuff+'_'+unit+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamarea_rad'+fsuff+'_'+unit+'.eps',transparent=False)
    

def plotbeam_sm(radarr,beam_sm,beam_min,beam_max,beam_sd,brad,fsuff):

    import matplotlib.pyplot as plot

    ##############################################################
    ##Plot az-averaged beams
    ##############################################################
    plot.figure(7)
    plot.clf()
    plot.subplot(3,1,1)
    plot.plot(radarr,beam_sm[:,0],'b-',label='Mean')
    plot.plot(radarr,beam_min[:,0],'b--',label='Min/Max')
    plot.plot(radarr,beam_max[:,0],'b--',label='_nolegend_')
    plot.plot(radarr,(beam_sm+beam_sd)[:,0],'b:',label='Std. Dev.')
    plot.plot(radarr,(beam_sm-beam_sd)[:,0],'b:',label='_nolegend_')
    #plot.title('Azimuthally averaged beams')
    plot.ylabel('PSW')
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    #plot.xscale('log')
    plot.legend(loc='upper right',ncol=3)
    
    plot.subplot(3,1,2)
    plot.plot(radarr,beam_sm[:,1],'g-',label='Mean')
    plot.plot(radarr,beam_min[:,1],'g--',label='Min/Max')
    plot.plot(radarr,beam_max[:,1],'g--',label='_nolegend_')
    plot.plot(radarr,(beam_sm+beam_sd)[:,1],'g:',label='Std. Dev.')
    plot.plot(radarr,(beam_sm-beam_sd)[:,1],'g:',label='_nolegend_')
    plot.ylabel('PMW')
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    #plot.xscale('log')
    
    plot.subplot(3,1,3)
    plot.plot(radarr,beam_sm[:,2],'r-',label='Mean')
    plot.plot(radarr,beam_min[:,2],'r--',label='Min/Max')
    plot.plot(radarr,beam_max[:,2],'r--',label='_nolegend_')
    plot.plot(radarr,(beam_sm+beam_sd)[:,2],'r:',label='Std. Dev.')
    plot.plot(radarr,(beam_sm-beam_sd)[:,2],'r:',label='_nolegend_')
    plot.ylabel('PLW')
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    #plot.xscale('log')
    plot.xlabel('Radius [arcsec]')
    
    plot.savefig('../Outputs/azbeam_range'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/azbeam_range'+fsuff+'.eps',transparent=False)
    
    ###################################################
    ##plot all three az-averaged beams
    plot.figure(8)
    plot.clf()
    plot.plot(radarr,beam_sm[:,0],'b-',label='PSW')
    plot.plot(radarr,beam_sm[:,1],'g-',label='PMW')
    plot.plot(radarr,beam_sm[:,2],'r-',label='PLW')
    #plot.plot(radarr,beam_max[:,0],'b--',label='Max val')
    #plot.plot(radarr,beam_max[:,1],'g--',label='_nolegend_')
    #plot.plot(radarr,beam_max[:,2],'r--',label='_nolegend_')
    
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    plot.ylabel('Beam profile')
    plot.xlabel('Radius [arcsec]')
    plot.legend(loc='upper right')
    
    plot.savefig('../Outputs/azbeam_all'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/azbeam_all'+fsuff+'.eps',transparent=False)
    
def plotK4(alpharr,K4,K4e,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(9)
    plot.clf()
    plot.plot(alpharr,K4[:,0],'b-',label=r'$K_\mathrm{4,pt}$ (PSW)')
    plot.plot(alpharr,K4[:,1],'g-',label=r'$K_\mathrm{4,pt}$ (PMW)')
    plot.plot(alpharr,K4[:,2],'r-',label=r'$K_\mathrm{4,pt}$ (PLW)')
    plot.plot(alpharr,K4e[:,0],'b--',label=r'$K_\mathrm{4,ext}$ (PSW)')
    plot.plot(alpharr,K4e[:,1],'g--',label=r'$K_\mathrm{4,ext}$ (PMW)')
    plot.plot(alpharr,K4e[:,2],'r--',label=r'$K_\mathrm{4,ext}$ (PLW)')
    
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'K$_4$ Correction factor')
    #plot.title('Current correction factors')
    plot.legend(ncol=2,loc='lower left')
    
    plot.savefig('../Outputs/K4corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/K4corr'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotK6(arr_a,arr_t,K6,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(10)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(K6[:,:,b])),max(log10(K6[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(K6[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1: plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Spectral index ($\alpha$)')
        plot.colorbar(format='%.1f')
    
    plot.savefig('../Outputs/K5corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/K5corr'+fsuff+'.eps',transparent=False)
    
    #############################################################

def plotKc3(arr_a,arr_t,Kc3,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(11)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(Kc3[:,:,b])),max(log10(Kc3[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(Kc3[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1:plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Spectral index ($\alpha$)')
        plot.colorbar(format='%.1f')
    
    plot.savefig('../Outputs/Kc2corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/Kc2corr'+fsuff+'.eps',transparent=False)
    
    #############################################################

def plotKc3t(arr_a,arr_t,Kc3t,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(12)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(Kc3t[:,:,b])),max(log10(Kc3t[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(Kc3t[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1:plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Spectral index ($\alpha$)')
        plot.colorbar()
    
    plot.savefig('../Outputs/Kc3corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/Kc3corr'+fsuff+'.eps',transparent=False)
    
    #############################################################


def plotareaa(alpharr,areaeff,areas,fsuff,unit='sr'):
    
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
    plot.savefig('../Outputs/beamarea_alpha'+fsuff+'_'+unit+'.eps',transparent=False)
    
def plotbeamcomp(radarr,beam_sm,beam_mod,area_rarr,area_rel,fsuff):
    
    fwhm=array((17.6,23.9,35.1))    
    plot.figure(14)
    plot.clf()

    xticksep=20
    xtickv=arange(0,100+xticksep,xticksep)
    nt=xtickv.size
    xtickl=[]
    xtickbl=[]
    for n in arange(nt): 
        xtickl.append(str(xtickv[n]))
        xtickbl.append('')    
    
    plot.subplot(3,1,1)
    legnep,=plot.plot(radarr,beam_sm[:,0],'k-',label=r'$B_\mathrm{Nep}$')
    sm0=plot.plot(radarr,beam_sm[:,0],'b-')
    sm1=plot.plot(radarr,beam_sm[:,1],'g-')
    sm2=plot.plot(radarr,beam_sm[:,2],'r-')
    legmod,=plot.plot(radarr,beam_mod[:,0],'k--',label=r'$B_\mathrm{mod}$')
    plot.plot(radarr,beam_mod[:,0],'b--')
    plot.plot(radarr,beam_mod[:,1],'g--')
    plot.plot(radarr,beam_mod[:,2],'r--')
    plot.ylabel('Beam response')
    plot.yscale('log')
    #plot.xlabel('Radius ["]')
    plot.ylim(1e-4,1.)
    plot.xlim(0.,100.)
    plot.xticks(xtickv,xtickbl)
    plot.axvline(x=fwhm[0]/2,color='b',linestyle=':')
    plot.axvline(x=fwhm[1]/2,color='g',linestyle=':')
    plot.axvline(x=fwhm[2]/2,color='r',linestyle=':')
    l1=plot.legend((legnep,legmod),(r'$B_\mathrm{Nep}$',r'$B_\mathrm{mod}$'),loc='upper right',ncol=2)
    plot.legend((sm0,sm1,sm2),('PSW','PMW','PLW'),loc='lower left',ncol=3)
    plot.gca().lines.remove(legnep)
    plot.gca().lines.remove(legmod)
    plot.gca().add_artist(l1)
    
    
    plot.subplot(3,1,2)
    plot.plot(radarr,beam_mod[:,0]-beam_sm[:,0],'b-',label='PSW')
    plot.plot(radarr,beam_mod[:,1]-beam_sm[:,1],'g-',label='PWM')
    plot.plot(radarr,beam_mod[:,2]-beam_sm[:,2],'r-',label='PLW')
    plot.ylim(-0.02,0.02)
    plot.xlim(0,100)
    plot.axvline(x=fwhm[0]/2,color='b',linestyle=':')
    plot.axvline(x=fwhm[1]/2,color='g',linestyle=':')
    plot.axvline(x=fwhm[2]/2,color='r',linestyle=':')
    plot.axhline(y=0.,color='k',linestyle='-')
    plot.ylabel(r'$B_\mathrm{mod} - B_\mathrm{Nep}$')
    plot.xticks(xtickv,xtickbl)
    plot.legend(loc='lower right',ncol=3)

    #print 'min/max (PSW):',min((beam_mod[:,0]-beam_sm[:,0])/beam_sm[:,0]),max((beam_mod[:,0]-beam_sm[:,0])/beam_sm[:,0])
    #print 'min/max (PMW):',min((beam_mod[:,1]-beam_sm[:,1])/beam_sm[:,1]),max((beam_mod[:,1]-beam_sm[:,1])/beam_sm[:,1])
    #print 'min/max (PLW):',min((beam_mod[:,2]-beam_sm[:,2])/beam_sm[:,2]),max((beam_mod[:,2]-beam_sm[:,2])/beam_sm[:,2])
    plot.subplot(3,1,3)
    plot.plot(area_rarr,area_rel[:,0],'b-',label='PSW')
    plot.plot(area_rarr,area_rel[:,1],'g-',label='PMW')
    plot.plot(area_rarr,area_rel[:,2],'r-',label='PLW')
    plot.axhline(y=0,color='k',linestyle='-')
    #plot.yscale('symlog')#,basey=10.,linthresh=0.01,subsy=[0,1,2,3,4,5,6,7,8,9])
    #plot.ylim(-10.,10.)

    plot.axvline(x=fwhm[0]/2,color='b',linestyle=':')
    plot.axvline(x=fwhm[1]/2,color='g',linestyle=':')
    plot.axvline(x=fwhm[2]/2,color='r',linestyle=':')    

    plot.xlabel(r'Radius, $\theta_r$ ["]')
    #plot.ylabel(r'$(\Omega_\mathrm{eff}\big|_{\theta<\theta_r}-\Omega_\mathrm{Nep}\big|_{\theta<\theta_r})$')
    #plot.ylabel(r'$(\Omega_\mathrm{eff}(\theta_r)-\Omega_\mathrm{Nep}(\theta_r)/\Omega_\mathrm{Nep}(\theta_r)$')
    plot.ylabel(r'$\left[(\Omega_\mathrm{eff}-\Omega_\mathrm{Nep})/\Omega_\mathrm{Nep}\right]_{\theta<\theta_r}$')
    plot.xlim(0.,100.)
    plot.xticks(xtickv,xtickl)
    #plot.yticks((-0.02,0,0.02))
    plot.legend(loc='lower right',ncol=3)
    
    plot.savefig('../Outputs/beamcomp'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamcomp'+fsuff+'.eps',transparent=False)

def plotareath(thetarr,areas,fsuff,unit='sr'):
    
    from numpy import array
    from beam import arcsec2sr    
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    area0=array([426.,771.,1626.])

    if unit=='sr':
        #convert to steradians
        area0=arcsec2sr(area0)*1.e8
        areas=arcsec2sr(areas)*1.e8
        ytit=r'Area [sr] $\times 10^8$'
    else: ytit=r'Area [sq-arcsec$^2$]'
      
    plot.figure(15)
    plot.clf()
    plot.plot(thetarr,areas[0,:],'b-',label='PSW')
    plot.plot(thetarr,areas[1,:],'g-',label='PMW')
    plot.plot(thetarr,areas[2,:],'r-',label='PLW')
    plot.axhline(y=area0[0],color='k',linestyle='--',label='Current beam area')    
    plot.axhline(y=area0[0],color='b',linestyle='--')
    plot.axhline(y=area0[1],color='g',linestyle='--')
    plot.axhline(y=area0[2],color='r',linestyle='--')
    #plot.title('Beam areas')
    plot.ylabel(ytit)
    plot.xlabel('Source FHWM ["]')
    plot.legend(loc='lower right',numpoints=1,ncol=2)
    
    plot.savefig('../Outputs/beamarea_theta'+fsuff+'_'+unit+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamarea_theta'+fsuff+'_'+unit+'.eps',transparent=False)

def plotK6_va(alpharr,thetarr,K6,K5,aplot,fsuff):

    plot.figure(16)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    for b in range(3):
        plot.plot(thetarr,K6[aplot,:,b]/1.e6,color=cols[b],label=r'$K_6$ (%s)'%(tit[b]))
    for b in range(3)        :
        plot.axhline(K5[aplot,b]/1.e6,linestyle='--',color=cols[b],label=r'$K_5$ (%s)'%(tit[b]))
    plot.xlabel(r'Source FWHM, $\theta_0$ ["]')
    plot.ylabel(r'Conversion $K_6$ (MJy/sr per Jy measured)')
    plot.legend(loc='upper right',ncol=2)
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(1.e1,1.e7)
    plot.xlim(1.e-1,1.e4)
    plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))
    
    filename='../Outputs/K6_a%.1f'%(alpharr[aplot])
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False,dpi=300)

def plotKc3_va(alpharr,thetarr,Kc3,Kc2,aplot,fsuff):

    plot.figure(17)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    for b in range(3):
        plot.plot(thetarr,Kc3[aplot,:,b]/1.e6,color=cols[b],label=r'$K_\mathrm{c3}$ (%s)'%(tit[b]))
    for b in range(3):        
        plot.axhline(Kc2[aplot,b]/1.e6,linestyle='--',color=cols[b],label=r'$K_\mathrm{c2}$ (%s)'%(tit[b]))
    plot.xlabel(r'Source FWHM, $\theta_0$ ["]')
    plot.ylabel(r'Correction $K_\mathrm{c3}$ (MJy/sr per Jy from pipeline)')
    plot.legend(loc='upper right',ncol=2)
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(1.e1,1.e7)
    plot.xlim(1.e-1,1.e4)
    plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))
    
    filename='../Outputs/Kc3_a%.1f'%(alpharr[aplot])
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False,dpi=300)

def plotKc3t_va(alpharr,thetarr,Kc3t,Kc,aplot,fsuff):

    plot.figure(18)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']    
    for b in range(3):
        plot.plot(thetarr,Kc3t[aplot,:,b],color=cols[b],label=r'$K_\mathrm{c3}^\mathrm{tot}$ (%s)'%(tit[b]))
    for b in range(3):        
        plot.axhline(Kc[aplot,b],linestyle='--',color=cols[b],label=r'$K_\mathrm{c,pt}$ (%s)'%(tit[b]))
    plot.xlabel('Source FWHM ["]')
    plot.ylabel(r'Correction $K_\mathrm{c3}^\mathrm{tot}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.legend(loc='upper left',ncol=2)
    plot.yscale('linear')
    plot.xscale('log')
    plot.ylim(0,10)
    plot.xlim(1,100)
    plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))    
    
    filename='../Outputs/Kc3t_a%.1f'%(alpharr[aplot])
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False)

def plotKc_all(alpharr,Kc,Kce,Kc3,Kc3t,fsuff):
   
   plot.figure(19)
   plot.clf()
   tit=['PSW','PMW','PLW']
   cols=['b','g','r']
   legKc,=plot.plot(alpharr,Kc[:,0],color='k',linestyle=':')
   legKce,=plot.plot(alpharr,Kce[:,0],color='k',linestyle='--')
   legKc3t,=plot.plot(alpharr,Kc3t[:,0],color='k',linestyle='-')
   legpsw,=plot.plot(alpharr,Kc3t[:,0],color='b',linestyle='-')
   legpmw,=plot.plot(alpharr,Kc3t[:,0],color='g',linestyle='-')
   legplw,=plot.plot(alpharr,Kc3t[:,0],color='r',linestyle='-')
   for b in range(3):
       plot.plot(alpharr,Kc[:,b],color=cols[b],linestyle=':')
       plot.plot(alpharr,Kce[:,b],color=cols[b],linestyle='--')
       #plot.plot(alpharr,Kc2[:,b],color=cols[b],linestyle='-')
       plot.plot(alpharr,Kc3t[:,b],color=cols[b],linestyle='-')
   plot.xlabel(r'Spectral index $\alpha$')
   plot.ylabel('Correction factor')
   plot.gca
   plot.legend((legKc,legKce,legKc3t,legpsw,legpmw,legplw),
               (r'$K_\mathrm{c,pt}$',r'$K_\mathrm{c,ext}$',r'$K_\mathrm{c3}^\mathrm{tot}$','PSW','PMW','PLW'),
                loc='lower left',ncol=2)
   plot.gca().lines.remove(legKc)
   plot.gca().lines.remove(legKce)
   plot.gca().lines.remove(legKc3t)
   plot.gca().lines.remove(legpsw)
   plot.gca().lines.remove(legpmw)
   plot.gca().lines.remove(legplw)

   plot.savefig('../Outputs/Kc_all'+fsuff+'.png',transparent=False,dpi=300)
   plot.savefig('../Outputs/Kc_all'+fsuff+'.eps',transparent=False)
     
def plotKcKc3t_thsmall(alpharr,Kc,Kc3t,thplot,fsuff):

    plot.figure(20)
    plot.clf()

    tit=['PSW','PMW','PLW']
    cols=['b','g','r']
    
    for b in range(3):
        plot.plot(alpharr,1.e3*(Kc[:,b]-Kc3t[:,b])/Kc[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$10^3 \times \ (K_\mathrm{c}-K_\mathrm{c3}^\mathrm{tot})/K_\mathrm{c}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Spectral index $\alpha$')
    plot.legend(loc='upper left',ncol=3)
    plot.ylim(-4,1)
    
    plot.savefig('../Outputs/KcKc3t_thsmall'+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/KcKc3t_thsmall'+fsuff+'.eps',transparent=False)
    
def plotKc2Kc3_thbig(alpharr,Kc3,Kc2,thplot,fsuff):

    plot.figure(21)
    plot.clf()
    tit=['PSW','PMW','PLW']
    cols=['b','g','r']
    for b in range(3):
        plot.plot(alpharr,1.e5*(Kc2[:,b]-Kc3[:,b])/Kc2[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$10^5 \times \ (K_\mathrm{c2}-K_\mathrm{c3})/K_\mathrm{c2}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Spectral index $\alpha$')

    plot.legend(loc='upper right',ncol=3)
    plot.savefig('../Outputs/Kc2Kc3_thbig'+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/Kc2Kc3_thbig'+fsuff+'.eps',transparent=False)

def plotKc2(alpharr,Kce,Kc2,area0,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(22)
    plot.clf()
    plot.subplot(2,1,1)
    
    lkc2,=plot.plot(alpharr,Kc2[:,0],'k-',label=r'$K_\mathrm{c2}$')
    lkce,=plot.plot(alpharr,Kce[:,0]/area0[0],'k--',label=r'$K_\mathrm{c,ext}/\Omega_0$')
    lpsw,=plot.plot(alpharr,Kc2[:,0],'b-',label=r'$K_\mathrm{c2}$ (PSW)')
    lpmw,=plot.plot(alpharr,Kc2[:,1],'g-',label=r'$K_\mathrm{c2}$ (PMW)')
    lplw,=plot.plot(alpharr,Kc2[:,2],'r-',label=r'$K_\mathrm{c2}$ (PLW)')
        
    plot.plot(alpharr,Kce[:,0]/area0[0],'b--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PSW)')
    plot.plot(alpharr,Kce[:,1]/area0[1],'g--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PMW)')
    plot.plot(alpharr,Kce[:,2]/area0[2],'r--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PLW)')
    
    #plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'Conversion factor (Jy/sr per Jy$_\mathrm{pip}$)')
    plot.ylim(ymax=1.3e8)
    #plot.title('Colour correction factors')        
    #plot.legend(ncol=2,loc='lower right')
    plot.legend((lkc2,lkce),(r'$K_\mathrm{c2}$',r'$K_\mathrm{c,ext}/\Omega_0$'),ncol=2,loc='upper right')

    plot.subplot(2,1,2)
    plot.plot(alpharr,(Kce[:,0]/area0[0])/Kc2[:,0],'b-',label='PSW')
    plot.plot(alpharr,(Kce[:,1]/area0[1])/Kc2[:,1],'g-',label='PMW')
    plot.plot(alpharr,(Kce[:,2]/area0[2])/Kc2[:,2],'r-',label='PLW')
    
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'$(K_\mathrm{c,ext}/\Omega_0)\, / \, K_\mathrm{c2}$')
    plot.ylim(1.0,1.12)
    plot.legend(ncol=3,loc='lower right')
        
    plot.savefig('../Outputs/colcorr2'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/colcorr2'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotKc2rel(alpharr,Kce,Kc2,area0,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(23)
    plot.clf()
    
    #lpsw,=plot.plot(alpharr,Kce[:,0]/Kc2[:,0],'b-')
    #lpmw,=plot.plot(alpharr,Kce[:,1]/Kc2[:,1],'g-')
    #lplw,=plot.plot(alpharr,Kce[:,2]/Kc2[:,2],'r-')

    #lext,=plot.plot(alpharr,Kce[:,0]/Kc2[:,0],'k:',label=r'$K_\mathrm{c,ext}/K_\mathrm{c2}^\mathrm{bm}$')    
    plot.plot(alpharr,(Kce[:,0]/area0[0])/Kc2[:,0],'b-',label='PSW')
    plot.plot(alpharr,(Kce[:,1]/area0[1])/Kc2[:,1],'g-',label='PMW')
    plot.plot(alpharr,(Kce[:,2]/area0[2])/Kc2[:,2],'r-',label='PLW')
        
    #lpt,=plot.plot(alpharr,Kc[:,0]/Kc2[:,0],'k--',label=r'$K_\mathrm{c,pt}/K_\mathrm{c2}^\mathrm{bm}$')
    #plot.plot(alpharr,Kc[:,0]/Kc2[:,0],'b--')
    #plot.plot(alpharr,Kc[:,1]/Kc2[:,1],'g--')
    #plot.plot(alpharr,Kc[:,2]/Kc2[:,2],'r--')
    
    #lpip,=plot.plot(alpharr,Kpip[0]/Kc2[:,0],'k-.',label=r'$K_\mathrm{pip}/K_\mathrm{c2}^\mathrm{bm}$')
    #plot.plot(alpharr,Kpip[0]/Kc2[:,0],'b-.')
    #plot.plot(alpharr,Kpip[1]/Kc2[:,1],'g-.')
    #plot.plot(alpharr,Kpip[2]/Kc2[:,2],'r-.')
    
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'$(K_\mathrm{c,ext}/\Omega_0)\, / \, K_\mathrm{c2}$')
    
    plot.legend(loc='lower left')
    #plot.title('Colour correction factors')
    #plot.legend(ncol=3,loc='upper right')
    #plot.legend((lext,lpt,lpsw,lpmw,lplw),
    #            (r'$K_\mathrm{c,ext}/K_\mathrm{c2}$',r'$K_\mathrm{c,pt}/K_\mathrm{c2}$',
    #             'PSW','PMW','PLW'),
    #             ncol=3,loc='upper right')
    
    #plot.gca().lines.remove(lext)
    #plot.gca().lines.remove(lpt)
    #plot.gca().lines.remove(lpsw)
    #plot.gca().lines.remove(lpmw)
    #plot.gca().lines.remove(lplw)
    
    plot.savefig('../Outputs/colcorr2rel'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/colcorr2rel'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotblim(radarr,beam_scl_lim,beam_cmb_lim,beam_fix,brad,fsuff):

    import matplotlib.pyplot as plot

    ##############################################################
    ##Plot az-averaged beams
    ##############################################################
    plot.figure(24)
    plot.clf()
    tit=['PSW','PMW','PLW']
    for b in range(3):
        plot.subplot(3,1,b+1)
        plot.plot(radarr,beam_scl_lim[:,0,b],'r:')
        plot.plot(radarr,beam_scl_lim[:,1,b],'g:')
        plot.plot(radarr,beam_scl_lim[:,2,b],'b:')
        
        plot.plot(radarr,beam_fix[:,b],'k--')
        
        plot.plot(radarr,beam_cmb_lim[:,0,b],'r-',label=r'$\nu_\mathrm{lo}$')
        plot.plot(radarr,beam_cmb_lim[:,1,b],'g-',label=r'$\nu_\mathrm{c}$')
        plot.plot(radarr,beam_cmb_lim[:,2,b],'b-',label=r'$\nu_\mathrm{hi}$')
                
        plot.yscale('log')
        plot.ylim(1.e-8,1)
        plot.xlim(0,brad)
        if b == 0: plot.legend(ncol=3,loc='upper right')
        plot.ylabel(tit[b])
        if b==2: plot.xlabel('radius [arcsec]')
    
    plot.savefig('../Outputs/beamlim'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamlim'+fsuff+'.eps',transparent=False)
    
def plotbcomb(radarr,beam_scl,beam_fix,beam_cmb,bsw,fsuff):

    plot.figure(30)
    plot.clf()    
    plot.plot(radarr,beam_scl[:,0],'b:')
    plot.plot(radarr,beam_scl[:,1],'g:')
    plot.plot(radarr,beam_scl[:,2],'r:')
    plot.plot(radarr,beam_fix[:,0],'b--')
    plot.plot(radarr,beam_fix[:,1],'g--')
    plot.plot(radarr,beam_fix[:,2],'r--')
    lpsw,=plot.plot(radarr,beam_cmb[:,0],'b-',label='PSW')
    lcmb,=plot.plot(radarr,beam_cmb[:,0],'k-',label='Combined')
    lpmw,=plot.plot(radarr,beam_cmb[:,1],'g-',label='PMW')
    lfix,=plot.plot(radarr,beam_cmb[:,0],'k--',label='Fixed')
    lplw,=plot.plot(radarr,beam_cmb[:,2],'r-',label='PLW')
    lscl,=plot.plot(radarr,beam_cmb[:,0],'k:',label='Scaled')
    plot.plot(bsw[0,:],beam_cmb[bsw[0,:],0],'bx')
    plot.plot(bsw[1,:],beam_cmb[bsw[1,:],1],'gx')
    plot.plot(bsw[2,:],beam_cmb[bsw[2,:],2],'rx')
    plot.yscale('log')
    plot.legend(ncol=3,loc='upper right')
    
    plot.gca().lines.remove(lcmb)
    plot.gca().lines.remove(lscl)
    plot.gca().lines.remove(lfix)
    
    plot.savefig('../Outputs/beamcomb'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamcomb'+fsuff+'.eps',transparent=False)
    
if __name__ == "__main__":
    maincode()