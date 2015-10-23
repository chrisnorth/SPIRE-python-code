#from __future__ import division
#General python modules
#import pylab
from numpy import zeros,arange,array,floor,ceil,log,log10,sum,max,min,pi,sqrt
from numpy import concatenate,inf
#import sys
#import matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plot
from scipy import where
#from scipy.interpolate import RectBivariateSpline
import argparse

#Specific modules
from beam import getbeam,beam_regrid,beam_azsm,get_effnu_az,beamarea
from beam import arcsec2sr,avgbeam,beamarea_eff,modbeam,modbeam_area
from beam import beamarea_az_nu,beamarea_az_th_nu,beamarea_az_th
from rsrf import getrsrf,getapf,bandedge
from calc_conv import calc_k4,calc_kc,calc_k5

def maincode():
    parser=argparse.ArgumentParser(description='Herschel-SPIRE Extended Emission Calibration')
    parser.add_argument('--beam',action='store',default='M',dest='beam',help='Type of beam [M(easued)|T(heoretical)|G(aussian)|E(lliptical Gaussian)]. Default=M')
    parser.add_argument('--brad',action='store',default=350.,type=float,dest='brad',help='Radius to integrate beam out to (arcsec). Default=350')
    parser.add_argument('--bwid',action='store',default=2000.,type=float,dest='bwid',help='Full width of beam read from file (arcsec). Default=2000')
    parser.add_argument('--bzlim',action='store',default=None,type=float,dest='bzlim',help='Zero-limit for beamfiles, below which the values are replaced by BZVAL. Default=None (i.e. don\'t set limit)')
    parser.add_argument('--bplim',action='store',default=None,type=float,dest='bplim',help='Logarithm (base 10) of lower limit on beam plots. Default=-8')
    parser.add_argument('--bzval',action='store',default=0.,type=float,dest='bzval',help='Value with which to replace beam values below the zero-limit (BZLIM). Default=0')
    parser.add_argument('--aunit',action='store',default="sr",dest='aunit',help='Unit of area to use for area plots [arcsec|sr]. Default=sr')
    parser.add_argument('--bcut', action='store_true',default=False,dest='bcut',help='Set to plot beam cuts')
    parser.add_argument('--bazsm', action='store_true',default=False,dest='bazsm',help='Set to plot azimuth-smoothed beams')
    parser.add_argument('--newa', action='store_true',default=False,dest='newa',help='Set to use new area in colour-correction calculations')

    parser.add_argument('--psrsf', action='store_true',default=False,dest='prsrf',help='Set to plot RSRF and Ap. Eff.')
    parser.add_argument('--pbeam', action='store_true',default=False,dest='pbeam',help='Set to plot beam maps (overridden by NOPBEAM)')
    parser.add_argument('--pbcomp', action='store_true',default=False,dest='pbcomp',help='Set to plot beam comparison profiles (measured vs. theoretical)')
    parser.add_argument('--nopbeam', action='store_true',default=False,dest='nopbeam',help='Set to NOT plot beam maps (overrides PBEAM)')
    parser.add_argument('--pcorr', action='store_true',default=False,dest='pcorr',help='Set to plot colour correction factors')
    parser.add_argument('--pareaf', action='store_true',default=False,dest='pareaf',help='Set to plot beam area against frequency')
    parser.add_argument('--parear', action='store_true',default=False,dest='parear',help='Set to plot beam area against radius')
    parser.add_argument('--pareaa', action='store_true',default=False,dest='pareaa',help='Set to plot beam area against spectral index')
    parser.add_argument('--pareath', action='store_true',default=False,dest='pareath',help='Set to plot beam area against source FWHM')
    parser.add_argument('--pall', action='store_true',default=False,dest='pall',help='Plot all plots (set NOPBEAM to exclude beam maps')

    args=parser.parse_args()
    bwid=args.bwid
    bzlim=args.bzlim
    bzval=args.bzval
    brad=args.brad
    pall=args.pall
    aunit=args.aunit
    beamtype=args.beam
    newa=args.newa
    
    if args.bplim: bplim=args.bplim
    else: bplim=None
    if pall:
        bcut=True
        bazsm=True
        prsrf=True
        pbeam=True
        pbcomp=True
        pcorr=True
        pareaf=True
        parear=True
        pareaa=True
        pareath=True
    else:
        bcut=args.bcut
        bazsm=args.bazsm
        prsrf=args.prsrf
        pbeam=args.pbeam
        pcorr=args.pcorr
        pareaf=args.pareaf
        parear=args.parear
        pareaa=args.pareaa
        pbcomp=args.pbcomp
        pareath=args.pareath
    
    if args.nopbeam: pbeam=False

    fsuff='_theta_B'+beamtype+'_Bw%d_Br%d'%(int(bwid),int(brad))
    if bzlim != None:
        fsuff=fsuff+'_Bl%1g_Bv%1g'%(bzlim,bzval)
    if bplim != None:
        fsuff=fsuff+'_Bpl%g'%(bplim)
    if newa == True:
        fsuff=fsuff+'_newArea'
    print 'File_suffix= '+fsuff
    
    plotarr=array([prsrf,pbeam,pcorr,pareaf,bazsm,bcut,parear,pareaa,pbcomp,pareath])    

    mpl.rcParams['legend.fontsize']='medium'

    
    c=299792458. #speed of light
    #set band centres
    wlc=(250.e-6,350.e-6,500.e-6)
    wlc=array(wlc)
    nuc=c/wlc
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
    ##  Read in beam  ##
    ####################
    
    #Beamtype options:
    #  G=Gaussian
    #  E=Elliptical Gaussian
    #  M=Measured beams
    #  T=Theoretical (modelled) beams)
    #beamtype='M' #read from command line
    print 'Getting beams of type "%s"...' % (beamtype)
    
    #Set beam grid options
    bgrid=1. #grid size in arcsec
    #bwid=1000 #beam width in arcsec
    #bzlim=1.e-8
    #bzval=0.
    npx=bwid/bgrid
    beams=zeros((npx,npx,3))
    areas=zeros(3)
    for b in band:
        beamx,beamy,beams[:,:,b]=getbeam(beamtype,b,bgrid,bwid,bzlim=bzlim,bzval=bzval)
        areas[b]=sum(beams[:,:,b])
        
    print 'Beam areas: [%.2f , %.2f , %.2f] sq.arcsec' % (areas[0],areas[1],areas[2])
    
    if pbeam:
        plotbeam(beamx,beamy,beams,brad,fsuff,plimin=bplim)
    
    #Set effective frequencies for the beams
    
    ####################
    ##  Read in RSRF  ##
    ####################
    
    #RSRFtype options:
    #  T=Top Hat
    #  M=Measured
    rsrftype='M'
    print 'Getting RSRF of type "%s"...' % (rsrftype)
    rsrfarr=zeros((nnu,3))
    ilim=zeros((3,2))
    nulim=zeros((3,2))
    ilim_2=zeros((3,2))
    nulim_2=zeros((3,2))
    for b in band:
        rsrfarr[:,b]=getrsrf(rsrftype,b,nuarr)
        (ilim[b,:],nulim[b,:])=bandedge(nuarr,rsrfarr[:,b])
        (ilim_2[b,:],nulim_2[b,:])=bandedge(nuarr,rsrfarr[:,b],fact=2.)
        print 'Band %d limits: [%.2f:%.2f] GHz ([%.2f:%.2f] um)' % (b,nulim[b,0]/1.e9,nulim[b,1]/1.e9,c/nulim[b,1]*1.e6,c/nulim[b,0]*1.e6)
        #print 'Band %d limits: [%.2f:%.2f] um' % (b,c/nulim[b,1]*1.e6,c/nulim[b,0]*1.e6)
    #Aperture Efficiency options:
    # R=Real
    # U=Uniform
    apftype='R'
    print 'Getting Aperture Efficiency of type "%s"...' % (apftype)
    apfarr=zeros((nnu,3))
    for b in band:
        apfarr[:,b]=getapf(apftype,b,nuarr)
    
    #Compute 'extended source' rsrf
    rsrfearr=zeros((nnu,3))
    for b in band:
        rsrfearr[:,b]=rsrfarr[:,b]*(nuc[b]/nuarr)**2
    
    if prsrf:
        plotrsrf(nuarr,rsrfarr,rsrfearr,apfarr,fsuff)
    
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
    
    
    ######################
    ##  Make beam cuts  ##
    ######################
    ##Will calculate this empirically eventually
    
    if bcut:
        print 'Making beam cuts...'
        a_bm=array([1.0184,1.0164,1.0204])
        nu0_bm=nuc*a_bm
        nuL_bm=array([1033.,741.,491.])*1.e9
        nuU_bm=array([1418.,1008.,732.])*1.e9
        
        ncut=3
        beam_cutx=zeros((npx,ncut,3))
        beam_cuty=zeros((npx,ncut,3))
        p0=int(floor(npx/2.))
        cutx=beamx[:,p0]
        cuty=beamy[p0,:]
            
        for b in band:
            beam_cutx[:,0,b]=beams[:,p0,b]
            beam_cuty[:,0,b]=beams[p0,:,b]
        
            scalbeam=beam_regrid(beamx,beamy,beams[:,:,b],nu0_bm[b],nuL_bm[b])
            beam_cutx[:,1,b]=scalbeam[:,p0]
            beam_cuty[:,1,b]=scalbeam[p0,:]
            
            scalbeam=beam_regrid(beamx,beamy,beams[:,:,b],nu0_bm[b],nuU_bm[b])
            beam_cutx[:,2,b]=scalbeam[:,p0]
            beam_cuty[:,2,b]=scalbeam[p0,:]

        plotbeam_cut(cutx,cuty,beam_cutx,beam_cuty,brad,fsuff)
            
    ################################
    ##  Azimuthally smooth beams  ##
    ################################
    print 'Making azimuthally smoothed beam...'
    radarr=arange(0,brad,1.0)
    nrad=radarr.size
    beam_sm=zeros((nrad,3))
    beam_min=zeros((nrad,3))
    beam_max=zeros((nrad,3))
    beam_sd=zeros((nrad,3))
    # for the azimuthally smoothed measured profile, with no stitching
    files=['../Inputs/PSW_Measured.txt','../Inputs/PMW_Measured.txt','../Inputs/PLW_Measured.txt']
    
    # for the azimuthally smoothed profile with theoretical component
    #files=['../Inputs/PSW_Alex's.txt','../Inputs/PMW_Alex's.txt','../Inputs/PLW_Alex's.txt']
    
    # for the azimuthally smoothed theoretical profile 
    #files=['../Inputs/PSW_Theoretical.txt','../Inputs/PMW_Theoretical.txt','../Inputs/PLW_Theoretical.txt']    
    for b in band:
        #(beam_sm[:,b],beam_min[:,b],beam_max[:,b],beam_sd[:,b])=beam_azsm(beamx,beamy,beams[:,:,b],radarr,retall=True)
        f = open(files[b], 'r')
        line = f.readlines()
        for i in range(brad):
            beam_sm[i,b]=(float(line[i]))
        
    
    #if bazsm:
     #   plotbeam_sm(r-adarr,beam_sm,beam_min,beam_max,beam_sd,brad,fsuff)
    #beam_sm=array(beam_sm)

    ##########################################
    ##  Calculate effective beam frequency  ##
    ##########################################
    #brad=350.
    print 'Calculating Effective frequencies...'
    aprec=0.0001
    nueff=get_effnu_az(radarr,beam_sm,rsrfarr,nuarr,nuc,brad,aprec=0.0001, verbose = None)
    print 'Effective wavelengths: [%.2f, %.2f, %.2f] um' % (c/nueff[0]*1.e6,c/nueff[1]*1.e6,c/nueff[2]*1.e6)
    
    thetarr=concatenate(([1.e-8],arange(1.,10.,1.)*1.e-3, \
                        arange(1.,10.,1.)*1.e-2, \
                        arange(1.,10.,1.)*1.e-1, \
                        arange(1.,10.,1.), \
                        arange(1.,10.,1.)*10.,\
                        arange(1.,11.,1.)*100.,[10000]))
    ithi=array(where(thetarr == max(thetarr)))[0]
    print 'ithi= %d'%ithi
    #thetarr=array((1.e-3,1000.))
    #ithi=1.

    
    nth=thetarr.size

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
    area_bmth=beamarea_az_th_nu(nuarr,thetarr,radarr,beam_sm,nueff,ilim_2,dnu_a=dnu_a,verbose=True)
    area_bm=beamarea_az_nu(nuarr,radarr,beam_sm,nueff,ilim_2,dnu_a=dnu_a,verbose=False)

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
    area_eff=zeros((nalph,nth,3))
    for t in arange(nth):
        for b in band:
            area_eff[:,t,b]=beamarea_eff(alpharr,nuc[b],nuarr,area_bmth[:,t,b],rsrfarr[:,b])
    
    if pareaa:
        print 'Plotting Area vs. alpha'
        plotareaa(alpharr,area_eff,areas,ithi,fsuff,unit=aunit)
    
    #convert area to sr    
    area_bmth_sr=arcsec2sr(area_bmth)
    area_bm_sr=arcsec2sr(area_bm)
    
    if pareath:
        print 'Plotting Area vs. theta'
        areas_th=zeros((3,nth))
        for t in arange(nth):
            for b in band:
                areas_th[b,t]=beamarea_az_th(radarr,beam_sm[:,b],thetarr[t])
            print '%.3g": %.2g , %.2g , %.2g'%(thetarr[t],areas_th[0,t],areas_th[1,t],areas_th[2,t])
        plotareath(thetarr,areas_th,fsuff,unit=aunit)                
    
    if parear or pbcomp:
        print 'Computing beam area profiles...'
        area_rarr=arange(0,brad,10)
        nrad_a=area_rarr.size
        areas_r=zeros((nrad_a,3))
        for b in band:
            print 'Band %d...' % b
            for r in arange(nrad_a):
                areas_r[r,b]=beamarea(beamx,beamy,beams[:,:,b],area_rarr[r])
    
        print 'Areas: [%.2f , %.2f , %.2f] sq-arcsec' % (areas_r[nrad_a-1,0],areas_r[nrad_a-1,1],areas_r[nrad_a-1,2])
        
        if parear:
            plotarear(area_rarr,areas_r,fsuff,unit=aunit)
    
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
            beam_mod[:,b]=modbeam(radarr,beam_sm[:,b],nuarr,rsrfarr[:,b],nueff[b],alpha_nep[b])
            beam_comp[:,b]=(beam_mod[:,b]-beam_sm[:,b])/areas[b]
            areas_nep[:,b]=modbeam_area(radarr,beam_sm[:,b],area_rarr2)
            areas_mod[:,b]=modbeam_area(radarr,beam_mod[:,b],area_rarr2)
            areas_rel[:,b]=(areas_mod[:,b]-areas_nep[:,b])/areas_nep[:,b]
            #areas_rel[:,b]=(modbeam_area[:,b]-areas_r[:,b])/areas_r[:,b]
            #print 'areas (Nep,Mod,Rel):',areas_nep[:,b],areas_mod[:,b],areas_rel[:,b]
            print 'min/max',b,min(areas_rel[:,b]),max(areas_rel[:,b])
        plotbeamcomp(radarr,beam_sm,beam_mod,area_rarr2,areas_rel,fsuff)
    
    ############################
    ##  Calculate K5 and Kc2  ##
    ############################
    
    Kc1=zeros((nalph,3))
    K5f=zeros((nalph,3))
    K5=zeros((nalph,nth,3))
    Kc2=zeros((nalph,nth,3))
    Kc3=zeros((nalph,nth,3))
    
    for b in band:
        print 'band %d'%b
        for a in arange(nalph):
            ##calculate fully extended sources
            K5f[a,b]=calc_k5(alpharr[a],nuarr,rsrfarr[:,b],apfarr[:,b],area_bm_sr[:,b],nuc[b])
            Kc1[a,b]=K5f[a,b]/K4pip[b]
            for t in arange(nth):
                K5[a,t,b]=calc_k5(alpharr[a],nuarr,rsrfarr[:,b],apfarr[:,b],area_bmth_sr[:,t,b],nuc[b])
                #print '(alpha %d,theta %.3f): %g , %f'%(alpharr[a],thetarr[t],K5[a,t,b],log10(K5[a,t,b]))
                Kc2[a,t,b]=K5[a,t,b]/(K4pip[b])#/arcsec2sr(areas[b]))
                Kc3[a,t,b]=K5[a,t,b]/K4pip[b]*srcarea1[t]
    
    ##################################################
    ## K5f: Jy/sr per Jy_meas for fully extended source [f(alpha)]
    ## Kc1: Jy/sr per Jy_pip for fully extended source [f(alpha)]
    ## K5: Jy/sr per Jy_meas for partial source [f(alpha,theta)]
    ## Kc2: Jy/sr per Jy_pip for partial source [f(alpha,theta)]
    ## Kc3: Jy per Jy_pip for partial source [f(alpha,theta)]
    ##################################################
    
    if pcorr:
        #plotK5(arr_a,arr_t,K5,fsuff)
        #plotKc2(arr_a,arr_t,Kc2,fsuff)

        aplot=6 #alpha=-1
        iths=0
        thsmall=thetarr[iths]
        ithb=nth-1
        thbig=thetarr[ithb]
        
        plotK5_va(alpharr,thetarr,K5,K5f,aplot,fsuff)
        plotKc2_va(alpharr,thetarr,Kc2,Kc1,aplot,fsuff)
        plotKc3_va(alpharr,thetarr,Kc3,K4,aplot,fsuff)
        plotKc_all(alpharr,Kc,Kce,Kc2[:,nth-1,:],Kc3[:,0,:],fsuff)
        plotKcKc3_thsmall(alpharr,Kc,Kc3[:,iths,:],thsmall,fsuff)
        plotKc1Kc2_thbig(alpharr,Kc1,Kc2[:,ithb,:],thbig,fsuff)
        #plotsrcarea(thetarr,srcarea1,srcarea2,fsuff)
        #plotKc3(arr_a,arr_t,Kc3,fsuff)        
        #plotKc2(alpharr,Kc,Kce,Kc2[:,:,ithi],fsuff)
        #plotKc2rel(alpharr,K4pip,Kc,Kce,Kc2[:,:,ithi],fsuff)
           
    ######################
    ##  Write to files  ##
    ######################
    print 'Writing K4,K4e,Kc,Kce,K5f,Kc1 to files...'
    
    file_K4='../Outputs/K4'+fsuff+'.csv'
    file_K4e='../Outputs/K4e'+fsuff+'.csv'
    file_Kc='../Outputs/Kc'+fsuff+'.csv'
    file_Kce='../Outputs/Kce'+fsuff+'.csv'
    file_K5f='../Outputs/K5f'+fsuff+'.csv'
    file_Kc1='../Outputs/Kc1'+fsuff+'.csv'
    file_tex='../Outputs/Kc_Kce'+fsuff+'_tex.csv'
    
    fK4=open(file_K4,'w')
    fK4e=open(file_K4e,'w')
    fKc=open(file_Kc,'w')
    fKce=open(file_Kce,'w')
    fK5f=open(file_K5f,'w')
    fKc1=open(file_Kc1,'w')
    ftex=open(file_tex,'w')
    
    line='#alpha , PSW , PMW , PLW \n'
    fK4.write(line)
    fK4e.write(line)
    fKc.write(line)
    fKce.write(line)
    fK5f.write(line)
    fKc1.write(line)
        
    for a in arange(nalph):
        line='%.1f , %.4f , %.4f , %.4f \n' % (alpharr[a],K4[a,0],K4[a,1],K4[a,2])
        fK4.write(line)
        line='%.1f , %.4f , %.4f , %.4f \n' % (alpharr[a],K4e[a,0],K4e[a,1],K4e[a,2])
        fK4e.write(line)
        line='%.1f , %.4f , %.4f , %.4f \n' % (alpharr[a],Kc[a,0],Kc[a,1],Kc[a,2])
        fKc.write(line)
        line='%.1f , %.4f , %.4f , %.4f \n' % (alpharr[a],Kce[a,0],Kce[a,1],Kce[a,2])
        fKce.write(line)
        line='%.1f , %.5g , %.5g , %.5g \n' % (alpharr[a],K5f[a,0],K5f[a,1],K5f[a,2])
        fK5f.write(line)
        line='%.1f , %.5g , %.5g , %.5g \n' % (alpharr[a],Kc1[a,0],Kc1[a,1],Kc1[a,2])
        fKc1.write(line)
        line='%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n' % (alpharr[a],Kc[a,0],Kc[a,1],Kc[a,2],Kce[a,0],Kce[a,1],Kce[a,2])
        ftex.write(line)
        
    fK4.close()
    fK4e.close()
    fKc.close()
    fKce.close()
    fK5f.close()
    fKc1.close()
    ftex.close()

    print 'Writing K5,Kc2,Kc3 to files...'
    
    bname=['PSW','PMW','PLW']
    for b in band:
        #file names
        file_K5_at='../Outputs/K5_'+bname[b]+fsuff+'.csv'
        file_Kc2_at='../Outputs/Kc2_'+bname[b]+fsuff+'.csv'
        file_Kc3_at='../Outputs/Kc3_'+bname[b]+fsuff+'.csv'

        #open files        
        fK5at=open(file_K5_at,'w')
        fKc2at=open(file_Kc2_at,'w')
        fKc3at=open(file_Kc3_at,'w')

        #column headings
        fK5at.write('#theta/alpha')
        fKc2at.write('#theta/alpha')
        fKc3at.write('#theta/alpha')

        for a in arange(nalph):
            fK5at.write(',%.1f'%(alpharr[a]))
            fKc2at.write(',%.1f'%(alpharr[a]))
            fKc3at.write(',%.1f'%(alpharr[a]))
        fK5at.write('\n') #end-of-line
        fKc2at.write('\n') #end-of-line
        fKc3at.write('\n') #end-of-line

        #rows
        for t in arange(nth):
            fK5at.write('%.3f'%thetarr[t]) #row headings
            fKc2at.write('%.3f'%thetarr[t]) #row headings
            fKc3at.write('%.3f'%thetarr[t]) #row headings
            for a in arange(nalph):
                fK5at.write(',%g'%K5[a,t,b]) #conv factors
                fKc2at.write(',%g'%Kc2[a,t,b]) #conv factors
                fKc3at.write(',%g'%Kc3[a,t,b]) #conv factors
            fK5at.write('\n')
            fKc2at.write('\n')
            fKc3at.write('\n')

        #close files
        fK5at.close()
        fKc2at.close()
        fKc3at.close()
    
    #Show plots
    print 'Displaying plots...'
    if plotarr.any():
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
    plot.plot(alpharr,Kc[:,0],'b-',label=r'$K_\mathrm{c,PS}$ (PSW)')
    plot.plot(alpharr,Kc[:,1],'g-',label=r'$K_\mathrm{c,PS}$ (PMW)')
    plot.plot(alpharr,Kc[:,2],'r-',label=r'$K_\mathrm{c,PS}$ (PLW)')
    plot.plot(alpharr,Kce[:,0],'b--',label=r'$K_\mathrm{c,ext}$ (PSW)')
    plot.plot(alpharr,Kce[:,1],'g--',label=r'$K_\mathrm{c,ext}$ (PMW)')
    plot.plot(alpharr,Kce[:,2],'r--',label=r'$K_\mathrm{c,ext}$ (PLW)')
        
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel('Colour Correction factor')
    plot.title('Colour correction factors')
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
    plot.legend(loc='lower right',numpoints=1,ncol=2)
    
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
    plot.plot(alpharr,K4[:,0],'b-',label='Pt Src. (PSW)')
    plot.plot(alpharr,K4[:,1],'g-',label='Pr Src. (PMW)')
    plot.plot(alpharr,K4[:,2],'r-',label='Pt Src. (PLW)')
    plot.plot(alpharr,K4e[:,0],'b--',label='Ext. (PSW)')
    plot.plot(alpharr,K4e[:,1],'g--',label='Ext. (PMW)')
    plot.plot(alpharr,K4e[:,2],'r--',label='Ext. (PLW)')
    
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'K$_4$ Correction factor')
    #plot.title('Current correction factors')
    plot.legend(ncol=2,loc='lower right')
    
    plot.savefig('../Outputs/K4corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/K4corr'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotK5(arr_a,arr_t,K5,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(10)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(K5[:,:,b])),max(log10(K5[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(K5[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1: plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Spectral index ($\alpha$)')
        plot.colorbar(format='%.1f')
    
    plot.savefig('../Outputs/K5corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/K5corr'+fsuff+'.eps',transparent=False)
    
    #############################################################

def plotKc2(arr_a,arr_t,Kc2,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(11)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(Kc2[:,:,b])),max(log10(Kc2[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(Kc2[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1:plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Spectral index ($\alpha$)')
        plot.colorbar(format='%.1f')
    
    plot.savefig('../Outputs/Kc2corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/Kc2corr'+fsuff+'.eps',transparent=False)
    
    #############################################################

def plotKc3(arr_a,arr_t,Kc3,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(12)
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
        plot.colorbar()
    
    plot.savefig('../Outputs/Kc3corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/Kc3corr'+fsuff+'.eps',transparent=False)
    
    #############################################################


def plotareaa(alpharr,areaeff,areas,ith,fsuff,unit='sr'):
    
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
        ytitu='[sr]'
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
    plot.plot(alpharr,areaeff[:,ith,0],'b-',label=r'$\Omega_\mathrm{eff}$')
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
    plot.plot(alpharr,areaeff[:,ith,1],'g-',label=r'$\Omega_\mathrm{eff}$')
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
    plot.plot(alpharr,areaeff[:,ith,2],'r-',label=r'$\Omega_\mathrm{eff}$')
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
        ytit='Area [sr]'
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

def plotK5_va(alpharr,thetarr,K5,K5f,aplot,fsuff):

    plot.figure(16)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    for b in range(3):
        plot.plot(thetarr,K5[aplot,:,b]/1.e6,color=cols[b],label=tit[b])
        plot.axhline(K5f[aplot,b]/1.e6,linestyle='--',color=cols[b],label=tit[b]+' (Fully Ext.)')
    plot.xlabel('Source FWHM ["]')
    plot.ylabel(r'Conversion $K_5$ (MJy/sr per Jy measured)')
    plot.legend(loc='upper right')
    plot.yscale('log')
    plot.xscale('log')
    plot.title(r'$\alpha='+str(alpharr[aplot]))
    
    filename='../Outputs/K5_a%.1f'%(aplot)
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False,dpi=300)

def plotKc2_va(alpharr,thetarr,Kc2,Kc1,aplot,fsuff):

    plot.figure(17)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    for b in range(3):
        plot.plot(thetarr,Kc2[aplot,:,b]/1.e6,color=cols[b],label=tit[b])
        plot.axhline(Kc1[aplot,b]/1.e6,linestyle='--',color=cols[b],label=tit[b]+' (Full Ext.)')
    plot.xlabel('Source FWHM ["]')
    plot.ylabel(r'Correction $K_\mathrm{c,2}$ (MJy/sr per Jy from pipeline)')
    plot.legend(loc='upper right')
    plot.yscale('log')
    plot.xscale('log')
    plot.title(r'$\alpha='+str(alpharr[aplot]))
    
    filename='../Outputs/Kc2_a%.1f'%(aplot)
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False,dpi=300)

def plotKc3_va(alpharr,thetarr,Kc3,K4,aplot,fsuff):

    plot.figure(18)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']    
    for b in range(3):
        plot.plot(thetarr,Kc3[aplot,:,b],color=cols[b],label=tit[b])
        plot.axhline(K4[aplot,b],linestyle='--',color=cols[b],label=tit[b]+' (Pt. Src.)')
    plot.xlabel('Source FWHM ["]')
    plot.ylabel(r'Correction $K_\mathrm{c,3}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.legend(loc='upper left')
    plot.yscale('linear')
    plot.xscale('log')
    plot.ylim(0,10)
    plot.xlim(1,100)
    plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))    
    
    filename='../Outputs/Kc3_a%.1f'%(aplot)
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False)

def plotKc_all(alpharr,Kc,Kce,Kc2,Kc3,fsuff):
   
   plot.figure(19)
   plot.clf()
   tit=['PSW','PMW','PLW']
   cols=['b','g','r']
   legKc,=plot.plot(alpharr,Kc[:,0],color='k',linestyle=':')
   legKce,=plot.plot(alpharr,Kce[:,0],color='k',linestyle='--')
   legKc3,=plot.plot(alpharr,Kc3[:,0],color='k',linestyle='-')
   legpsw,=plot.plot(alpharr,Kc3[:,0],color='b',linestyle='-')
   legpmw,=plot.plot(alpharr,Kc3[:,0],color='g',linestyle='-')
   legplw,=plot.plot(alpharr,Kc3[:,0],color='r',linestyle='-')
   for b in range(3):
       plot.plot(alpharr,Kc[:,b],color=cols[b],linestyle=':')
       plot.plot(alpharr,Kce[:,b],color=cols[b],linestyle='--')
       #plot.plot(alpharr,Kc2[:,b],color=cols[b],linestyle='-')
       plot.plot(alpharr,Kc3[:,b],color=cols[b],linestyle='-')
   plot.xlabel(r'Spectral index $\alpha')
   plot.ylabel('Correction factor')
   plot.gca
   plot.legend((legKc,legKce,legKc3,legpsw,legpmw,legplw),
               (r'$K_\mathrm{c}$ (pt src)',r'$K_\mathrm{c,e}$ (fake-ext)',r'$K_\mathrm{c,3}$ (fully-ext)','PSW','PMW','PLW'),
                loc='lower center',ncol=2)
   plot.gca().lines.remove(legKc)
   plot.gca().lines.remove(legKce)
   plot.gca().lines.remove(legKc3)
   plot.gca().lines.remove(legpsw)
   plot.gca().lines.remove(legpmw)
   plot.gca().lines.remove(legplw)

   plot.savefig('../Outputs/Kc_all'+fsuff+'.png',transparent=False,dpi=300)
   plot.savefig('../Outputs/Kc_all'+fsuff+'.eps',transparent=False)
     
def plotKcKc3_thsmall(alpharr,Kc,Kc3,thplot,fsuff):

    plot.figure(20)
    plot.clf()
    plot.subplot(2,1,1)
    tit=['PSW','PMW','PLW']
    cols=['b','g','r']
    for b in range(3):
        plot.plot(alpharr,Kc[:,b],color=cols[b],linestyle='-',label=tit[b]+r' ($\theta=0$)')
        plot.plot(alpharr,Kc3[:,b],color=cols[b],linestyle='--',label=tit[b]+r'( $\theta=%g$)'%(thplot))
    plot.ylabel(r'Correction factors $K_\mathrm{c}$,$K_\mathrm{c,3}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Spectral index $\alpha$')
    plot.legend(loc='lower center',ncol=3)
    
    plot.subplot(2,1,2)
    for b in range(3):
        plot.plot(alpharr,(Kc[:,b]-Kc3[:,b])/Kc[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$(K_\mathrm{c}-K_\mathrm{c,3})/K_\mathrm{c}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Spectral index $\alpha$')
    
    plot.savefig('../Outputs/KcKc3_thsmall'+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/KcKc3_thsmall'+fsuff+'.eps',transparent=False)
    
def plotKc1Kc2_thbig(alpharr,Kc2,Kc1,thplot,fsuff):

    plot.figure(21)
    plot.clf()
    plot.subplot(2,1,1)
    tit=['PSW','PMW','PLW']
    cols=['b','g','r']
    for b in range(3):
        plot.plot(alpharr,Kc2[:,b]/1.e6,color=cols[b],linestyle='-',label=tit[b]+r' ($\theta=\infty$)')
        plot.plot(alpharr,Kc1[:,b]/1.e6,color=cols[b],linestyle='--',label=tit[b]+r'( $\theta=%g$)'%(thplot))
    plot.ylabel(r'Correction factors $K_\mathrm{c,1}$,$K_\mathrm{c,2}$ (MJy/sr per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Spectral index $\alpha$')
    plot.legend(loc='lower center',ncol=3)
    
    plot.subplot(2,1,2)
    for b in range(3):
        plot.plot(alpharr,(Kc1[:,b]-Kc2[:,b])/Kc1[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$(K_\mathrm{c,1}-K_\mathrm{c,2})/K_\mathrm{c,1}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Spectral index $\alpha$')
        
    plot.savefig('../Outputs/Kc1Kc2_thbig'+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/Kc1Kc2_thbig'+fsuff+'.eps',transparent=False)

if __name__ == "__main__":
    maincode()