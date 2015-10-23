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
from beam import getnewbeamprofile,comb_beam,get_effnu_az_newbeam,beamprof_nu
from beam import arcsec2sr,avgbeam,modbeam_new,modbeam_area,measbeamtemp_new
from beam import beamarea_az,beamarea_az_nu_new,beamarea_az_th_nu_new,beamarea_az_th
from rsrf import getrsrf,getapf,bandedge
from calc_conv import calc_k4,calc_k4temp,calc_k5,calc_k5temp

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
    parser.add_argument('--beta', action='store',default=2.,dest='beta',type=float,help='Greybody temperature spectral index')
    parser.add_argument('--nepmod', action='store',default='ESA4',dest='nepmod',help='Model to use for Neptune spectrum ["ESA2"|"ESA4"]')

    parser.add_argument('--pcorrt', action='store_true',default=False,dest='pcorr',help='Set to plot colour correction factors')
    parser.add_argument('--pareat', action='store_true',default=False,dest='pareat',help='Set to plot beam area against temperature')
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
    beta=args.beta
    nepmod=args.nepmod
    
    if pall:
        pcorr=True
        pareat=True
    else:
        pcorr=args.pcorr
        pareat=args.pareat

    if nepmod=="ESA2":
        alpha_nep=array((1.26,1.39,1.45))
    elif nepmod=="ESA4":
        alpha_nep=array((1.29,1.42,1.47))
    else:
        print 'Invalid value for NEPMOD: "%s". Must be "ESA2" or "ESA4"'
        sys.exit()
    
    fsuff='_theta_newBprofBS_temp_beta%.1f_Br%d_Ind%.2f_%s'%(beta,int(brad),ind,nepmod)
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
    
    lsum='Greybody spectral index: %.1f \n'%beta
    fsumm.write(lsum)

    
    plotarr=array([pcorr,pareat])

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
            print 'Band %d spectral resolution (nuc/dnu): %.2f'%(b,nuc[b]/(nulim[b,1]-nulim[b,0]))
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
    sys.exit()
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

    lsum='\nPassbands:\n'
    fsumm.write(lsum)
    lsum='Bandedge (PSW) [THz]: %.4f, %.4f \n' % (nulim[0,0]*1.e-12,nulim[0,1]*1.e-12)
    fsumm.write(lsum)
    lsum='Bandedge (PMW) [THz]: %.4f, %.4f \n' % (nulim[1,0]*1.e-12,nulim[1,1]*1.e-12)
    fsumm.write(lsum)
    lsum='Bandedge (PLW) [THz]: %.4f, %.4f \n' % (nulim[2,0]*1.e-12,nulim[2,1]*1.e-12)
    fsumm.write(lsum)

    ####################
    ##  Calculate K4  ##
    ####################
    
    print 'Calculating K4t,Kct...'
    temparr=arange(1.,41.,0.1)
    ntemp=temparr.size
    
    K4t=zeros((ntemp,3))
    K4te=zeros((ntemp,3))
    
    K4pip=zeros(3)
    alphapip=-1.
    
    Kct=zeros((ntemp,3))
    Kcte=zeros((ntemp,3))
    KmonP=zeros((ntemp,3))
    KcolP=zeros((ntemp,3))
    
    for b in band:
        #Calculate K4 (pipeline)
        K4pip[b] = calc_k4(alphapip,rsrfarr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
        #Calculate K4 (extended,pipeline)
        #K4epip[b] = calc_k4(alphapip,rsrfearr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
        
        for tt in range(ntemp):
            #Calculate K4t (point source)
            K4t[tt,b] = calc_k4temp(temparr[tt],beta,rsrfarr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
            #Calculate k4 (fake-extended)
            K4te[tt,b] = calc_k4temp(temparr[tt],beta,rsrfearr[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
            
            #Calculate correction factor (poitn & fake-extended)
            Kct[tt,b] = K4t[tt,b] / K4pip[b]
            Kcte[tt,b] = K4te[tt,b] / K4pip[b]
            
            KmonP[tt,b]=K4t[tt,b]
            KcolP[tt,b]=Kct[tt,b]
            
    if pcorr:
        plotK4t(temparr,K4t,K4te,beta,fsuff)
        plotKct(temparr,Kct,Kcte,beta,fsuff)

    #####################################################
    ## K4t: Jy per Jy_meas for point source [f(alpha)]
    ## K4pip: Jy_pip per Jy_meas for point source (assumes alpha=-1)
    ## K4te: Jy per Jy_meas for extended source [f(alpha)]
    ## Kct: Jy per Jy_pip for point source [f(alpha)]
    ## Kcte: Jy per Jy_pip for extended source [f(alpha)]
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
    
    print 'Beam areas: [%.2f,%.2f,%.2f]'%(areas[0],areas[1],areas[2])
    
    lsum='\nBeams:\n'
    fsumm.write(lsum)
    lsum='Beam Areas [sq. arcsec]: %.2f, %.2f, %.2f \n'%(areas[0], areas[1], areas[2])
    fsumm.write(lsum)

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
        nueff=get_effnu_az_newbeam(radarr,beam_scl,beam_fix,rsrfarr,nuarr,nuc,brad,alpha_nep,aprec=aprec,verbose=True,ind=ind)
    print 'Effective wavelengths:  [%.2f, %.2f, %.2f] um' % (c/nueff[0]*1.e6,c/nueff[1]*1.e6,c/nueff[2]*1.e6)
    #print 'Calculating Effective frequencies (beam map)...'
    #nueff2=get_effnu(beamx,beamy,beams,rsrfarr,nuarr,nuc,brad,aprec=aprec,verbose=True)
    #print 'Effective wavelengths: [%.2f, %.2f, %.2f] um' % (c/nueff2[0]*1.e6,c/nueff2[1]*1.e6,c/nueff2[2]*1.e6)
    
    lsum='Effective frequencies [THz]: %.5f, %.5f, %.5f \n'%(nueff[0]*1.e-12,nueff[1]*1.e-12,nueff[2]*1.e-12)    
    fsumm.write(lsum)
    lsum='Effective/nominal frequencies: %.5f, %.5f, %.5f \n'%(nueff[0]/nuc[0], nueff[1]/nuc[1],nueff[2]/nuc[2])
    fsumm.write(lsum)
    
    ### thetarr is in FWHM
    ###########################
    ### Long list of thetas ###
    ###########################
#    thetarr=concatenate(([1.e-8],arange(1.,10.,0.2)*1.e-3, \
#                        arange(1.,10.,0.2)*1.e-2, \
#                        arange(1.,10.,0.2)*1.e-1, \
#                        arange(1.,10.,0.2), \
#                        arange(1.,10.,0.2)*10.,\
#                        arange(1.,11.,0.2)*100.,[10000]))
#    #theta values to print to file
#    ithp=array((0,28,32,37,38,39,41,46,55,56)) ##for long theta list

    ############################
    ### Short list of thetas ###
    ############################
    thetarr=array((1.e-8,1.e-3,1.e-2,1.e-1,1.e0,1.e1,1.e2,1.e3,1.e4))
    #theta values to print to file
    ithp=arange(len(thetarr)) ##for short list

    ithi=array(where(thetarr == max(thetarr)))[0]
    
    print 'ithi= %f'%ithi
    #thetarr=array((1.e-3,1000.))
    #ithi=1.
    
    nthp=ithp.size
    nth=thetarr.size

    #for t in range(nth):
    #    print t,thetarr[t]
    #sys.exit()

    srcarea1=arcsec2sr(pi*(thetarr/(2.*sqrt(log(2.))))**2)
    #srcarea2=zeros(nth)
    #for t in arange(nth):
    #    srcarea2[t]=srcarea(thetarr[t])
    
    nth=thetarr.size
    arr_a=zeros((ntemp,nth))
    arr_t=zeros((ntemp,nth))
    for tt in arange(ntemp):
        for th in arange(nth):
            arr_a[tt,th]=temparr[tt]
            arr_t[tt,th]=thetarr[th]

    
    print 'Calculating beam area over frequency...'
    #Set frequency arrays for beam area calculations
    dnu_a=10
    area_bmth=beamarea_az_th_nu_new(nuarr,thetarr,radarr,beam_scl,beam_fix,nueff,ilim_2,dnu_a=dnu_a,verbose=True,ind=ind)
    area_bm=beamarea_az_nu_new(nuarr,radarr,beam_scl,beam_fix,nueff,ilim_2,dnu_a=dnu_a,verbose=False,ind=ind)

#    if pareaf:
#        avgnu_bm=zeros((nth,3))
#        avgarea_bm=zeros((nth,3))
#        print 'area_bmth:',area_bmth.shape
#        for th in arange(nth):
#            for b in band:
#                (avgnu_bm[tt,b],avgarea_bm[th,b])=avgbeam(nuarr,area_bmth[:,th,b])
#        print 'Plotting Area vs. freq'
#        plotareaf(nuarr,area_bmth,ilim_2,avgnu_bm,avgarea_bm,ithi,fsuff,unit=aunit)

    print 'Calculating effective beam area...'
    area_efft=zeros((ntemp,3))
    for b in band:
        print 'Band %d...'%(b)
        for tt in arange(ntemp):
            area_efft[tt,b]=measbeamtemp_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],temparr[tt],beta,ind=ind)
    filebeam='../Outputs/beamarea'+fsuff+'.csv'
    fbm=open(filebeam,'w')
    line='#temp , PSW , PMW , PLW \n'
    fbm.write(line)
    for tt in range(ntemp):
        line='%.2f, %.5f, %.5f, %.5f\n'%(temparr[tt],area_efft[tt,0],area_efft[tt,1],area_efft[tt,2])
        fbm.write(line)
    fbm.close()
    
    if pareat:
        print 'Plotting Area vs. alpha'
        plotareat(temparr,area_efft,areas,beta,fsuff,unit=aunit)
        
    #convert area to sr    
    area_bmth_sr=arcsec2sr(area_bmth)
    area_bm_sr=arcsec2sr(area_bm)
    area_efft_sr=arcsec2sr(area_efft)
    
#    if pareath:
#        print 'Plotting Area vs. theta'
#        areas_th=zeros((3,nth))
#        for t in arange(nth):
#            for b in band:
#                areas_th[b,t]=beamarea_az_th(radarr,beam_cmb[:,b],thetarr[t])
#            print '%.3g": %.2g , %.2g , %.2g'%(thetarr[t],areas_th[0,t],areas_th[1,t],areas_th[2,t])
#        plotareath(thetarr,areas_th,fsuff,unit=aunit)                
    
#    if parear or pbcomp:
#        print 'Computing beam area profiles...'
#        area_rarr=arange(0,brad+10,10)
#        nrad_a=area_rarr.size
#        areas_r=zeros((nrad_a,3))
#        for b in band:
#            print 'Band %d...' % b
#            for r in arange(nrad_a):
#                #areas_r[r,b]=beamarea_az(radarr,beam_sm[:,b],brad=area_rarr[r])
#                radx=float(area_rarr[r])
#                areas_r[r,b]=beamarea_az(radarr,beam_cmb[:,b],brad=radx)
#                #print '  Radius %.1f: Area %.2f'%(radx,areas_r[r,b])
#            print min(beam_cmb[:,b]),max(beam_cmb[:,b])
#        print 'Areas: [%.2f , %.2f , %.2f] sq-arcsec' % (areas_r[nrad_a-1,0],areas_r[nrad_a-1,1],areas_r[nrad_a-1,2])
#        
#        if parear:
#            plotarear(area_rarr,areas_r,fsuff,unit=aunit)
#
#        plot.figure(30)
#        plot.plot(radarr,beam_scl[:,0],'b-')
#        plot.plot(radarr,beam_scl[:,1],'g-')
#        plot.plot(radarr,beam_scl[:,2],'r-')
#        plot.plot(radarr,beam_fix[:,0],'b--')
#        plot.plot(radarr,beam_fix[:,1],'g--')
#        plot.plot(radarr,beam_fix[:,2],'r--')
#        plot.yscale('log')
        
#    if pbcomp:
#        print 'Comparing modelled and measured areas...'
#        beam_mod=zeros((nrad,3))
#        beam_comp=zeros((nrad,3))
#        #areas_mod=zeros((nrad_a,3))
#        area_rarr2=arange(1,101,1)
#        nrad_a2=area_rarr2.size
#        areas_nep=zeros((nrad_a2,3))
#        areas_mod=zeros((nrad_a2,3))
#        areas_rel=zeros((nrad_a2,3))
#        for b in band:
#            print 'Band %d...' % b
#            beam_mod[:,b]=modbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpha_nep[b],ind=ind)
#            beam_comp[:,b]=(beam_mod[:,b]-beam_cmb[:,b])/areas[b]
#            areas_nep[:,b]=modbeam_area(radarr,beam_cmb[:,b],area_rarr2)
#            areas_mod[:,b]=modbeam_area(radarr,beam_mod[:,b],area_rarr2)
#            areas_rel[:,b]=(areas_mod[:,b]-areas_nep[:,b])/areas_nep[:,b]
#            #areas_rel[:,b]=(modbeam_area[:,b]-areas_r[:,b])/areas_r[:,b]
#            #print 'areas (Nep,Mod,Rel):',areas_nep[:,b],areas_mod[:,b],areas_rel[:,b]
#            print 'min/max',b,min(areas_rel[:,b]),max(areas_rel[:,b])
#        plotbeamcomp(radarr,beam_cmb,beam_mod,area_rarr2,areas_rel,fsuff)
    
    ############################
    ##  Calculate K5 and Kc2  ##
    ############################
    
    K5t=zeros((ntemp,3))
    KpipE=zeros(3)
    Kct2=zeros((ntemp,3))
    Kct2b=zeros((ntemp,3))
    K6t=zeros((ntemp,nth,3))
    Kct3=zeros((ntemp,nth,3))
    Kct3t=zeros((ntemp,nth,3))
    KmonE=zeros((ntemp,nth,3))
    Kuni=zeros((ntemp,3))
    KcolE=zeros((ntemp,nth,3))
    KcolEuni=zeros((ntemp,3))
    
    for b in band:
        print 'band %d'%b       
        KpipE[b]=calc_k5(alphapip,nuarr,rsrfarr[:,b],apfarr[:,b],area_bm_sr[:,b],nuc[b])
        for tt in range(ntemp):
            K5t[tt,b]=calc_k5temp(temparr[tt],beta,nuarr,rsrfarr[:,b],apfarr[:,b],area_bm_sr[:,b],nuc[b])
            Kct2[tt,b]=K5t[tt,b]/K4pip[b]
            Kct2b[tt,b]=Kct2[tt,b]*area_efft_sr[tt,b]
            Kuni[tt,b]=K5t[tt,b]
            KcolEuni[tt,b]=K5t[tt,b]/KpipE[b]
            for th in arange(nth):
                K6t[tt,th,b]=calc_k5temp(temparr[tt],beta,nuarr,rsrfarr[:,b],apfarr[:,b],area_bmth_sr[:,th,b],nuc[b])
                #print '(alpha %d,theta %.3f): %g , %f'%(alpharr[a],thetarr[t],K5[a,t,b],log10(K5[a,t,b]))
                Kct3[tt,th,b]=K6t[tt,th,b]/(K4pip[b])#/arcsec2sr(areas[b]))
                Kct3t[tt,th,b]=K6t[tt,th,b]/K4pip[b]*srcarea1[th]
                KmonE[tt,th,b]=K6t[tt,th,b]
                KcolE[tt,th,b]=K6t[tt,th,b]/KpipE[b]

            
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
        
        plotK6t_vt(temparr,thetarr,K6t,K5t,aplot,beta,fsuff)
        plotKct3_vt(temparr,thetarr,Kct3,Kct2,aplot,beta,fsuff)
        plotKct3t_vt(temparr,thetarr,Kct3t,Kct,aplot,beta,fsuff)
        plotKct_all(temparr,Kct,Kcte,Kct3[:,ithb,:],Kct3t[:,0,:],beta,fsuff)
        plotKctKct3t_thsmall(temparr,Kct,Kct3t[:,iths,:],thsmall,beta,fsuff)
        plotKct2Kct3_thbig(temparr,Kct2,Kct3[:,ithb,:],thbig,beta,fsuff)
        #plotsrcarea(thetarr,srcarea1,srcarea2,fsuff)
        #plotKc3(arr_a,arr_t,Kc3,fsuff)        
        plotKct2(temparr,Kcte,Kct2,area0,beta,fsuff)
        #plotKc2rel(alpharr,Kce,Kc2,area0,fsuff)
    
    #print Kct3[:,nth-1,0]
    #print Kct2[:,0]
    #print (Kct3[:,nth-1,0]-Kct2[:,0])/Kct2[:,0]
    
    lsum='\nOutputs:\n'
    fsumm.write(lsum)
    tsum1=4 #alpha=5K
    tsum2=24 #alpha=25K
    #lsum='alpha 1: %.1f \n'%(alpharr[alpharr[tsum1]])
    #fsumm.write(line)
    #lsum='alpha 2: %.1f \n'%(alpharr[tsum2])
    #fsumm.write(line)
    lsum='K4t(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum1],K4t[tsum1,0],K4t[tsum1,1],K4t[tsum1,2])
    fsumm.write(lsum)
    lsum='K4te(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum1],K4te[tsum1,0],K4te[tsum1,1],K4te[tsum1,2])
    fsumm.write(lsum)
    lsum='Kct(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum1],Kct[tsum1,0],Kct[tsum1,1],Kct[tsum1,2])
    fsumm.write(lsum)
    lsum='Kcte(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum1],Kcte[tsum1,0],Kcte[tsum1,1],Kcte[tsum1,2])
    fsumm.write(lsum)
    lsum='K5t(%dK): %.5g, %.5g, %.5g \n'% (temparr[tsum1],K5t[tsum1,0], K5t[tsum1,1], K5t[tsum1,2])
    fsumm.write(lsum)
    lsum='Kct2(%dK): %.5g, %.5g, %.5g \n'% (temparr[tsum1],Kct2[tsum1,0], Kct2[tsum1,1], Kct2[tsum1,2])
    fsumm.write(lsum)
    lsum='Kct2b(%dK): %.5f, %.5f, %.5f \n'% (temparr[tsum1],Kct2b[tsum1,0], Kct2b[tsum1,1], Kct2b[tsum1,2])
    fsumm.write(lsum)
    lsum='K4t(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum2],K4t[tsum2,0],K4t[tsum2,1],K4t[tsum2,2])
    fsumm.write(lsum)
    lsum='K4te(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum2],K4te[tsum2,0],K4te[tsum2,1],K4te[tsum2,2])
    fsumm.write(lsum)
    lsum='Kct(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum2],Kct[tsum2,0],Kct[tsum2,1],Kct[tsum2,2])
    fsumm.write(lsum)
    lsum='Kcte(%dK): %.5f, %.5f, %.5f \n'%(temparr[tsum2],Kcte[tsum2,0],Kcte[tsum2,1],Kcte[tsum2,2])
    fsumm.write(lsum)
    lsum='K5t(%dK): %.5g, %.5g, %.5g \n'% (temparr[tsum2],K5t[tsum2,0], K5t[tsum2,1], K5t[tsum2,2])
    fsumm.write(lsum)
    lsum='Kct2(%dK): %.5g, %.5g, %.5g \n'% (temparr[tsum2],Kct2[tsum2,0], Kct2[tsum2,1], Kct2[tsum2,2])
    fsumm.write(lsum)
    lsum='Kct2b(%dK): %.5f, %.5f, %.5f \n'% (temparr[tsum2],Kct2b[tsum2,0], Kct2b[tsum2,1], Kct2b[tsum2,2])
    fsumm.write(lsum)
       
    ######################
    ##  Write to files  ##
    ######################
    print 'Writing K4t,K4te,Kct,Kcte,K5t,Kct2,Kct2b to files...'
    
    ##set filenames
    file_K4t='../Outputs/K4t'+fsuff+'.csv'
    file_K4te='../Outputs/K4te'+fsuff+'.csv'
    file_Kct='../Outputs/Kct'+fsuff+'.csv'
    file_Kcte='../Outputs/Kcte'+fsuff+'.csv'
    file_K5t='../Outputs/K5t'+fsuff+'.csv'
    file_Kct2='../Outputs/Kct2'+fsuff+'.csv'
    file_Kct2b='../Outputs/Kct2b'+fsuff+'.csv'
    
    ##params as named in documents
    file_KmonP='../Outputs/KmonP-t'+fsuff+'.csv'
    file_KcolP='../Outputs/KcolP-t'+fsuff+'.csv'
    file_Kuni='../Outputs/Kuni-t'+fsuff+'.csv'
    file_KcolEuni='../Outputs/KcolEuni-t'+fsuff+'.csv'
    
    ##other parameters
    file_beam='../Outputs/beamarea-t'+fsuff+'.csv'
    #file_fwhm='../Outputs/beamfwhm'+fsuff+'.csv'
    
    ##TeX tables
    file_tex='../Docs/Tables/Kct_Kcte'+fsuff+'.tex'
    file_tex2='../Docs/Tables/Kct2_Kct2b'+fsuff+'.tex'
    file_tex3='../Docs/Tables/Kct_Kct2'+fsuff+'.tex'
    
    fK4t=open(file_K4t,'w')
    fK4te=open(file_K4te,'w')
    fKct=open(file_Kct,'w')
    fKcte=open(file_Kcte,'w')
    fK5t=open(file_K5t,'w')
    fKct2=open(file_Kct2,'w')
    fKct2b=open(file_Kct2b,'w')
    fKmonP=open(file_KmonP,'w')
    fKcolP=open(file_KcolP,'w')
    fKuni=open(file_Kuni,'w')
    fKcolEuni=open(file_KcolEuni,'w')
    ftex=open(file_tex,'w')
    ftex2=open(file_tex2,'w')
    ftex3=open(file_tex3,'w')
    
    line='#temp , PSW , PMW , PLW \n'
    fK4t.write(line)
    fK4te.write(line)
    fKct.write(line)
    fKcte.write(line)
    fK5t.write(line)
    fKct2.write(line)
    fKct2b.write(line)
    fKmonP.write(line)
    fKcolP.write(line)
    fKuni.write(line)
    fKcolEuni.write(line)

    line='\\begin{tabular}{c|ccc|ccc}\n'
    ftex.write(line)
    ftex2.write(line)
    line='& \\multicolumn{3}{c|}{$K_\\mathrm{c,pt}$} & \\multicolumn{3}{c}{$K_\\mathrm{c,ext}$} \\\\\n'
    ftex.write(line)
    line='& \\multicolumn{3}{c|}{$K_\\mathrm{c2} \\times 10^{-7}$} & \\multicolumn{3}{c}{$K_\\mathrm{c2}^\\mathrm{bm}$} \\\\\n'
    ftex2.write(line)
    line='& \\multicolumn{3}{c|}{$K_\\mathrm{CP}$} & \\multicolumn{3}{c}{$K_\\mathrm{CE}$} \\\\\n'
    ftex3.write(line)
    line='Temp (K) & PSW & PMW & PLW & PSW & PMW & PLW \\\\\n'
    ftex.write(line)
    ftex2.write(line)
    ftex3.write(line)
    line='\\hline\n'
    ftex.write(line)
    ftex2.write(line)
    ftex3.write(line)

    Kct2_p=Kct2
    Kct2_exp=floor(log10(Kct2_p))
    Kct2_man=Kct2_p/10.**Kct2_exp
    Kct2b_p=Kct2b
    Kct2b_exp=floor(log10(Kct2b_p))
    Kct2b_man=Kct2b_p/10.**Kct2b_exp

        
    for tt in arange(ntemp):
        line='%.1f , %.9f , %.9f , %.9f \n' % (temparr[tt],K4t[tt,0],K4t[tt,1],K4t[tt,2])
        fK4t.write(line)
        line='%.1f , %.9f , %.9f , %.9f \n' % (temparr[tt],K4te[tt,0],K4te[tt,1],K4te[tt,2])
        fK4te.write(line)
        line='%.1f , %.9f , %.9f , %.9f \n' % (temparr[tt],Kct[tt,0],Kct[tt,1],Kct[tt,2])
        fKct.write(line)
        line='%.1f , %.9f , %.9f , %.9f \n' % (temparr[tt],Kcte[tt,0],Kcte[tt,1],Kcte[tt,2])
        fKcte.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (temparr[tt],K5t[tt,0],K5t[tt,1],K5t[tt,2])
        fK5t.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (temparr[tt],Kct2[tt,0],Kct2[tt,1],Kct2[tt,2])
        fKct2.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (temparr[tt],Kct2b[tt,0],Kct2b[tt,1],Kct2b[tt,2])
        fKct2b.write(line)

        line='%.1f , %.9g , %.9g , %.9g \n' % (temparr[tt],KmonP[tt,0],KmonP[tt,1],KmonP[tt,2])
        fKmonP.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (temparr[tt],KcolP[tt,0],KcolP[tt,1],KcolP[tt,2])
        fKcolP.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (temparr[tt],Kuni[tt,0],Kuni[tt,1],Kuni[tt,2])
        fKuni.write(line)
        line='%.1f , %.9g , %.9g , %.9g \n' % (temparr[tt],KcolEuni[tt,0],KcolEuni[tt,1],KcolEuni[tt,2])
        fKcolEuni.write(line)
        
        line='%.1f & %.9f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n' % (temparr[tt],Kct[tt,0],Kct[tt,1],Kct[tt,2],Kcte[tt,0],Kcte[tt,1],Kcte[tt,2])
        ftex.write(line)
        line='%.1f & %.4g & %.4g & %.4g & %.4f & %.4f & %.4f \\\\ \n' % (temparr[tt],Kct2[tt,0]/1.e7,Kct2[tt,1]/1.e7,Kct2[tt,2]/1.e7,Kct2b[tt,0],Kct2b[tt,1],Kct2b[tt,2])
        ftex2.write(line)
        line='%.1f & %.4f & %.4f & %.4f & %.4g & %.4g & %.4g \\\\ \n' % (temparr[tt],Kct[tt,0],Kct[tt,1],Kct[tt,2],Kct2[tt,0]/1.e-6,Kct2[tt,1]/1.e-6,Kct2[tt,2]/1.e-6)
        ftex3.write(line)
        
    fK4t.close()
    fK4te.close()
    fKct.close()
    fKcte.close()
    fK5t.close()
    fKct2.close()
    fKct2b.close()
    fKmonP.close()
    fKcolP.close()
    fKuni.close()
    fKcolEuni.close()

    line='\\end{tabular}\n'
    ftex.write(line)
    ftex2.write(line)
    ftex3.write(line)
    ftex.close()
    ftex2.close()
    ftex3.close()


    print 'Writing K6,Kc3,Kc3i to files...'
    
    th_p=thetarr
    th_exp=floor(log10(th_p))
    th_man=th_p/10.**th_exp
    Kct3_p=Kct3
    Kct3_exp=floor(log10(Kct3_p))
    Kct3_man=Kct3_p/10.**Kct3_exp
    Kct3t_p=Kct3t
    Kct3t_exp=floor(log10(Kct3t_p))
    Kct3t_man=Kct3t_p/10.**Kct3t_exp


    bname=['PSW','PMW','PLW']
    for b in band:
        #file names
        file_K6t_at='../Outputs/K6t_'+bname[b]+fsuff+'.csv'
        file_Kct3_at='../Outputs/Kct3_'+bname[b]+fsuff+'.csv'
        file_Kct3t_at='../Outputs/Kct3t_'+bname[b]+fsuff+'.csv'
        file_Kct3_tex='../Docs/Tables/Kct3_'+bname[b]+fsuff+'.tex'
        file_Kct3t_tex='../Docs/Tables/Kct3t_'+bname[b]+fsuff+'.tex'
        #open files        
        fK6tat=open(file_K6t_at,'w')
        fKct3at=open(file_Kct3_at,'w')
        fKct3tat=open(file_Kct3t_at,'w')
        fKct3_tex=open(file_Kct3_tex,'w')
        fKct3t_tex=open(file_Kct3t_tex,'w')

        #column headings
        fK6tat.write('#theta/temp')
        fKct3at.write('#theta/temp')
        fKct3tat.write('#theta/atemp')

        for tt in arange(ntemp):
            fK6tat.write(',%.1f'%(temparr[tt]))
            fKct3at.write(',%.1f'%(temparr[tt]))
            fKct3tat.write(',%.1f'%(temparr[tt]))
        fK6tat.write('\n') #end-of-line
        fKct3at.write('\n') #end-of-line
        fKct3tat.write('\n') #end-of-line

        #rows
        for t in arange(nth):
            fK6tat.write('%.3f'%thetarr[t]) #row headings
            fKct3at.write('%.3f'%thetarr[t]) #row headings
            fKct3tat.write('%.3f'%thetarr[t]) #row headings
            for tt in arange(ntemp):
                fK6tat.write(',%.9g'%K6t[tt,t,b]) #conv factors
                fKct3at.write(',%.9g'%Kct3[tt,t,b]) #conv factors
                fKct3tat.write(',%.9g'%Kct3t[tt,t,b]) #conv factors
            fK6tat.write('\n')
            fKct3at.write('\n')
            fKct3tat.write('\n')


        line='\\begin{tabular}{c|ccccccccccc}\n'
        fKct3_tex.write(line)
        fKct3t_tex.write(line)
        line='& \\multicolumn{10}{c}{Source FWHM, $\\theta_0$ (arcsec)} \\\\\n'
        fKct3_tex.write(line)
        fKct3t_tex.write(line)
        fKct3_tex.write('Temp (K) ')
        fKct3t_tex.write('Temp (K) ')
        for p in range(nthp):
            if thetarr[ithp[p]] < 0.1:
                fKct3_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
                fKct3t_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
            elif thetarr[ithp[p]] >= 10000.:
                fKct3_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
                fKct3t_tex.write('& $%.1f \\times 10^{%d}$ '%(th_man[ithp[p]],th_exp[ithp[p]]))
            else:
                fKct3_tex.write('& %.1f '%thetarr[ithp[p]])
                fKct3t_tex.write('& %.1f '%thetarr[ithp[p]])
        fKct3_tex.write('\\\\\n')
        fKct3t_tex.write('\\\\\n')
        fKct3_tex.write('\\hline\n')
        fKct3t_tex.write('\\hline\n')

        for tt in range(ntemp):
            fKct3_tex.write('%.1f '%temparr[tt])
            fKct3t_tex.write('%.1f '%temparr[tt])
            for p in range(nthp):
                if Kct3[tt,ithp[p],b] >= 10:
                    fKct3_tex.write('& $%.2f \\times 10^{%d}$ '%(Kct3_man[tt,ithp[p],b],Kct3_exp[tt,ithp[p],b]))
                elif Kct3[tt,ithp[p],b] <= 0.01:
                    fKct3_tex.write('& $%.2f \\times 10^{%d}$ '%(Kct3_man[tt,ithp[p],b],Kct3_exp[tt,ithp[p],b]))
                else:
                    fKct3_tex.write('& %.4f '%Kct3[tt,ithp[p],b])
                if Kct3t[tt,ithp[p],b] >= 10:
                    fKct3t_tex.write('& $%.2f \\times 10^{%d}$ '%(Kct3t_man[tt,ithp[p],b],Kct3t_exp[tt,ithp[p],b]))
                elif Kct3t[tt,ithp[p],b] <= 0.01:
                    fKct3t_tex.write('& $%.2f \\times 10^{%d}$ '%(Kct3t_man[tt,ithp[p],b],Kct3t_exp[tt,ithp[p],b]))
                else:
                    fKct3t_tex.write('& %.4f '%Kct3t[tt,ithp[p],b])
            fKct3_tex.write('\\\\\n')
            fKct3t_tex.write('\\\\\n')
        
        fKct3_tex.write('\\end{tabular}\n')
        fKct3t_tex.write('\\end{tabular}\n')

        #close files
        fK6tat.close()
        fKct3at.close()
        fKct3tat.close()
        fKct3_tex.close()
        fKct3t_tex.close()

    fsumm.close()    
    
    #Show plots
    if plotarr.any():
        print 'Displaying plots...'
        plot.show()
    
    
#===============================================================================
#     Plotting routines
#===============================================================================

    
def plotKct(temparr,Kct,Kcte,beta,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(3)
    plot.clf()
    plot.plot(temparr,Kct[:,0],'b-',label=r'$K_\mathrm{c,pt}$ (PSW)')
    plot.plot(temparr,Kct[:,1],'g-',label=r'$K_\mathrm{c,pt}$ (PMW)')
    plot.plot(temparr,Kct[:,2],'r-',label=r'$K_\mathrm{c,pt}$ (PLW)')
    plot.plot(temparr,Kcte[:,0],'b--',label=r'$K_\mathrm{c,ext}$ (PSW)')
    plot.plot(temparr,Kcte[:,1],'g--',label=r'$K_\mathrm{c,ext}$ (PMW)')
    plot.plot(temparr,Kcte[:,2],'r--',label=r'$K_\mathrm{c,ext}$ (PLW)')
        
    plot.xlabel(r'Temperature (K)')
    plot.ylabel('Colour Correction factor')
    #plot.title('Colour correction factors')
    plot.legend(ncol=2,loc='lower left')
    
    plot.savefig('../Outputs/colcorrt'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/colcorrt'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotK4t(temparr,K4t,K4te,beta,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(9)
    plot.clf()
    plot.plot(temparr,K4t[:,0],'b-',label=r'$K_\mathrm{4,pt}$ (PSW)')
    plot.plot(temparr,K4t[:,1],'g-',label=r'$K_\mathrm{4,pt}$ (PMW)')
    plot.plot(temparr,K4t[:,2],'r-',label=r'$K_\mathrm{4,pt}$ (PLW)')
    plot.plot(temparr,K4te[:,0],'b--',label=r'$K_\mathrm{4,ext}$ (PSW)')
    plot.plot(temparr,K4te[:,1],'g--',label=r'$K_\mathrm{4,ext}$ (PMW)')
    plot.plot(temparr,K4te[:,2],'r--',label=r'$K_\mathrm{4,ext}$ (PLW)')
    
    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'K$_4$ Correction factor')
    #plot.title('Current correction factors')
    plot.legend(ncol=2,loc='lower left')
    
    plot.savefig('../Outputs/K4tcorr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/K4tcorr'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotK6t(arr_a,arr_t,K6t,beta,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(10)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(K6t[:,:,b])),max(log10(K6t[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(K6t[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1: plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Temperature (K)')
        plot.colorbar(format='%.1f')
    
    plot.savefig('../Outputs/K6tcorr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/K6tcorr'+fsuff+'.eps',transparent=False)
    
    #############################################################

def plotKct3(arr_a,arr_t,Kct3,beta,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(11)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(Kct3[:,:,b])),max(log10(Kct3[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(Kct3[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1:plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Temperature (K)')
        plot.colorbar(format='%.1f')
    
    plot.savefig('../Outputs/Kct3corr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/Kct3corr'+fsuff+'.eps',transparent=False)
    
    #############################################################

def plotKct3t(arr_a,arr_t,Kct3t,beta,fsuff):
    ##############################################################
    ##plot colour correction factors
    
    tit=('PSW','PMW','PLW')
    
    plot.figure(12)
    plot.clf()
    for b in arange(3):
        plim=[min(log10(Kct3t[:,:,b])),max(log10(Kct3t[:,:,b]))]
        pr=plim[1]-plim[0]
        ncol=32
        clevs=arange(plim[0],plim[1]+pr/ncol,pr/ncol)
        plot.subplot(2,2,b+1)
        plot.contourf(arr_a[:,:-1],arr_t[:,:-1],log10(Kct3t[:,:-1,b]),levels=clevs)
        plot.title(tit[b])
        plot.yscale('log')
        if b != 1:plot.ylabel(r'$\theta_\mathrm{FWHM}$ ["]')
        if b != 0:plot.xlabel(r'Spectral index ($\alpha$)')
        plot.colorbar()
    
    plot.savefig('../Outputs/Kct3tcorr'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/Kct3tcorr'+fsuff+'.eps',transparent=False)
    
    #############################################################


def plotareat(temparr,areaefft,areas,beta,fsuff,alpha_nep,unit='sr'):
    
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
        areaefft=arcsec2sr(areaefft)*1.e8
        ytitu=r'[sr] $\times 10^8$'
    else: ytitu=r'[arcsec${^2}]$'
    
    #xtickv=arange(-4,5,1)
    #nt=xtickv.size
    #xtickl=[]
    #xtickbl=[]
    #for n in arange(nt): 
    #    xtickl.append(str(xtickv[n]))
    #    xtickbl.append('')
    
    plot.figure(13)
    plot.clf()
    plot.subplot(3,1,1)
    plot.plot(temparr,areaefft[:,0],'b-',label=r'$\Omega_\mathrm{eff}$')
    plot.axhline(y=areas[0],color='b',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    plot.axhline(y=area0[0],color='b',linestyle=':',label=r'$\Omega_0$')
    #plot.axvline(x=alpha_nep[0],color='b',linestyle=':')
    #plot.ylabel('Beam Area '+ytitu)
    #plot.xticks(xtickv,xtickbl)
    #plot.title('PSW')
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    plot.text(x0+0.05*xr,y0+0.05*yr,'PSW')
    plot.legend(loc='upper right',ncol=3)

    plot.subplot(3,1,2)
    plot.plot(temparr,areaefft[:,1],'g-',label=r'$\Omega_\mathrm{eff}$')
    plot.axhline(y=areas[1],color='g',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    plot.axhline(y=area0[1],color='g',linestyle=':',label=r'$\Omega_0$')
    plot.axvline(x=alpha_nep[1],color='g',linestyle=':')
    plot.ylabel('Beam Area '+ytitu)
    #plot.xticks(xtickv,xtickbl)
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    plot.text(x0+0.05*xr,y0+0.05*yr,'PMW')
    #plot.title('PMW')
    
    plot.subplot(3,1,3)
    plot.plot(temparr,areaefft[:,2],'r-',label=r'$\Omega_\mathrm{eff}$')
    plot.axhline(y=areas[2],color='r',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    plot.axhline(y=area0[2],color='r',linestyle=':',label=r'$\Omega_0$')
    plot.axvline(x=alpha_nep[2],color='r',linestyle=':')
    #plot.ylabel('Beam Area '+ytitu)
    plot.xlabel(r'Temperature (K)')
    #plot.xticks(xtickv,xtickl)
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    plot.text(x0+0.05*xr,y0+0.05*yr,'PLW')
    #plot.title('PLW')
    
    plot.savefig('../Outputs/beamarea_temp'+fsuff+'_'+unit+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/beamarea_temp'+fsuff+'_'+unit+'.eps',transparent=False)
    
def plotK6t_vt(temparr,thetarr,K6t,K5t,tplot,beta,fsuff):

    plot.figure(16)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    for b in range(3):
        plot.plot(thetarr,K6t[tplot,:,b]/1.e6,color=cols[b],label=r'$K_6$ (%s)'%(tit[b]))
    for b in range(3)        :
        plot.axhline(K5t[tplot,b]/1.e6,linestyle='--',color=cols[b],label=r'$K_5$ (%s)'%(tit[b]))
    plot.xlabel(r'Source FWHM, $\theta_0$ ["]')
    plot.ylabel(r'Conversion $K_6$ (MJy/sr per Jy measured)')
    plot.legend(loc='upper right',ncol=2)
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(1.e1,1.e7)
    plot.xlim(1.e-1,1.e4)
    plot.title(r'Temp=%.1fK'%(temparr[tplot]))
    
    filename='../Outputs/K6t_a%.1f'%(temparr[tplot])
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False,dpi=300)

def plotKct3_vt(temparr,thetarr,Kct3,Kct2,tplot,beta,fsuff):

    plot.figure(17)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']
    for b in range(3):
        plot.plot(thetarr,Kct3[tplot,:,b]/1.e6,color=cols[b],label=r'$K_\mathrm{c3}$ (%s)'%(tit[b]))
    for b in range(3):        
        plot.axhline(Kct2[tplot,b]/1.e6,linestyle='--',color=cols[b],label=r'$K_\mathrm{c2}$ (%s)'%(tit[b]))
    plot.xlabel(r'Source FWHM, $\theta_0$ ["]')
    plot.ylabel(r'Correction $K_\mathrm{c3}$ (MJy/sr per Jy from pipeline)')
    plot.legend(loc='upper right',ncol=2)
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(1.e1,1.e7)
    plot.xlim(1.e-1,1.e4)
    plot.title(r'Temp=%.1fK'%(temparr[tplot]))
    
    filename='../Outputs/Kct3_a%.1f'%(temparr[tplot])
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False,dpi=300)

def plotKct3t_vt(temparr,thetarr,Kct3t,Kct,tplot,beta,fsuff):

    plot.figure(18)
    plot.clf()
    cols=['b','g','r']
    tit=['PSW','PMW','PLW']    
    for b in range(3):
        plot.plot(thetarr,Kct3t[tplot,:,b],color=cols[b],label=r'$K_\mathrm{c3}^\mathrm{tot}$ (%s)'%(tit[b]))
    for b in range(3):        
        plot.axhline(Kct[tplot,b],linestyle='--',color=cols[b],label=r'$K_\mathrm{c,pt}$ (%s)'%(tit[b]))
    plot.xlabel('Source FWHM ["]')
    plot.ylabel(r'Correction $K_\mathrm{c3}^\mathrm{tot}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.legend(loc='upper left',ncol=2)
    plot.yscale('linear')
    plot.xscale('log')
    plot.ylim(0,10)
    plot.xlim(1,100)
    plot.title(r'$\alpha=%.1f$'%(temparr[tplot]))    
    
    filename='../Outputs/Kct3t_a%.1f'%(temparr[tplot])
    plot.savefig(filename+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig(filename+fsuff+'.eps',transparent=False)

def plotKct_all(temparr,Kct,Kcte,Kct3,Kct3t,beta,fsuff):
   
   plot.figure(19)
   plot.clf()
   tit=['PSW','PMW','PLW']
   cols=['b','g','r']
   legKc,=plot.plot(temparr,Kct[:,0],color='k',linestyle=':')
   legKce,=plot.plot(temparr,Kcte[:,0],color='k',linestyle='--')
   legKc3t,=plot.plot(temparr,Kct3t[:,0],color='k',linestyle='-')
   legpsw,=plot.plot(temparr,Kct3t[:,0],color='b',linestyle='-')
   legpmw,=plot.plot(temparr,Kct3t[:,0],color='g',linestyle='-')
   legplw,=plot.plot(temparr,Kct3t[:,0],color='r',linestyle='-')
   for b in range(3):
       plot.plot(temparr,Kct[:,b],color=cols[b],linestyle=':')
       plot.plot(temparr,Kcte[:,b],color=cols[b],linestyle='--')
       #plot.plot(alpharr,Kc2[:,b],color=cols[b],linestyle='-')
       plot.plot(temparr,Kct3t[:,b],color=cols[b],linestyle='-')
   plot.xlabel(r'Temperature (K)')
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

   plot.savefig('../Outputs/Kct_all'+fsuff+'.png',transparent=False,dpi=300)
   plot.savefig('../Outputs/Kct_all'+fsuff+'.eps',transparent=False)
     
def plotKctKct3t_thsmall(temparr,Kct,Kct3t,thplot,beta,fsuff):

    plot.figure(20)
    plot.clf()

    tit=['PSW','PMW','PLW']
    cols=['b','g','r']
    
    for b in range(3):
        plot.plot(temparr,1.e3*(Kct[:,b]-Kct3t[:,b])/Kct[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$10^3 \times \ (K_\mathrm{c}-K_\mathrm{c3}^\mathrm{tot})/K_\mathrm{c}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Temperature (K)')
    plot.legend(loc='upper left',ncol=3)
    plot.ylim(-4,1)
    
    plot.savefig('../Outputs/KctKct3t_thsmall'+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/KctKct3t_thsmall'+fsuff+'.eps',transparent=False)
    
def plotKct2Kct3_thbig(temparr,Kct3,Kct2,thplot,beta,fsuff):

    plot.figure(21)
    plot.clf()
    tit=['PSW','PMW','PLW']
    cols=['b','g','r']
    for b in range(3):
        plot.plot(temparr,1.e5*(Kct2[:,b]-Kct3[:,b])/Kct2[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$10^5 \times \ (K_\mathrm{c2}-K_\mathrm{c3})/K_\mathrm{c2}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Temperature (K)')

    plot.legend(loc='upper right',ncol=3)
    plot.savefig('../Outputs/Kct2Kct3_thbig'+fsuff+'.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/Kct2Kct3_thbig'+fsuff+'.eps',transparent=False)

def plotKct2(temparr,Kcte,Kct2,area0,beta,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(22)
    plot.clf()
    plot.subplot(2,1,1)
    
    lkc2,=plot.plot(temparr,Kct2[:,0],'k-',label=r'$K_\mathrm{c2}$')
    lkce,=plot.plot(temparr,Kcte[:,0]/area0[0],'k--',label=r'$K_\mathrm{c,ext}/\Omega_0$')
    lpsw,=plot.plot(temparr,Kct2[:,0],'b-',label=r'$K_\mathrm{c2}$ (PSW)')
    lpmw,=plot.plot(temparr,Kct2[:,1],'g-',label=r'$K_\mathrm{c2}$ (PMW)')
    lplw,=plot.plot(temparr,Kct2[:,2],'r-',label=r'$K_\mathrm{c2}$ (PLW)')
        
    plot.plot(temparr,Kcte[:,0]/area0[0],'b--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PSW)')
    plot.plot(temparr,Kcte[:,1]/area0[1],'g--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PMW)')
    plot.plot(temparr,Kcte[:,2]/area0[2],'r--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PLW)')
    
    #plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'Conversion factor (Jy/sr per Jy$_\mathrm{pip}$)')
    plot.ylim(ymax=1.3e8)
    #plot.title('Colour correction factors')        
    #plot.legend(ncol=2,loc='lower right')
    plot.legend((lkc2,lkce),(r'$K_\mathrm{c2}$',r'$K_\mathrm{c,ext}/\Omega_0$'),ncol=2,loc='upper right')

    plot.subplot(2,1,2)
    plot.plot(temparr,(Kcte[:,0]/area0[0])/Kct2[:,0],'b-',label='PSW')
    plot.plot(temparr,(Kcte[:,1]/area0[1])/Kct2[:,1],'g-',label='PMW')
    plot.plot(temparr,(Kcte[:,2]/area0[2])/Kct2[:,2],'r-',label='PLW')
    
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'$(K_\mathrm{c,ext}/\Omega_0)\, / \, K_\mathrm{c2}$')
    plot.ylim(1.0,1.12)
    plot.legend(ncol=3,loc='lower right')
        
    plot.savefig('../Outputs/colcorrt2'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/colcorrt2'+fsuff+'.eps',transparent=False)
    
    ###############################################################

def plotKct2rel(temparr,Kcte,Kct2,area0,beta,fsuff):
    ##############################################################
    ##plot colour correction factors
    plot.figure(23)
    plot.clf()
    
    #lpsw,=plot.plot(alpharr,Kce[:,0]/Kc2[:,0],'b-')
    #lpmw,=plot.plot(alpharr,Kce[:,1]/Kc2[:,1],'g-')
    #lplw,=plot.plot(alpharr,Kce[:,2]/Kc2[:,2],'r-')

    #lext,=plot.plot(alpharr,Kce[:,0]/Kc2[:,0],'k:',label=r'$K_\mathrm{c,ext}/K_\mathrm{c2}^\mathrm{bm}$')    
    plot.plot(temparr,(Kcte[:,0]/area0[0])/Kct2[:,0],'b-',label='PSW')
    plot.plot(temparr,(Kcte[:,1]/area0[1])/Kct2[:,1],'g-',label='PMW')
    plot.plot(temparr,(Kcte[:,2]/area0[2])/Kct2[:,2],'r-',label='PLW')
        
    #lpt,=plot.plot(alpharr,Kc[:,0]/Kc2[:,0],'k--',label=r'$K_\mathrm{c,pt}/K_\mathrm{c2}^\mathrm{bm}$')
    #plot.plot(alpharr,Kc[:,0]/Kc2[:,0],'b--')
    #plot.plot(alpharr,Kc[:,1]/Kc2[:,1],'g--')
    #plot.plot(alpharr,Kc[:,2]/Kc2[:,2],'r--')
    
    #lpip,=plot.plot(alpharr,Kpip[0]/Kc2[:,0],'k-.',label=r'$K_\mathrm{pip}/K_\mathrm{c2}^\mathrm{bm}$')
    #plot.plot(alpharr,Kpip[0]/Kc2[:,0],'b-.')
    #plot.plot(alpharr,Kpip[1]/Kc2[:,1],'g-.')
    #plot.plot(alpharr,Kpip[2]/Kc2[:,2],'r-.')
    
    plot.xlabel(r'Temperature (K)')
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
    
    plot.savefig('../Outputs/colcorrt2rel'+fsuff+'.png',transparent=False,dpi=300.)
    plot.savefig('../Outputs/colcorrt2rel'+fsuff+'.eps',transparent=False)
    
    ###############################################################
   
if __name__ == "__main__":
    maincode()