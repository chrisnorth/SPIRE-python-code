#from __future__ import division
#General python modules
#import pylab
from numpy import zeros,ones,arange,array,floor,ceil,log,log10,sum,max,min,pi,sqrt
from numpy import concatenate,inf,size,zeros_like
from numpy import polyfit
from scipy.stats import linregress
import sys
#import matplotlib
import matplotlib as mpl

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['legend.fontsize']='x-large'
mpl.rcParams['lines.linewidth']=2
mpl.rcParams['axes.linewidth']=2
mpl.rcParams['axes.labelsize']='x-large'
mpl.rcParams['xtick.labelsize']='x-large'
mpl.rcParams['ytick.labelsize']='x-large'
mpl.rcParams['font.family']='serif'
mpl.rcParams['text.usetex']=True
mpl.rcParams['font.weight']='bold'
mpl.rcParams['font.serif']=['Times', 'Palatino', 'New Century Schoolbook', 'Bookman', 'Computer Modern Roman']
mpl.rcParams['font.sans-serif']=['Helvetica', 'Avant Garde', 'Computer Modern Sans serif']
mpl.rcParams['font.cursive']=['Zapf Chancery']
mpl.rcParams['font.monospace']=['Courier', 'Computer Modern Typewriter']

import matplotlib.pyplot as plot
from scipy import where
#from scipy.interpolate import RectBivariateSpline
import argparse
import os.path
import string
from csv import reader

#Specific modules
from beam import getnewbeamprofile,comb_beam,get_effnu_az_newbeam,beamprof_nu
from beam import arcsec2sr,avgbeam,modbeam_new,modbeam_area,measbeam_new
from beam import beamarea_az,beamarea_az_nu_new,beamarea_az_th_nu_new,beamarea_az_th
from rsrf import getrsrf,getapf,bandedge
from calc_conv import calc_k4,calc_kc,calc_k5

def maincode():
    parser=argparse.ArgumentParser(description='Herschel-SPIRE Extended Emission Calibration')
    parser.add_argument('--brad',action='store',default=600.,type=float,dest='brad',help='Radius to integrate beam out to (arcsec). Default=350')
    parser.add_argument('--ind',action='store',default=0.85,type=float,dest='ind',help='Power law index with which beam FWHM changes with frequency (FWHM \\propto \\nu^{ind}). Default=0.65')
    parser.add_argument('--bzlim',action='store',default=None,type=float,dest='bzlim',help='Zero-limit for beamfiles, below which the values are replaced by BZVAL. Default=None (i.e. don\'t set limit)')
    parser.add_argument('--bzval',action='store',default=0.,type=float,dest='bzval',help='Value with which to replace beam values below the zero-limit (BZLIM). Default=0')
    parser.add_argument('--aunit',action='store',default="arcsec",dest='aunit',help='Unit of area to use for area plots [arcsec|sr]. Default=sr')
    parser.add_argument('--newa', action='store_true',default=False,dest='newa',help='Set to use new area in colour-correction calculations')
    parser.add_argument('--rsrf', action='store',default='M',dest='rsrftype',help='RSRF type to use [M=Measured | T=Top-hat]')
    parser.add_argument('--apeff', action='store',default='R',dest='apftype',help='Aperture Efficiency type to use [R=Real | U=Uniform]')
    parser.add_argument('--beta1',action='store',default=2.,type=float,dest='beta1',help='Greybody 1 spectral index')
    parser.add_argument('--beta2',action='store',default=1.5,type=float,dest='beta2',help='Greybody 2 spectral index')
    parser.add_argument('--afita',action='store',default=1,type=int,dest='afita',help='Degree of polynomial fit to area with spectral index')
    parser.add_argument('--afitt',action='store',default=2,type=int,dest='afitt',help='Degree of polynomial fit to area with temperature')
    parser.add_argument('--nepmod', action='store',default='ESA4',dest='nepmod',help='Model to use for Neptune spectrum ["ESA2"|"ESA4"]')

    parser.add_argument('--prsrf', action='store_true',default=False,dest='prsrf',help='Set to plot RSRF and Ap. Eff.')
    parser.add_argument('--pblim', action='store_true',default=False,dest='pblim',help='Set to plot beam profiles at band limits')
    parser.add_argument('--pbcmb', action='store_true',default=False,dest='pbcmb',help='Set to plot combination of scaled and fixed beams')
    parser.add_argument('--pbcomp', action='store_true',default=False,dest='pbcomp',help='Set to plot beam comparison profiles (measured vs. theoretical)')
    parser.add_argument('--nopbeam', action='store_true',default=False,dest='nopbeam',help='Set to NOT plot beam maps (overrides PBEAM)')
    parser.add_argument('--pbmap', action='store_true',default=False,dest='pbmap',help='Set to plot beam maps (overridden by NOPBEAM)')
    parser.add_argument('--pcorr', action='store_true',default=False,dest='pcorr',help='Set to plot colour correction factors')
    parser.add_argument('--pcort', action='store_true',default=False,dest='pcort',help='Set to plot colour correction factors over temperature')
    parser.add_argument('--pareaf', action='store_true',default=False,dest='pareaf',help='Set to plot beam area against frequency')
    parser.add_argument('--parear', action='store_true',default=False,dest='parear',help='Set to plot beam area against radius')
    parser.add_argument('--pareaa', action='store_true',default=False,dest='pareaa',help='Set to plot beam area against spectral index')
    parser.add_argument('--pareat', action='store_true',default=False,dest='pareat',help='Set to plot beam area against temperature')
    parser.add_argument('--pareath', action='store_true',default=False,dest='pareath',help='Set to plot beam area against source FWHM')
    parser.add_argument('--pfwhm', action='store_true',default=False,dest='pfwhm',help='Set to plot edge-taper and FWHM against wavelength.')
    parser.add_argument('--pbrad', action='store_true',default=False,dest='pbrad',help='Set to plot correction parameters as function of beam radius.')
    parser.add_argument('--puni', action='store_true',default=False,dest='puni',help='Set to plot correction parameters for uniform bands')
    parser.add_argument('--pall', action='store_true',default=False,dest='pall',help='Plot all plots (set NOPBEAM to exclude beam maps)')
    parser.add_argument('--newext', action='store_true',default=False,dest='newext',help='Set to use new extended source method (pipeline in MJy/sr)')

    
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
    beta1=args.beta1
    beta2=args.beta2
    newext=args.newext
    afita=args.afita
    afitt=args.afitt
    nepmod=args.nepmod
    
    if pall:
        prsrf=True
        pblim=True
        pbcmb=True
        pbcomp=True
        pcorr=True
        pcort=True
        pareaf=True
        parear=True
        pareaa=True
        pareat=True
        pareath=True
        pfwhm=True
        pbrad=True
        puni=True
        pbmap=True
    else:
        prsrf=args.prsrf
        pblim=args.pblim
        pbcmb=args.pbcmb
        pcorr=args.pcorr
        pcort=args.pcort
        pareaf=args.pareaf
        parear=args.parear
        pareaa=args.pareaa
        pareat=args.pareat
        pbcomp=args.pbcomp
        pareath=args.pareath
        pfwhm=args.pfwhm
        pbrad=args.pbrad
        puni=args.puni
        pbmap=args.pbmap
    
    if args.nopbeam == True:
        pbmap==False
        
    plotarr=array([prsrf,pblim,pcorr,pcort,pareaf,parear,pareaa,pareat,pbcomp,pareath,pfwhm,pbrad,puni,pbmap])
    
    ##########################
    ### Set figure numbers and labels
    ##########################
    figAntEdge='1a'
    figAntAng='1b'
    figAntFWHM='2a'
    figAntArea='2b'
    figAntApEff='3'
    figAbsAreaApEff='4'
    figuniAntVarR='5a'
    figuniAntVarNu0='5b'
    figuniAntExt='5c'
    figuniAbsVarR='6a'
    figuniAbsVarNu0='6b'
    figuniAbsExt='6c'
    figRSRF='7'
    figKc='8a'
    figKct='8b'
    figbmappsw='9a'
    figbmappmw='9b'
    figbmapplw='9c'
    figbcomb='10'
    figblim='11'
    figbeamcomp='12'
    figareaa='13a'
    figareat='13b'
    figKc2='14a'
    figKct2='14b'
    figKc3va='15'
    figKc3tva='16'
    figKc3tva_zoom=''
    figK4='-'
    figK6va='-'
    figKcall='-'
    figKcKc3tsmall='-'
    figKcKc3tbig='-'
    figSrcArea='-'    
    figKc3='-'
    figKcr='-'
    figKcKc2='-'
    figKctKct2='-'
    figarear='-'

    if newext:
        figKc2rel=figKc2 #only of newext set
        figKct2rel=figKct2 #only of newext set
        figKc3varel=figKc3va #only of newext set
        figKc3tvarel=figKc3tva #only of newext set
        figKc3tvarel_zoom=figKc3tva_zoom #only of newext set
        figKcKc3tbigrel=figKcKc3tbig #only of newext set
        figKcKc2rel=figKcKc2 #only of newext set
        figKctKct2rel=figKctKct2 #only of newext set
    
    labAntEdge='a'
    labAntAng='b'
    labAntFWHM='a'
    labAntArea='b'
    labAntApEff=''
    labAbsAreaApEff=''
    labuniAntVarR='a'
    labuniAntVarNu0='b'
    labuniAntExt='c'
    labuniAbsVarR='a'
    labuniAbsVarNu0='b'
    labuniAbsExt='c'
    labRSRF=''
    labKc='a'
    labKct='b'
    labbmappsw='a'
    labbmappmw='b'
    labbmapplw='c'
    labbcomb=''
    labblim=''
    labbeamcomp=''
    labareaa='a'
    labareat='b'
    labKc2='a'
    labKct2='b'
    labKc3va=''
    labKc3tva=''
    labKc3tva_zoom=''
    labK4=''
    labK6va=''
    labKcall=''
    labKcKc3tsmall=''
    labKcKc3tbig=''
    labSrcArea=''    
    labKc3=''
    labKcr=''
    labKcKc2=''
    labKctKct2=''
    labarear=''

    if newext:
        labKc2rel=labKc2 #only of newext set
        labKct2rel=labKct2 #only of newext set
        labKc3varel=labKc3va #only of newext set
        labKc3tvarel=labKc3tva #only of newext set
        labKc3tvarel_zoom=labKc3tva_zoom #only of newext set
        labKcKc3tbigrel=labKcKc3tbig #only of newext set
        labKcKc2rel=labKcKc2 #only of newext set
        labKctKct2rel=labKctKct2 #only of newext set

        
    ### Plot figures for uniform bands
    if puni:

        fileuniAnt='../Outputs/Kc_Uniform_Antenna.csv'
        funiAnt=open(fileuniAnt,'r').readlines()
        nAlpha=len(funiAnt)-2
        nRes=len(string.split(funiAnt[0],','))-1
        resList=zeros(nRes)
        nu0List=zeros(nRes)
        alphaList=zeros(nAlpha)
        KcUniAnt=zeros((nRes,nAlpha))
        for r in range(nRes):
            resList[r]=float(string.split(funiAnt[0],',')[r+1])
            nu0List[r]=float(string.split(funiAnt[1],',')[r+1])
            for a in range(nAlpha):
                alphaList[a]=float(string.split(funiAnt[a+2],',')[0])
                KcUniAnt[r,a]=float(string.split(funiAnt[a+2],',')[r+1])
        
        resPlot=[0,3,4]
        nu0Plot=[0,1,2]
        plotuniant(resList,nu0List,alphaList,KcUniAnt,resPlot,nu0Plot,figuniAntVarR,figuniAntVarNu0,labuniAntVarR,labuniAntVarNu0)

        fileuniAbs='../Outputs/Kc_Uniform_Absorber.csv'
        funiAbs=open(fileuniAbs,'r').readlines()
        nAlpha=len(funiAbs)-2
        nRes=len(string.split(funiAbs[0],','))-1
        resList=zeros(nRes)
        nu0List=zeros(nRes)
        alphaList=zeros(nAlpha)
        KcUniAbs=zeros((nRes,nAlpha))
        for r in range(nRes):
            resList[r]=float(string.split(funiAbs[0],',')[r+1])
            nu0List[r]=float(string.split(funiAbs[1],',')[r+1])
            for a in range(nAlpha):
                alphaList[a]=float(string.split(funiAbs[a+2],',')[0])
                KcUniAbs[r,a]=float(string.split(funiAbs[a+2],',')[r+1])
        
        resPlot=[0,3,4,5]
        nu0Plot=[0,1,2]
        plotuniabs(resList,nu0List,alphaList,KcUniAbs,resPlot,nu0Plot,figuniAbsVarR,figuniAbsVarNu0,labuniAbsVarR,labuniAbsVarNu0)
    
        fileuniAntExt='../Outputs/Kc2_Uniform_Antenna.csv'
        funiAntExt=open(fileuniAntExt,'r').readlines()
        nAlpha=len(funiAntExt)-2
        nRes=len(string.split(funiAntExt[0],','))-1
        resList=zeros(nRes)
        nu0List=zeros(nRes)
        alphaList=zeros(nAlpha)
        KcUniAntExt=zeros((nRes,nAlpha))
        for r in range(nRes):
            resList[r]=float(string.split(funiAntExt[0],',')[r+1])
            nu0List[r]=float(string.split(funiAntExt[1],',')[r+1])
            for a in range(nAlpha):
                alphaList[a]=float(string.split(funiAntExt[a+2],',')[0])
                KcUniAntExt[r,a]=float(string.split(funiAntExt[a+2],',')[r+1])
        
        resPlot=[0]
        nu0Plot=[0]
        plotuniantext(resList,nu0List,alphaList,KcUniAntExt,resPlot,nu0Plot,figuniAntExt,labuniAntExt)
        
        fileuniAbsExt='../Outputs/Kc2_Uniform_Absorber.csv'
        funiAbsExt=open(fileuniAbsExt,'r').readlines()
        nAlpha=len(funiAbsExt)-2
        nRes=len(string.split(funiAbsExt[0],','))-1
        resList=zeros(nRes)
        nu0List=zeros(nRes)
        alphaList=zeros(nAlpha)
        KcUniAbsExt=zeros((nRes,nAlpha))
        for r in range(nRes):
            resList[r]=float(string.split(funiAbsExt[0],',')[r+1])
            nu0List[r]=float(string.split(funiAbsExt[1],',')[r+1])
            for a in range(nAlpha):
                alphaList[a]=float(string.split(funiAbsExt[a+2],',')[0])
                KcUniAbsExt[r,a]=float(string.split(funiAbsExt[a+2],',')[r+1])
        
        resPlot=[0]
        nu0Plot=[0]
        plotuniabsext(resList,nu0List,alphaList,KcUniAbsExt,resPlot,nu0Plot,figuniAbsExt,labuniAbsExt)
        
    ###Plot Matt's figures    
    if pfwhm:
        file1a='../Inputs/freq_edge-taper.csv'
        f1a=reader(open(file1a))
        nu1a=[]
        et1a=[]
        for row in f1a:
            if string.find(row[0],'#') < 0:
                nu1a.append(float(row[0]))
                et1a.append(float(row[1]))
        nu1a=array(nu1a)
        et1a=array(et1a)
        
        file1b='../Inputs/ang_beam-response.csv'
        f1b=reader(open(file1b))
        ang1b_1=[]
        beam1b_1=[]
        ang1b_2=[]
        beam1b_2=[]
        ang1b_3=[]
        beam1b_3=[]
        for row in f1b:
            if string.find(row[0],'#') < 0:
                ang1b_1.append(float(row[0]))
                beam1b_1.append(float(row[1]))
                ang1b_2.append(float(row[2]))
                beam1b_2.append(float(row[3]))
                ang1b_3.append(float(row[4]))
                beam1b_3.append(float(row[5]))
        ang1b_1=array(ang1b_1)
        beam1b_1=array(beam1b_1)
        ang1b_2=array(ang1b_2)
        beam1b_2=array(beam1b_2)        
        ang1b_3=array(ang1b_3)
        beam1b_3=array(beam1b_3)
        
        plot1(nu1a,et1a,ang1b_1,beam1b_1,ang1b_2,beam1b_2,ang1b_3,beam1b_3,figAntEdge,figAntAng,labAntEdge,labAntAng)
        
        file2a='../Inputs/freq_fwhm.csv'
        f2a=reader(open(file2a))
        nu2a=[]
        fwhm2a_1=[]
        fwhm2a_2=[]
        for row in f2a:
            if string.find(row[0],'#') < 0:
                nu2a.append(float(row[0]))
                fwhm2a_1.append(float(row[1]))
                fwhm2a_2.append(float(row[2]))
        nu2a=array(nu2a)
        fwhm2a_1=array(fwhm2a_1)
        fwhm2a_2=array(fwhm2a_2)
        
        file2b='../Inputs/freq_beam-area.csv'
        f2b=reader(open(file2b))
        nu2b=[]
        area2b_1=[]
        area2b_2=[]
        for row in f2b:
            if string.find(row[0],'#') < 0:
                nu2b.append(float(row[0]))
                area2b_1.append(float(row[1]))
                area2b_2.append(float(row[2]))
        nu2b=array(nu2b)
        area2b_1=array(area2b_1)
        area2b_2=array(area2b_2)
        
        plot2(nu2a,fwhm2a_1,fwhm2a_2,nu2b,area2b_1,area2b_2,figAntFWHM,figAntArea,labAntFWHM,labAntArea)
        
        file3='../Inputs/freq_ap-eff.csv'
        f3=reader(open(file3))
        nu3=[]
        ap3=[]
        for row in f3:
            if string.find(row[0],'#') < 0:
                nu3.append(float(row[0]))
                ap3.append(float(row[1]))
        nu3=array(nu3)
        ap3=array(ap3)
        
        plot3(nu3,ap3,figAntApEff,labAntApEff)
        
        file4='../Inputs/freq_beam-area_ap-eff.csv'
        f4=reader(open(file4))
        nu4=[]
        area4=[]
        ap4=[]
        for row in f4:
            if string.find(row[0],'#') < 0:
                nu4.append(float(row[0]))
                area4.append(float(row[1]))
                ap4.append(float(row[2]))
        nu4=array(nu4)
        area4=array(area4)
        ap4=array(ap4)

        plot4(nu4,area4,ap4,figAbsAreaApEff,labAbsAreaApEff)
        
####    ##Plot FWHM and edge-taper figures
####    if pfwhm:
####        ##Read in Fig 1a (edge taper vs wavelength)
####        file1a='../Inputs/edge-taper_wavelength.csv'
####        f1a=reader(open(file1a))
####        wl1a=[]
####        et1a=[]
####        for row in f1a:
####            if string.find(row[0],'#') < 0:
####                wl1a.append(float(row[0]))
####                et1a.append(float(row[1]))
####        wl1a=array(wl1a)
####        et1a=array(et1a)
####        
####        ##Read in Fig 1b (FWHM vs edge taper)
####        file1b='../Inputs/fwhm_edge-taper.csv'
####        f1b=reader(open(file1b))
####        et1b=[]
####        fw1b=[]
####        for row in f1b:
####            if string.find(row[0],'#') < 0:
####                et1b.append(float(row[0]))
####                fw1b.append(float(row[1]))
####        et1b=array(et1b)
####        fw1b=array(fw1b)
####
####        plot1(wl1a,et1a,et1b,fw1b)
####        
####        file2ff='../Inputs/far-field-fwhm_wavelength.csv'
####        f2ff=reader(open(file2ff))
####        wl2ff=[]
####        fw2ff=[]
####        for row in f2ff:
####            if string.find(row[0],'#') < 0:
####                wl2ff.append(float(row[0]))
####                fw2ff.append(float(row[1]))
####        wl2ff=array(wl2ff)
####        fw2ff=array(fw2ff)
####        
####        file2cet='../Inputs/const-edge-taper-fwhm_wavelength.csv'
####        f2cet=reader(open(file2cet))
####        wl2cet=[]
####        fw2cet=[]
####        for row in f2cet:
####            if string.find(row[0],'#') < 0:
####                wl2cet.append(float(row[0]))
####                fw2cet.append(float(row[1]))
####        wl2cet=array(wl2cet)
####        fw2cet=array(fw2cet)
####        
####        file2b085='../Inputs/fwhm_b0.85_wavelength.csv'
####        f2b085=reader(open(file2b085))
####        wl2b085=[]
####        fw2b085=[]
####        for row in f2b085:
####            if string.find(row[0],'#') < 0:
####                wl2b085.append(float(row[0]))
####                fw2b085.append(float(row[1]))
####        wl2b085=array(wl2b085)
####        fw2b085=array(fw2b085)
####        
####        file2b078='../Inputs/fwhm_b0.78_wavelength.csv'
####        f2b078=reader(open(file2b078))
####        wl2b078=[]
####        fw2b078=[]
####        for row in f2b078:
####            if string.find(row[0],'#') < 0:
####                wl2b078.append(float(row[0]))
####                fw2b078.append(float(row[1]))
####        wl2b078=array(wl2b078)
####        fw2b078=array(fw2b078)
####
####        plot2(wl2ff,fw2ff,wl2cet,fw2cet,wl2b085,fw2b085,wl2b078,fw2b078)
####
####        file3='../Inputs/fwhm_pixel-size.csv'
####        f3=reader(open(file3))
####        px3=[]
####        fw3=[]
####        for row in f3:
####            if string.find(row[0],'#') < 0:
####                px3.append(float(row[0]))
####                fw3.append(float(row[1]))
####        px3=array(px3)
####        fw3=array(fw3)
####    
####        plot3(px3,fw3)

    if nepmod=="ESA2":
        alpha_nep=array((1.26,1.39,1.45))
    elif nepmod=="ESA4":
        alpha_nep=array((1.29,1.42,1.47))
    elif nepmod=="":
        alpha_nep=array((1.29,1.42,1.47))
    else:
        print 'Invalid value for NEPMOD: "%s". Must be "ESA2" or "ESA4"'
        sys.exit()
            
    fsuff='_theta_final_newBprofBS_Br%d_Ind%.2f'%(int(brad),ind)
    if nepmod != '':
        fsuff=fsuff+'_%s'%(nepmod)
    if bzlim != None:
        fsuff=fsuff+'_Bl%1g_Bv%1g'%(bzlim,bzval)
    if rsrftype == 'T':
        fsuff=fsuff+'_noRSRF'
    if apftype == 'U':
        fsuff=fsuff+'_noApEff'
    if newa == True:
        fsuff=fsuff+'_newArea'
    print 'File_suffix= '+fsuff

    ##########################
    ## Read in summary file ##
    ##########################
    
    ##Check whether summary file exists
    file_summ='../Outputs/Summ'+fsuff+'.dat'
    exists=os.path.isfile(file_summ)
    if not exists:
        print 'Files do not exist. Please run extcal_theta_newbeam with the same options.'
        sys.exit()
    
    fsumm=open(file_summ,'r')
    l_sum=fsumm.readlines()
    nl=size(l_sum)
    fsumm.close()
    fsumm=open(file_summ,'r')
    #print l_sum
    nulim=zeros((3,2))
    for l in range(nl):
        line=fsumm.readline()
        #print line
        s0=string.find(line,': ')
        #print string.find(line,'Date')
        if string.find(line,'Date') >=0:
            datestr=line[s0+2:-2]
        if string.find(line,'Beam Areas') >=0:
            area_str=string.split(line[s0+2:-2],',')
            areas=array([float(area_str[0]),float(area_str[1]),float(area_str[2])])
        if string.find(line,'Effective frequencies') >=0:
            nueff_str=string.split(line[s0+2:-2],',')
            nueff=array([float(nueff_str[0])*1.e12,float(nueff_str[1])*1.e12,float(nueff_str[2])*1.e12])
        if string.find(line,'Bandedge (PSW)') >=0:
            nulim[0,0]=float(string.split(line[s0+2:-2],',')[0])*1.e12
            nulim[0,1]=float(string.split(line[s0+2:-2],',')[1])*1.e12
        if string.find(line,'Bandedge (PMW)') >=0:
            nulim[1,0]=float(string.split(line[s0+2:-2],',')[0])*1.e12
            nulim[1,1]=float(string.split(line[s0+2:-2],',')[1])*1.e12
        if string.find(line,'Bandedge (PLW)') >=0:
            nulim[2,0]=float(string.split(line[s0+2:-2],',')[0])*1.e12
            nulim[2,1]=float(string.split(line[s0+2:-2],',')[1])*1.e12
        if string.find(line,'Beam freq-dependent power-law') >= 0:
            ind=float(line[s0+2:])
    print 'Date:',datestr
    
    #print 'nu_eff:',nueff
    #print 'nulim:',nulim
    fsumm.close()
        
    c=299792458. #speed of light

    #set band centres
    wlc=(250.e-6,350.e-6,500.e-6)
    wlc=array(wlc)

    nuc=c/wlc
    #print nuc/1.e12
    band=arange(3)
    alpharr=None
    #########################
    ##  Source properties  ##
    #########################
    
    #Set source brightness and spectral index
    
    area0=arcsec2sr(array([426.,771.,1626.])) #area in steradians
    
    ####################
    ##  Read in RSRF  ##
    ####################
    
    print 'Reading in RSRF...'
    file='../Outputs/RSRF'+fsuff+'.csv'
    rsrf_input=reader(open(file))
    nulist=[]
    rsrf_psw=[]
    rsrf_pmw=[]
    rsrf_plw=[]
    for row in rsrf_input:
        if string.find(row[0],'#') < 0:
            nulist.append(float(row[0])*1.e12)
            rsrf_psw.append(float(row[1]))
            rsrf_pmw.append(float(row[2]))
            rsrf_plw.append(float(row[3]))

    nuarr=array(nulist)
    nnu=size(nuarr)
    rsrfarr=zeros((nnu,3))
    rsrfarr[:,0]=rsrf_psw
    rsrfarr[:,1]=rsrf_pmw
    rsrfarr[:,2]=rsrf_plw

    ilim=zeros((3,2))
    ilim_2=zeros((3,2))
    nulim_2=zeros((3,2))
    nulim_3=zeros((3,2))
    for b in range(3):
        (ilim[b,:],nulim_2[b,:])=bandedge(nuarr,rsrfarr[:,b])
        (ilim_2[b,:],nulim_3[b,:])=bandedge(nuarr,rsrfarr[:,b],fact=2.)
    
    if prsrf:

        ####################
        ##  Read in APF   ##
        ####################

        print 'Reading in Ap.Eff...'    
        file='../Outputs/ApEff'+fsuff+'.csv'
        apf_input=reader(open(file))
        apf_psw=[]
        apf_pmw=[]
        apf_plw=[]
        for row in apf_input:
            if string.find(row[0],'#') < 0:
                apf_psw.append(float(row[1]))
                apf_pmw.append(float(row[2]))
                apf_plw.append(float(row[3]))
        
        apfarr=zeros((nnu,3))
        apfarr[:,0]=apf_psw
        apfarr[:,1]=apf_pmw
        apfarr[:,2]=apf_plw
        
        #Compute 'extended source' rsrf
        rsrfearr=zeros((nnu,3))
        for b in band:
            rsrfearr[:,b]=rsrfarr[:,b]*(nuc[b]/nuarr)**2
        
        
        plotrsrf(nuarr,rsrfarr,rsrfearr,apfarr,fsuff,figRSRF,labRSRF)

    if pcorr:
        
        ########################################
        ##  Read in K4,K4e,Kc,Kce,K5,Kc2,Kc2b ##
        ########################################
        
        print 'Reading in K4...'
        
        ## Read in K4    
        
        file='../Outputs/K4'+fsuff+'.csv'
        k4_input=reader(open(file))
        alist=[]
        k4_psw=[]
        k4_pmw=[]
        k4_plw=[]
        Kpip=zeros(3)
        for row in k4_input:
            if string.find(row[0],'#') < 0:
                alist.append(float(row[0]))
                k4_psw.append(float(row[1]))
                k4_pmw.append(float(row[2]))
                k4_plw.append(float(row[3]))
                if float(row[0])==-1:
                    Kpip[0]=float(row[1])
                    Kpip[1]=float(row[2])
                    Kpip[2]=float(row[3])
    
        alpharr=array(alist)
        nalph=size(alpharr)
        K4=zeros((nalph,3))
        K4[:,0]=k4_psw
        K4[:,1]=k4_pmw
        K4[:,2]=k4_plw
    
        ## Read in K4e
        print 'Reading in K4e...'
        K4e=zeros((nalph,3))
        r=0
        file='../Outputs/K4e'+fsuff+'.csv'
        k4e_input=reader(open(file))
        for row in k4e_input:
            if string.find(row[0],'#') < 0:
                K4e[r,0]=float(row[1])
                K4e[r,1]=float(row[2])
                K4e[r,2]=float(row[3])
                r=r+1
       
        ## Read in Kc    
        print 'Reading in Kc...'
        Kc=zeros((nalph,3))
        r=0
        file='../Outputs/Kc'+fsuff+'.csv'
        kc_input=reader(open(file))
        for row in kc_input:
            if string.find(row[0],'#') < 0:
                Kc[r,0]=float(row[1])
                Kc[r,1]=float(row[2])
                Kc[r,2]=float(row[3])
                r=r+1
        
        ## Read in Kce    
        print 'Reading in Kce...'
        Kce=zeros((nalph,3))
        r=0
        file='../Outputs/Kce'+fsuff+'.csv'
        kce_input=reader(open(file))
        for row in kce_input:
            if string.find(row[0],'#') < 0:
                Kce[r,0]=float(row[1])
                Kce[r,1]=float(row[2])
                Kce[r,2]=float(row[3])
                r=r+1
    
        ## Read in K5    
        print 'Reading in K5...'
        K5=zeros((nalph,3))
        r=0
        file='../Outputs/K5'+fsuff+'.csv'
        k5_input=reader(open(file))
        for row in k5_input:
            if string.find(row[0],'#') < 0:
                K5[r,0]=float(row[1])
                K5[r,1]=float(row[2])
                K5[r,2]=float(row[3])
                r=r+1
    
        ## Read in Kc2
        print 'Reading in Kc2...'
        Kc2=zeros((nalph,3))
        Kc2rel=zeros((nalph,3))
        KpipE=zeros(3)
        r=0
        file='../Outputs/Kc2'+fsuff+'.csv'
        kc2_input=reader(open(file))
        for row in kc2_input:
            if string.find(row[0],'#') < 0:
                Kc2[r,0]=float(row[1])
                Kc2[r,1]=float(row[2])
                Kc2[r,2]=float(row[3])
                if float(row[0])==-1:
                    KpipE[0]=float(row[1])
                    KpipE[1]=float(row[2])
                    KpipE[2]=float(row[3])
                r=r+1
    
        for b in band:
            Kc2rel[:,b]=Kc2[:,b]/KpipE[b]
        
        ## Read in Kc2b
        print 'Reading in Kc2b...'
        Kc2b=zeros((nalph,3))
        Kc2brel=zeros((nalph,3))
        r=0
        file='../Outputs/Kc2b'+fsuff+'.csv'
        kc2b_input=reader(open(file))
        for row in kc2b_input:
            if string.find(row[0],'#') < 0:
                Kc2b[r,0]=float(row[1])
                Kc2b[r,1]=float(row[2])
                Kc2b[r,2]=float(row[3])
                r=r+1

        for b in band:
            Kc2brel[:,b]=Kc2b[:,b]/KpipE[b]

        plotK4(alpharr,K4,K4e,fsuff,figK4,labK4)
        plotKc(alpharr,Kc,Kce,fsuff,figKc,labKc)

        ##################################################
        ## K4: Jy/beam per Jy_meas for point source [f(alpha)]
        ## Kpip: Jy/beam per Jy_meas for point source with alpha=-1
        ## Kc: Jy/beam per Jy_pip for point source [f(alpha)]
        ## K4e: Jy/beam per Jy_meas for fully-extended source (old method) [f(alpha)]
        ## Kce: Jy/beam per Jy_pip for fully-extended source (old method) [f(alpha)]
        ## K5: Jy/sr per Jy_meas for fully-extended source [f(alpha)]
        ## Kc2: Jy/sr per Jy_pip for fully-extended source [f(alpha)]
        ## KpipE: Jy/sr per Jy_pip for fully extended source with alpha=-1
        ## Kc2rel=Kc2/KpipE [f(alpha)]
        ## Kc2b: Jy/beam per Jy_pip for fully-extended source [f(alpha)]
        ## Kc2brel: Kc2b/KpipE [f(alpha)]
        ##################################################


    if pcorr:
        #plotK5(arr_a,arr_t,K5,fsuff)
        #plotKc2(arr_a,arr_t,Kc2,fsuff)

        ###########################
        ## Read in K6, Kc3, Kc3t ##
        ###########################
        file='../Outputs/K6_PSW'+fsuff+'.csv'
        k6psw_input=reader(open(file))
        thlist=[]
        for row in k6psw_input:
            if string.find(row[0],'#') < 0:
                thlist.append(float(row[0]))
        thetarr=array(thlist)
        nth=thetarr.size
        
        K6=zeros((nalph,nth,3))
        Kc3=zeros((nalph,nth,3))
        Kc3t=zeros((nalph,nth,3))
        Kc3rel=zeros((nalph,nth,3))
        Kc3trel=zeros((nalph,nth,3))
        
        print 'Reading in K6...'
        r=0
        for row in k6psw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    K6[a,r,0]=float(row[a+1])
                r=r+1
        
        file='../Outputs/K6_PMW'+fsuff+'.csv'
        k6pmw_input=reader(open(file))
        r=0
        for row in k6pmw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    K6[a,r,1]=float(row[a+1])
                r=r+1
        
        file='../Outputs/K6_PLW'+fsuff+'.csv'
        k6plw_input=reader(open(file))
        r=0
        for row in k6plw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    K6[a,r,2]=float(row[a+1])
                r=r+1
        
        ## Read in Kc3
        print 'Reading in Kc3...'
        
        file='../Outputs/Kc3_PSW'+fsuff+'.csv'
        kc3psw_input=reader(open(file))
        r=0
        for row in kc3psw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    Kc3[a,r,0]=float(row[a+1])
                Kc3rel[:,r,0]=Kc3[:,r,0]/KpipE[0]
                r=r+1
        
        file='../Outputs/Kc3_PMW'+fsuff+'.csv'
        kc3pmw_input=reader(open(file))
        r=0
        for row in kc3pmw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    Kc3[a,r,1]=float(row[a+1])
                Kc3rel[:,r,1]=Kc3[:,r,1]/KpipE[1]
                r=r+1
        
        file='../Outputs/Kc3_PLW'+fsuff+'.csv'
        kc3plw_input=reader(open(file))
        r=0
        for row in kc3plw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    Kc3[a,r,2]=float(row[a+1])
                Kc3rel[:,r,2]=Kc3[:,r,2]/KpipE[2]
                r=r+1
        
        ## Read in Kc3t
        print 'Reading in Kc3t...'
        
        file='../Outputs/Kc3t_PSW'+fsuff+'.csv'
        kc3tpsw_input=reader(open(file))
        r=0
        for row in kc3tpsw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    Kc3t[a,r,0]=float(row[a+1])
                Kc3trel[:,r,0]=Kc3t[:,r,0]/KpipE[0]
                r=r+1
        
        file='../Outputs/Kc3t_PMW'+fsuff+'.csv'
        kc3tpmw_input=reader(open(file))
        r=0
        for row in kc3tpmw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    Kc3t[a,r,1]=float(row[a+1])
                Kc3trel[:,r,1]=Kc3t[:,r,1]/KpipE[1]
                r=r+1
        
        file='../Outputs/Kc3t_PLW'+fsuff+'.csv'
        kc3tplw_input=reader(open(file))
        r=0
        for row in kc3tplw_input:
            if string.find(row[0],'#') < 0:
                for a in range(nalph):
                    Kc3t[a,r,2]=float(row[a+1])
                Kc3trel[:,r,2]=Kc3t[:,r,2]/KpipE[2]
                r=r+1
                    
        ##################################################
        ## K6: Jy/sr per Jy_meas for partial source [f(alpha,theta)]
        ## Kc3: Jy/sr per Jy_pip for partial source [f(alpha,theta)]
        ## Kc3t: Jy per Jy_pip for partial source [f(alpha,theta)]
        ## Kc3rel: Kc3/KpipE [f(alpha,theta)]
        ## Kc3trel: Kc3t/KpipE [f(alpha,theta)]
        ##################################################

        aplot=14 #alpha=3
        #aplot=6 #alpha=-1
        iths=0
        thsmall=thetarr[iths]
        ithb=nth-1
        thbig=thetarr[ithb]
        

        
        plotK6_va(alpharr,thetarr,K6,K5,aplot,fsuff,figK6va,labK6va)
        if newext:
            plotKc3_varel(alpharr,thetarr,Kc3rel,Kc2rel,aplot,fsuff,figKc3varel,labKc3varel)
            plotKc3t_varel(alpharr,thetarr,Kc3trel,Kc,aplot,fsuff,figKc3tvarel,figKc3tvarel_zoom,labKc3tvarel,labKc3tvarel_zoom)
            #plotKc_allrel(alpharr,Kc,Kce,Kc3rel[:,ithb,:],Kc3trel[:,0,:],fsuff,figKcallrel,labKcallrel)
            #plotKcKc3t_thsmallrel(alpharr,Kc,Kc3trel[:,iths,:],thsmall,fsuff,figKcKc3tsmallrel,labKcKc3tsmallrel)
            plotKc2Kc3_thbigrel(alpharr,Kc2rel,Kc3rel[:,ithb,:],thbig,fsuff,figKcKc3tbigrel,labKcKc3tbigrel)
            plotKc2rel(alpharr,Kce,Kc2rel,area0,fsuff,figKc2rel,labKc2rel)
            plotKcKc2rel(alpharr,Kc,Kc2rel,fsuff,figKcKc2rel,labKcKc2rel)
        else:
            plotKc3_va(alpharr,thetarr,Kc3,Kc2,aplot,fsuff,figKc3va,labKc3va)
            plotKc3t_va(alpharr,thetarr,Kc3t,Kc,aplot,fsuff,figKc3tva,figKc3tva_zoom,labKc3tva,labKc3tva_zoom)
            plotKc_all(alpharr,Kc,Kce,Kc3[:,ithb,:],Kc3t[:,0,:],fsuff,figKcall,labKcall)
            plotKcKc3t_thsmall(alpharr,Kc,Kc3t[:,iths,:],thsmall,fsuff,figKcKc3tsmall,labKcKc3tsmall)
            plotKc2Kc3_thbig(alpharr,Kc2,Kc3[:,ithb,:],thbig,fsuff,figKcKc3tbig,labKcKc3tbig)
            #plotsrcarea(thetarr,srcarea1,srcarea2,fsuff,figSrcArea)
            #plotKc3(arr_a,arr_t,Kc3,fsuff,figKc3)        
            plotKc2(alpharr,Kce,Kc2,area0,fsuff,figKc2,labKc2)
            plotKcKc2(alpharr,Kc,Kc2,fsuff,figKcKc2,labKcKc2)
 
    ###################################
    ## Read in corrections over range of temperatures
    ##  for two different beta
    ###################################
 
    if pcort or pareat:
        fsufft1='_theta_newBprofBS_temp_beta%.1f_Br%d_Ind%.2f'%(beta1,int(brad),ind)
        fsufft2='_theta_newBprofBS_temp_beta%.1f_Br%d_Ind%.2f'%(beta2,int(brad),ind)
        if nepmod!="":
            fsufft1=fsufft1+'_%s'%(nepmod)
            fsufft2=fsufft2+'_%s'%(nepmod)
        if bzlim != None:
            fsufft1=fsufft1+'_Bl%1g_Bv%1g'%(bzlim,bzval)
            fsufft2=fsufft2+'_Bl%1g_Bv%1g'%(bzlim,bzval)
        if rsrftype == 'T':
            fsufft1=fsufft1+'_noRSRF'
            fsufft2=fsufft2+'_noRSRF'
        if apftype == 'U':
            fsufft1=fsufft1+'_noApEff'
            fsufft2=fsufft2+'_noApEff'
        if newa == True:
            fsufft1=fsufft1+'_newArea'
            fsufft2=fsufft2+'_newArea'
    
    if pcort:
        ## Read in Kct1
        tlist1=[]
        Kct1_psw=[]
        Kct1_pmw=[]
        Kct1_plw=[]
        filet1='../Outputs/Kct'+fsufft1+'.csv'
        kct1_input=reader(open(filet1))
        for row in kct1_input:
            if string.find(row[0],'#') < 0:
                tlist1.append(float(row[0]))
                Kct1_psw.append(float(row[1]))
                Kct1_pmw.append(float(row[2]))
                Kct1_plw.append(float(row[3]))
        temparr1=array(tlist1)
        ntemp1=temparr1.size
        Kct1=zeros((ntemp1,3))
        Kct1[:,0]=Kct1_psw
        Kct1[:,1]=Kct1_pmw
        Kct1[:,2]=Kct1_plw

        ##Read in Kct21       
        Kct21=zeros((ntemp1,3))
        r=0
        filet1='../Outputs/Kct2'+fsufft1+'.csv'
        kct21_input=reader(open(filet1))
        for row in kct21_input:
            if string.find(row[0],'#') < 0:
                Kct21[r,0]=float(row[1])
                Kct21[r,1]=float(row[2])
                Kct21[r,2]=float(row[3])
                r=r+1
        
        ## Read in Kct2
        tlist2=[]
        Kct2_psw=[]
        Kct2_pmw=[]
        Kct2_plw=[]
        filet2='../Outputs/Kct'+fsufft2+'.csv'
        kct2_input=reader(open(filet2))
        for row in kct2_input:
            if string.find(row[0],'#') < 0:
                tlist2.append(float(row[0]))
                Kct2_psw.append(float(row[1]))
                Kct2_pmw.append(float(row[2]))
                Kct2_plw.append(float(row[3]))
        temparr2=array(tlist2)
        ntemp2=temparr2.size
        Kct2=zeros((ntemp2,3))
        Kct2[:,0]=Kct2_psw
        Kct2[:,1]=Kct2_pmw
        Kct2[:,2]=Kct2_plw

        ##Read in Kct22
        Kct22=zeros((ntemp2,3))
        r=0
        filet2='../Outputs/Kct2'+fsufft2+'.csv'
        kct22_input=reader(open(filet2))
        for row in kct22_input:
            if string.find(row[0],'#') < 0:
                Kct22[r,0]=float(row[1])
                Kct22[r,1]=float(row[2])
                Kct22[r,2]=float(row[3])
                r=r+1

        if newext:
            Kct21rel=Kct21/KpipE
            Kct22rel=Kct22/KpipE
        
        plotKct(temparr1,Kct1,beta1,temparr2,Kct2,beta2,fsuff,figKct,labKct)
        if newext:
            plotKct2rel(temparr1,Kct21rel,beta1,temparr2,Kct22rel,beta2,fsuff,figKct2rel,labKct2rel)
            plotKctKct2rel(temparr1,Kct1,Kct21rel,beta1,temparr2,Kct2,Kct22rel,beta2,fsuff,figKctKct2rel,labKctKct2rel)
        else:
            plotKct2(temparr1,Kct21,beta1,temparr2,Kct22,beta2,fsuff,figKct2,labKct2)
            plotKctKct2(temparr1,Kct1,Kct21,beta1,temparr2,Kct2,Kct22,beta2,fsuff,figKctKct2,labKctKct2)
        

    
        ##################################################
        ## Kct1: Jy/beam per Jy_pip for point source with beta=beta1 [f(temp)]
        ## Kct2: Jy/beam per Jy_pip for point source with beta=beta2 [f(temp)]
        ## Kct21: Jy/sr per Jy_pip for fully-extended source with beta=beta1 [f(temp)]
        ## Kct22: Jy/sr per Jy_pip for fully-extended source with beta=beta2 [f(temp)]
        ## Kct21rel: Kct21/KpipE [f(temp)]
        ## Kct22rel: Kct22/KpipE [f(temp)]
        ##################################################
        
    ###########################################
    ## Read in corrections for various radii ##
    ###########################################
    
###    if pbrad:
###        ##Read in many band widths
###        radlist=array([25.,50.,60.,75.,100.,125.,150.,200.,250.,300.,500.,1000.])
###        nrad=len(radlist)
###        alphaplot=3.
###        Kc2rad=zeros((nrad,3))
###        Kc2radrel=zeros((nrad,3))
###
###        bsw=array([75,100,150])
###
###        for r in range(nrad):
###            fsuffr='_theta_newBprofBS_Br%d_Ind%.2f'%(int(radlist[r]),ind)
###            if bzlim != None:
###                fsuffr=fsuffr+'_Bl%1g_Bv%1g'%(bzlim,bzval)
###            if rsrftype == 'T':
###                fsuffr=fsuffr+'_noRSRF'
###            if apftype == 'U':
###                fsuffr=fsuffr+'_noApEff'
###            if newa == True:
###                fsuffr=fsuffr+'_newArea'
###            print 'Input file: '+fsuffr
###                
###            ##Read in Kc2
###            fileKc2='../Outputs/Kc2'+fsuffr+'.csv'
###            fKc2Input=reader(open(fileKc2,'r'))
###            for row in fKc2Input:
###                if string.find(row[0],'#') < 0:
###                    if float(row[0]) == alphaplot:
###                            Kc2rad[r,0]=float(row[1])
###                            Kc2rad[r,1]=float(row[2])
###                            Kc2rad[r,2]=float(row[3])
###
###        for r in range(nrad):
###            for b in range(3):
###                Kc2radrel[r,b]=(Kc2rad[r,b]-Kc2rad[nrad-1,b])/Kc2rad[nrad-1,b]
###                
###        print Kc2rad
###        print Kc2radrel               
###        plotKcr(radlist,alphaplot,Kc2radrel,bsw,fsuff,figKcr)
        
    ####################
    ##  Read in beam  ##
    ####################
    
    #Beamtype options:
    #  G=Gaussian
    #  E=Elliptical Gaussian
    #  M=Measured beams
    #  T=Theoretical (modelled) beams)
    #  S=Spliced beams
    
    if pbmap:
        print 'Plotting PSW beam map'
        fitsInPsw='../Inputs/spire_beams_measured/psw_beam_1arcsec.fits'
        plotprettybeam(figbmappsw,labbmappsw,fitsIn=fitsInPsw,band='psw')

        print 'Plotting PMW beam map'
        fitsInPmw='../Inputs/spire_beams_measured/pmw_beam_1arcsec.fits'
        plotprettybeam(figbmappmw,labbmappmw,fitsIn=fitsInPmw,band='pmw')
        
        print 'Plotting PLW beam map'
        fitsInPlw='../Inputs/spire_beams_measured/plw_beam_1arcsec.fits'
        plotprettybeam(figbmapplw,labbmapplw,fitsIn=fitsInPlw,band='plw')
        
    if pblim or pbcmb or parear or pbcomp or pareaa or pareat:
        print 'Reading in beams...'
        
        radarr=arange(brad)
        nrad=radarr.size
        beam_scl=zeros((nrad,3))
        beam_fix=zeros((nrad,3))
        beam_cmb=zeros((nrad,3))
        bsw=array([[70,76,78,87],[95,103,110,130],[135,145,169,180]])
        for b in band:
            beam_scl[:,b],beam_fix[:,b]=getnewbeamprofile(b,radarr,bzlim=bzlim,bzval=bzval)
            beam_cmb[:,b]=comb_beam(beam_scl[:,b],beam_fix[:,b])
            
        print 'Beam areas: [%.2f,%.2f,%.2f]'%(areas[0],areas[1],areas[2])
        
    if pblim:
        beam_scl_lim=zeros((nrad,3,3)) #rad,nu,band
        beam_cmb_lim=zeros((nrad,3,3)) #rad,nu,band
        for b in band:
            nulist=array([nulim[b,0],nuc[b],nulim[b,1]])
            beam_cmb_lim[:,:,b]=beamprof_nu(radarr,beam_scl[:,b],beam_fix[:,b],nulist,nuc[b],ind=ind)
            beam_scl_lim[:,:,b]=beamprof_nu(radarr,beam_scl[:,b],zeros(nrad),nulist,nuc[b],ind=ind)
        
        plotblim(radarr,beam_scl_lim,beam_cmb_lim,beam_fix,brad,fsuff,figblim,labblim)

    if pbcmb:
        plotbcomb(radarr,beam_scl,beam_fix,beam_cmb,bsw,fsuff,figbcomb,labbcomb)

    if parear or pbcomp:
        print 'Computing beam area profiles...'
        area_rarr=arange(0,brad+10,10)
        nrad_a=area_rarr.size
        areas_r=zeros((nrad_a,3))
        for b in band:
            for r in arange(nrad_a):
                #areas_r[r,b]=beamarea_az(radarr,beam_sm[:,b],brad=area_rarr[r])
                radx=float(area_rarr[r])
                areas_r[r,b]=beamarea_az(radarr,beam_cmb[:,b],brad=radx)
                #print '  Radius %.1f: Area %.2f'%(radx,areas_r[r,b])
            #print min(beam_cmb[:,b]),max(beam_cmb[:,b])
        #print 'Areas: [%.2f , %.2f , %.2f] sq-arcsec' % (areas_r[nrad_a-1,0],areas_r[nrad_a-1,1],areas_r[nrad_a-1,2])
        
        if parear:
            plotarear(area_rarr,areas_r,fsuff,figarear,labarear,unit=aunit)

        #plot.figure(30)
        #plot.clf()
        #plot.plot(radarr,beam_scl[:,0],'b-')
        #plot.plot(radarr,beam_scl[:,1],'g-')
        #plot.plot(radarr,beam_scl[:,2],'r-')
        #plot.plot(radarr,beam_fix[:,0],'b--')
        #plot.plot(radarr,beam_fix[:,1],'g--')
        #plot.plot(radarr,beam_fix[:,2],'r--')
        #plot.yscale('log')

    if pareaa or pareat:
        ## Read in Beam areas
        alphabm=[]
        areaa_psw=[]
        areaa_pmw=[]
        areaa_plw=[]
        areapip=zeros(3)
        fileareaa='../Outputs/beamarea'+fsuff+'.csv'
        areaa_input=reader(open(fileareaa))
        for row in areaa_input:
            if string.find(row[0],'#') < 0:
                alphabm.append(float(row[0]))
                areaa_psw.append(float(row[1]))
                areaa_pmw.append(float(row[2]))
                areaa_plw.append(float(row[3]))
                if float(row[0])==-1:
                    areapip[0]=float(row[1])
                    areapip[1]=float(row[2])
                    areapip[2]=float(row[3])
                    print areapip
        alphabm=array(alphabm)
        nabm=alphabm.size
        area_eff=zeros((nabm,3))
        area_eff[:,0]=areaa_psw
        area_eff[:,1]=areaa_pmw
        area_eff[:,2]=areaa_plw
        
        if pareaa:
            ##Fit beam areas
            print 'Fitting beam areas against spectral index...'
            fitpar=zeros((afita+1,3))
            areaa_fit=zeros((nabm,3))
            for b in band:
                print 'Band %d...'%(b)
                #for a in arange(nalph):
                #    area_eff[a,b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpharr[a],ind=ind)
                ##fit power law
                #fitpar[0,b],fitpar[1,b],rval,pval,stderr=linregress(log10(alphabm-alpha_nep),log10(area_eff[:,b]/areas[b])
                fitpar[:,b]=polyfit(alphabm,area_eff[:,b]/areapip[b],deg=afita)
                for d in range(0,afita+1):
                    apow=afita-d
                    print d,apow,fitpar[d,b]
                    areaa_fit[:,b]=areaa_fit[:,b] + fitpar[d,b]*alphabm[:]**(apow)
                areaa_fit[:,b]=areaa_fit[:,b]*areapip[b]
                print 'max offset (%): ',max(100.*abs(areaa_fit[:,b]-area_eff[:,b])/area_eff[:,b])
                #areaa_fit[:,b]=areapip[b]*(fitpar[0,b]*alphabm**2 + fitpar[1,b]*alphabm + fitpar[2,b])
                
            print 'Pipeline beam areas:',areapip
            plotareaa(alphabm,area_eff,areapip,areaa_fit,fitpar,fsuff,figareaa,labareaa)
        
        if pareat:
            ## Read in Beam areas
            tempbm=[]
            tmin=5.
            tmax=40.
            areat1_psw=[]
            areat1_pmw=[]
            areat1_plw=[]
            fileareat1='../Outputs/beamarea'+fsufft1+'.csv'
            areat1_input=reader(open(fileareat1))
            for row in areat1_input:
                if string.find(row[0],'#') < 0:
                    if float(row[0])>tmin and float(row[0])<tmax:
                        tempbm.append(float(row[0]))
                        areat1_psw.append(float(row[1]))
                        areat1_pmw.append(float(row[2]))
                        areat1_plw.append(float(row[3]))
                    
            tempbm=array(tempbm)
            ntbm=tempbm.size
            area_efft1=zeros((ntbm,3))
            area_efft1[:,0]=areat1_psw
            area_efft1[:,1]=areat1_pmw
            area_efft1[:,2]=areat1_plw
            
            areat2_psw=[]
            areat2_pmw=[]
            areat2_plw=[]
            fileareat2='../Outputs/beamarea'+fsufft2+'.csv'
            #fileareat2=fileareat1
            areat2_input=reader(open(fileareat2))
            for row in areat2_input:
                if string.find(row[0],'#') < 0:
                    if float(row[0])>tmin and float(row[0])<tmax:
                        areat2_psw.append(float(row[1]))
                        areat2_pmw.append(float(row[2]))
                        areat2_plw.append(float(row[3]))
            
            area_efft2=zeros((ntbm,3))
            area_efft2[:,0]=areat2_psw
            area_efft2[:,1]=areat2_pmw
            area_efft2[:,2]=areat2_plw
                        
            ##Fit beam areas
            print 'Fitting beam areas against temperature...'
            fitpart1=zeros((afitt+1,3))
            fitpart2=zeros((afitt+1,3))
            areat1_fit=zeros((ntbm,3))
            areat2_fit=zeros((ntbm,3))
            for b in band:
                print 'Band %d...'%(b)
                #for a in arange(nalph):
                #    area_eff[a,b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpharr[a],ind=ind)
                ##fit power law
                #fitpar[0,b],fitpar[1,b],rval,pval,stderr=linregress(log10(alphabm-alpha_nep),log10(area_eff[:,b]/areas[b])
                fitpart1[:,b]=polyfit(tempbm,area_efft1[:,b]/areapip[b],deg=afitt)
                fitpart2[:,b]=polyfit(tempbm,area_efft2[:,b]/areapip[b],deg=afitt)
                for d in range(0,afitt+1):
                    tpow=afitt-d
                    print d,tpow,fitpart1[d,b]
                    areat1_fit[:,b]=areat1_fit[:,b] + fitpart1[d,b]*tempbm[:]**(tpow)
                    areat2_fit[:,b]=areat2_fit[:,b] + fitpart2[d,b]*tempbm[:]**(tpow)
                areat1_fit[:,b]=areat1_fit[:,b]*areapip[b]
                areat2_fit[:,b]=areat2_fit[:,b]*areapip[b]
                print 'max offset [beta1] (%): ',max(100.*abs(areat1_fit[:,b]-area_efft1[:,b])/area_efft1[:,b])
                print 'max offset [beta2] (%): ',max(100.*abs(areat2_fit[:,b]-area_efft2[:,b])/area_efft2[:,b])
                #areaa_fit[:,b]=areapip[b]*(fitpar[0,b]*alphabm**2 + fitpar[1,b]*alphabm + fitpar[2,b])
                
            plotareat(tempbm,beta1,beta2,area_efft1,area_efft2,areapip,areat1_fit,areat2_fit,fitpart1,fitpart2,fsuff,figareat,labareat)
    
    if pbcomp:
        print 'Comparing modelled and measured areas...'
        #alpha_nep=array((1.26,1.39,1.45))
        beam_mod=zeros((nrad,3))
        beam_comp=zeros((nrad,3))
        #areas_mod=zeros((nrad_a,3))
        area_rarr2=arange(1,101,1)
        nrad_a2=area_rarr2.size
        areas_nep=zeros((nrad_a2,3))
        areas_mod=zeros((nrad_a2,3))
        areas_rel=zeros((nrad_a2,3))
        for b in band:
            beam_mod[:,b]=modbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpha_nep[b],ind=ind)
            beam_comp[:,b]=(beam_mod[:,b]-beam_cmb[:,b])/areas[b]
            areas_nep[:,b]=modbeam_area(radarr,beam_cmb[:,b],area_rarr2)
            areas_mod[:,b]=modbeam_area(radarr,beam_mod[:,b],area_rarr2)
            areas_rel[:,b]=(areas_mod[:,b]-areas_nep[:,b])/areas_nep[:,b]

        plot.figure(31)
        plot.clf()
        plot.plot(radarr,beam_cmb[:,0],'b-')
        plot.plot(radarr,beam_cmb[:,1],'g-')
        plot.plot(radarr,beam_cmb[:,2],'r-')
        plot.plot(radarr,beam_mod[:,0],'b--')
        plot.plot(radarr,beam_mod[:,1],'g--')
        plot.plot(radarr,beam_mod[:,2],'r--')
        plot.yscale('log')
        plotbeamcomp(radarr,beam_cmb,beam_mod,area_rarr2,areas_rel,fsuff,figbeamcomp,labbeamcomp)
       
    #Show plots
    if plotarr.any():
        print 'Displaying plots...'
        plot.show()
    
    
#===============================================================================
#     Plotting routines
#===============================================================================
        
def plotrsrf(nuarr,rsrfarr,rsrfearr,apfarr,fsuff,figRSRF,figlab):
    #########################################################
    nuplt=nuarr/1.e12
    #Plot rsrf and aperture function
    plot.figure(1)
    plot.clf()
    #plot rsrf
    plot.plot(nuplt,rsrfarr[:,0],'k-',label=r'Spectral Response')
    plot.plot(nuplt,rsrfarr[:,1],'k-',label=r'_nolegend_')
    plot.plot(nuplt,rsrfarr[:,2],'k-',label=r'_nolegend_')
    
    #plot nu^2.RSRF
    #plot.plot(nuplt,rsrfearr[:,0],'b--',label=r'RSRF x $\nu^{-2}$')
    #plot.plot(nuplt,rsrfearr[:,1],'g--',label='_nolegend_')
    #plot.plot(nuplt,rsrfearr[:,2],'r--',label='_nolegend_')
    
    #plot aperture efficieny
    plot.plot(nuplt,apfarr[:,0],'k--',label=r'Aperture Efficiency')
    plot.plot(nuplt,apfarr[:,1],'k--',label=r'_nolegend_')
    plot.plot(nuplt,apfarr[:,2],'k--',label=r'_nolegend_')
    
    plot.ylim(0.,1.)
    plot.xlim(0.3,1.7)
    plot.xlabel(r'Frequency [THz]')
    plot.ylabel(r'Spectral Response and Aperture Efficiency [arb. units]')
    #plot.title('RSRF and Aperture Efficiency')
    plot.legend(frameon=False)

    plot.box=False
    plot.annotate(r'250$\,\mu$m',(1.250,0.05),ha='center',color='k',weight='bold',size='x-large',family='sans-serif',va='center')
    plot.annotate(r'350$\,\mu$m',(0.900,0.05),ha='center',color='k',weight='bold',size='x-large',family='sans-serif',va='center')
    plot.annotate(r'500$\,\mu$m',(0.600,0.05),ha='center',color='k',weight='bold',size='x-large',family='sans-serif',va='center')
    plot.annotate(r'250$\,\mu$m',(1.250,0.78),ha='center',color='k',weight='bold',size='x-large',family='sans-serif',va='center')
    plot.annotate(r'350$\,\mu$m',(0.900,0.78),ha='center',color='k',weight='bold',size='x-large',family='sans-serif',va='center')
    plot.annotate(r'500$\,\mu$m',(0.600,0.78),ha='center',color='k',weight='bold',size='x-large',family='sans-serif',va='center')

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if figRSRF =='-':
        figname='../Docs/Paper/Figs/rsrf'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+figRSRF+'_rsrf'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    ########################################################

def plotbeam(beamx,beamy,beams,brad,fsuff,fignum,figlab,plimin=None):
    
    from numpy import arange,cos,sin,pi
    
    print 'Plotting beams...'
    print 'WARNING: this can be slow'
    #plot beams
    plot.figure(2)
    plot.clf()
    plot.hot()

    titles=array([r'PSM Beam',r'PMW Beam',r'PLW Beam'])

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
            plot.xlabel(r'x ["]')
        if b == 1:
            plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:]))),('','',''),rotation=90)
        else:
            plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:]))),rotation=90)
            plot.ylabel(r'y ["]')
        cb=plot.colorbar(format='%.1f')#,ticks=tlevs)
        cb.set_label(r'Log scale')
        
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

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/beams'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beams'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    #plot.savefig(figname+'.eps',transparent=False)
    
def plotKc(alpharr,Kc,Kce,fsuff,figKc,figlab):
    ##############################################################
    ##plot colour correction factors
    plot.figure(3)
    plot.clf()
    plot.plot(alpharr,Kc[:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(alpharr,Kc[:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(alpharr,Kc[:,2],'k:',label=r'500$\,\mu$m')
    #plot.plot(alpharr,Kce[:,0],'b--',label=r'$K_\mathrm{c,ext}$ (PSW)')
    #plot.plot(alpharr,Kce[:,1],'g--',label=r'$K_\mathrm{c,ext}$ (PMW)')
    #plot.plot(alpharr,Kce[:,2],'r--',label=r'$K_\mathrm{c,ext}$ (PLW)')
        
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{ColP}(\alpha)$')
    #plot.title('Colour correction factors')
    plot.legend(ncol=3,loc='lower left',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if figKc =='-':
        figname='../Docs/Paper/Figs/colcorr'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+figKc+'_colcorr'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    ###############################################################

def plotbeam_cut(cutx,cuty,beam_cutx,beam_cuty,brad,fsuff,fignum,figlab):
    ##plot beam cuts
    plot.figure(4)
    plot.clf()
    plot.subplot(3,1,1)
    plot.plot(cutx,beam_cutx[:,0,0],'g-',label=r'x (mid)')
    plot.plot(cutx,beam_cutx[:,1,0],'r-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,2,0],'b-',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,0,0],'g:',label=r'y (mid)')
    plot.plot(cuty,beam_cuty[:,1,0],'r:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,2,0],'b:',label='_nolegend_')
    plot.yscale('log')
    plot.ylabelr('PSW Beam [Log]')
    plot.xlim(0,brad)
    plot.ylim(1.e-8,1.)
    plot.legend(loc='upper right',ncol=2)
    
    plot.subplot(3,1,2)
    plot.plot(cutx,beam_cutx[:,0,1],'g-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,1,1],'r-',label=r'x (low)')
    plot.plot(cutx,beam_cutx[:,2,1],'b-',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,0,1],'g:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,1,1],'r:',label=r'y (low)')
    plot.plot(cuty,beam_cuty[:,2,1],'b:',label='_nolegend_')
    plot.yscale('log')
    plot.ylabel(r'PMW Beam (Log)')
    plot.xlim(0,brad)
    plot.ylim(1.e-8,1.)
    plot.legend(loc='upper right',ncol=2)
    
    plot.subplot(3,1,3)
    plot.plot(cutx,beam_cutx[:,0,2],'g-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,1,2],'r-',label='_nolegend_')
    plot.plot(cutx,beam_cutx[:,2,2],'b-',label=r'x (high)')
    plot.plot(cuty,beam_cuty[:,0,2],'g:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,1,2],'r:',label='_nolegend_')
    plot.plot(cuty,beam_cuty[:,2,2],'b:',label=r'y (high)')
    plot.yscale('log')
    plot.ylabel(r'PLW Beam [Log]')
    plot.xlabel(r'Arcsec')
    plot.xlim(0,brad)
    plot.ylim(1.e-8,1.)
    plot.legend(loc='upper right',ncol=2)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamcuts'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamcuts'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
def plotareaf(nuarr,area_bm,ilim,avgnu,avgarea,ith,fsuff,fignum,figlab,unit='sr'):
    
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
        ytit=r'Area [sr]'
    else: ytit=r'Area [arcsec$^2$]'
    
    plot.figure(5)
    plot.clf()
    
    plot.plot(nuarr[ilim[0,0]:ilim[0,1]]/1.e9,area_bm[ilim[0,0]:ilim[0,1],ith,0],'b-',label=r'PSW')
    plot.plot(nuarr[ilim[1,0]:ilim[1,1]]/1.e9,area_bm[ilim[1,0]:ilim[1,1],ith,1],'g-',label=r'PMW')
    plot.plot(nuarr[ilim[2,0]:ilim[2,1]]/1.e9,area_bm[ilim[2,0]:ilim[2,1],ith,2],'r-',label=r'PLW')
    plot.plot(nu0,area0,markeredgecolor='k',marker='^',linestyle='None',label=r'Current beam areas')
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
    plot.xlabel(r'Frequency [GHz]')
    plot.legend(loc='upper right',numpoints=1)
    
    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamarea_nu'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamarea_nu'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotarear(radarr,areas,fsuff,fignum,figlab,unit='sr'):
    
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
        ytit=r'Area [sr]'
    else: ytit=r'Area [sq-arcsec$^2$]'
      
    plot.figure(6)
    plot.clf()
    plot.subplot(2,1,1)    
    plot.plot(radarr,areas[:,0],'b-',label=r'PSW')
    plot.plot(radarr,areas[:,1],'g-',label=r'PMW')
    plot.plot(radarr,areas[:,2],'r-',label=r'PLW')
    plot.axhline(y=area0[0],color='k',linestyle='--',label='Current beam area')    
    plot.axhline(y=area0[0],color='b',linestyle='--')
    plot.axhline(y=area0[1],color='g',linestyle='--')
    plot.axhline(y=area0[2],color='r',linestyle='--')
    #plot.title('Beam areas')
    plot.ylabel(ytit)
    plot.xlabel(r'Radius ["]')
    #plot.legend(loc='lower right',numpoints=1,ncol=2)
    
    na=areas[:,0].size
    plot.subplot(2,1,2)
    plot.plot(radarr,(areas[:,0]-areas[na-1,0])/areas[na-1,0],'b-',label=r'PSW')
    plot.plot(radarr,(areas[:,1]-areas[na-1,1])/areas[na-1,1],'g-',label=r'PMW')
    plot.plot(radarr,(areas[:,2]-areas[na-1,2])/areas[na-1,2],'r-',label=r'PLW')
    
    plot.axhline((area0[0]-areas[na-1,0])/areas[na-1,0],color='k',linestyle='--',label=r'Current beam area')
    plot.axhline((area0[0]-areas[na-1,0])/areas[na-1,0],color='r',linestyle='--')
    plot.axhline((area0[1]-areas[na-1,1])/areas[na-1,1],color='g',linestyle='--')
    plot.axhline((area0[2]-areas[na-1,2])/areas[na-1,2],color='b',linestyle='--')
    plot.axhline(0,color='k',linestyle='-')
    plot.xlabel(r'Radius ["]')
    plot.ylabel(r'$\frac{Area(r)}{Area(r_\mathrm{max})} - 1$')
    plot.legend(loc='upper right',ncol=2)
    plot.ylim(-0.05,0.05)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/beamarea_rad'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamarea_rad'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False) 
    

def plotbeam_sm(radarr,beam_sm,beam_min,beam_max,beam_sd,brad,fsuff,fignum,fignum_all,figlab,figlab_all):

    import matplotlib.pyplot as plot

    ##############################################################
    ##Plot az-averaged beams
    ##############################################################
    plot.figure(7)
    plot.clf()
    plot.subplot(3,1,1)
    plot.plot(radarr,beam_sm[:,0],'b-',label=r'Mean')
    plot.plot(radarr,beam_min[:,0],'b--',label=r'Min/Max')
    plot.plot(radarr,beam_max[:,0],'b--',label=r'_nolegend_')
    plot.plot(radarr,(beam_sm+beam_sd)[:,0],'b:',label=r'Std. Dev.')
    plot.plot(radarr,(beam_sm-beam_sd)[:,0],'b:',label=r'_nolegend_')
    #plot.title('Azimuthally averaged beams')
    plot.ylabelr('PSW')
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    #plot.xscale('log')
    plot.legend(loc='upper right',ncol=3)
    
    plot.subplot(3,1,2)
    plot.plot(radarr,beam_sm[:,1],'g-',label=r'Mean')
    plot.plot(radarr,beam_min[:,1],'g--',label=r'Min/Max')
    plot.plot(radarr,beam_max[:,1],'g--',label=r'_nolegend_')
    plot.plot(radarr,(beam_sm+beam_sd)[:,1],'g:',label=r'Std. Dev.')
    plot.plot(radarr,(beam_sm-beam_sd)[:,1],'g:',label=r'_nolegend_')
    plot.ylabel(r'PMW')
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    #plot.xscale('log')
    
    plot.subplot(3,1,3)
    plot.plot(radarr,beam_sm[:,2],'r-',label=r'Mean')
    plot.plot(radarr,beam_min[:,2],'r--',label=r'Min/Max')
    plot.plot(radarr,beam_max[:,2],'r--',label=r'_nolegend_')
    plot.plot(radarr,(beam_sm+beam_sd)[:,2],'r:',label=r'Std. Dev.')
    plot.plot(radarr,(beam_sm-beam_sd)[:,2],'r:',label=r'_nolegend_')
    plot.ylabel(r'PLW')
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    #plot.xscale('log')
    plot.xlabel(r'Radius [arcsec]')

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/azbeam_range'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_azbeam_range'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False) 

    
    ###################################################
    ##plot all three az-averaged beams
    plot.figure(8)
    plot.clf()
    plot.plot(radarr,beam_sm[:,0],'b-',label=r'PSW')
    plot.plot(radarr,beam_sm[:,1],'g-',label=r'PMW')
    plot.plot(radarr,beam_sm[:,2],'r-',label=r'PLW')
    #plot.plot(radarr,beam_max[:,0],'b--',label='Max val')
    #plot.plot(radarr,beam_max[:,1],'g--',label='_nolegend_')
    #plot.plot(radarr,beam_max[:,2],'r--',label='_nolegend_')
    
    plot.yscale('log')
    plot.xlim(0.,brad)
    plot.ylim(1.e-9,1.)
    plot.ylabel(r'Beam profile')
    plot.xlabel(r'Radius [arcsec]')
    plot.legend(loc='upper right')

    if figlab_all != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab_all,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum_all =='-':
        figname='../Docs/Paper/Figs/azbeam_all'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum_all+'_azbeam_all'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False) 

    
def plotK4(alpharr,K4,K4e,fsuff,figK4,figlab):
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

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if figK4 =='-':
        figname='../Docs/Paper/Figs/K4corr'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+figK4+'_K4corr'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    ###############################################################

def plotK6(arr_a,arr_t,K6,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    
    tit=(r'PSW',r'PMW',r'PLW')
    
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

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/K6corr'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_K6corr'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
        
    #############################################################

def plotKc3(arr_a,arr_t,Kc3,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    
    tit=(r'PSW',r'PMW',r'PLW')
    
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

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc3corr'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kc3corr'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
  
    #############################################################

def plotKc3t(arr_a,arr_t,Kc3t,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    
    tit=(r'PSW',r'PMW',r'PLW')
    
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

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc3tcorr'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kc3tcorr'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)    
    
    #############################################################


def plotareaa(alpharr,areaeff,areapip,areaa_fit,fitpar,fsuff,fignum,figlab):
      
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    #area0=array([426.,771.,1626.])
    #alpha_nep=array((1.26,1.39,1.45))
    
    xtickv=arange(-4,5,1)
    nt=xtickv.size
    xtickl=[]
    xtickbl=[]
    for n in arange(nt): 
        xtickl.append(str(xtickv[n]))
        xtickbl.append('')
    
    bandStr=[r'$250\,\mu$m',r'$350\,\mu$m',r'$500\,\mu$m']
    #cols=['b','g','r']
    styles=['-','--',':']
    plot.figure(13)
    plot.clf()
    for b in range(3):
        plot.plot(alpharr,areaeff[:,b]/areapip[b],c='k',ls=styles[b],label=bandStr[b])
        #plot.plot(alpha_nep[b],areas[b],c=cols[b],ls='x')
    #plot.axhline(y=areas[0],color='b',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    #plot.axhline(y=area0[0],color='b',linestyle=':',label=r'$\Omega_0$')
    #plot.axvline(x=alpha_nep[0],color='b',linestyle=':')
    plot.ylabel('Relative Beam Sold Angle')
    plot.xlabel(r'Spectral index, $\alpha$')
    
    plot.legend(loc='lower left',ncol=1,frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamarea_alpha'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamarea_alpha'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    plot.gca().legend_=None

    #plot.xticks(xtickv,xtickl)
    #plot.title('PSW')
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    afita=len(fitpar[:,0])
    print afita
    for b in range(3):
        plot.plot(alpharr,areaa_fit[:,b]/areapip[b],c='r',ls=styles[b],lw=1)
        lab=r'%s: $\Omega/\Omega_\mathrm{pip}='%bandStr[b]
        apow=afita-1
        print 0,apow,fitpar[0,b]
        for d in range(0,afita):
            apow=afita-d-1
            print d,apow,fitpar[d,b]
            if fitpar[d,b]<0:
                lab=lab+'-'
            elif d>0:
                lab=lab+'+'
            lab=lab+r'%.3g'%(abs(fitpar[d,b]))
            if apow>0:
                lab=lab+r'\alpha'
            if apow>1:
                lab=lab+r'^%d'%(apow)
        lab=lab+r'$'
        lab = lab+' ; $\Omega_\mathrm{pip}=%.1f$ arcsec$^2$'%(areapip[b])
        #plot.text(x0+xr*0.05,y0+yr*0.05*(1+b),r'%s: $\Omega/\Omega_0=%.3g\alpha^2 + %.3g\alpha + %.3g$'%(bandStr[b],fitpar[0,b],fitpar[1,b],fitpar[2,b]),color=cols[b])
        plot.text(x0+xr*0.05,y0+yr*0.05*(1+b),lab,color='k')
    plot.legend(loc='upper right',ncol=1,frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamarea_alpha_withfitdeg'+str(afita-1)+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamarea_alpha_withfitdeg'+str(afita-1)+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotareat(temparr,beta1,beta2,area_efft1,area_efft2,areapip,areat1_fit,areat2_fit,fitpart1,fitpart2,fsuff,fignum,figlab):
      
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    #area0=array([426.,771.,1626.])
    #alpha_nep=array((1.26,1.39,1.45))
    
    xtickv=arange(-4,5,1)
    nt=xtickv.size
    xtickl=[]
    xtickbl=[]
    for n in arange(nt): 
        xtickl.append(str(xtickv[n]))
        xtickbl.append('')
    
    bandStr=[r'$250\,\mu$m',r'$350\,\mu$m',r'$500\,\mu$m']
    #cols=['b','g','r']
    plot.figure(54)
    plot.clf()
    
    plot.plot(temparr,area_efft1[:,0]/areapip[0],'k-',label=r'%s ($\beta=%.1f$)'%(bandStr[0],beta1))
    plot.plot(temparr,area_efft1[:,1]/areapip[1],'k--',label=r'%s ($\beta=%.1f$)'%(bandStr[1],beta1))
    plot.plot(temparr,area_efft1[:,2]/areapip[2],'k:',label=r'%s ($\beta=%.1f$)'%(bandStr[2],beta1))
    plot.plot(temparr,area_efft2[:,0]/areapip[0],'ko',markevery=10,label=r'%s ($\beta=%.1f$)'%(bandStr[0],beta2))
    plot.plot(temparr,area_efft2[:,1]/areapip[1],'ks',mfc='None',markevery=10,label=r'%s ($\beta=%.1f$)'%(bandStr[1],beta2))
    plot.plot(temparr,area_efft2[:,2]/areapip[2],'ko',mfc='None',markevery=10,label=r'%s ($\beta=%.1f$)'%(bandStr[2],beta2))
    
    plot.ylabel('Relative Beam Solid Angle')
    plot.xlabel(r'Temperature [K]')
    
    plot.legend(loc='upper right',ncol=1,frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamarea_temp'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamarea_temp'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    plot.gca().legend_=None

    #plot.xticks(xtickv,xtickl)
    #plot.title('PSW')
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    afitt=len(fitpart1[:,0])
    print afitt
    
    plot.plot(temparr,areat1_fit[:,0]/areapip[0],'r-',lw=1)
    plot.plot(temparr,areat1_fit[:,1]/areapip[1],'r--',lw=1)
    plot.plot(temparr,areat1_fit[:,2]/areapip[2],'r:',lw=1)
    plot.plot(temparr,areat2_fit[:,0]/areapip[0],'ro',markevery=10,lw=1)
    plot.plot(temparr,areat2_fit[:,1]/areapip[1],'rs',mfc=None,markevery=10,lw=1)
    plot.plot(temparr,areat2_fit[:,2]/areapip[2],'ro',mfc=None,markevery=10,lw=1)
            
    for b in range(3):
        lab=r'%s: $\Omega/\Omega_\mathrm{pip}='%bandStr[b]
        for d in range(0,afitt):
            tpow=afitt-d-1
            print d,tpow,fitpart1[d,b]
            if fitpart1[d,b]<0:
                lab=lab+'-'
            elif d>0:
                lab=lab+'+'
            lab=lab+r'%.3g'%(abs(fitpart1[d,b]))
            if tpow>0:
                lab=lab+r'T'
            if tpow>1:
                lab=lab+r'^%d'%(tpow)
        lab=lab+r'$'
        #lab = lab+' ; $\Omega_\mathrm{pip}=%.1f$ arcsec$^2$'%(areapip[b])
        #plot.text(x0+xr*0.05,y0+yr*0.05*(1+b),r'%s: $\Omega/\Omega_0=%.3g\alpha^2 + %.3g\alpha + %.3g$'%(bandStr[b],fitpar[0,b],fitpar[1,b],fitpar[2,b]),color=cols[b])
        plot.text(x0+xr*0.05,y0+yr*0.05*(1+b),lab,color='k')
    plot.legend(loc='upper right',ncol=1,frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamarea_temp_withfitdeg'+str(afitt-1)+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamarea_temp_withfitdeg'+str(afitt-1)+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotbeamcomp(radarr,beam_sm,beam_mod,area_rarr,area_rel,fsuff,fignum,figlab):
    
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
    
    #plot.subplot(2,1,1)
    legnep,=plot.plot(radarr,beam_sm[:,0],'k-',label=r'$B_\mathrm{Nep}$')
    sm0=plot.plot(radarr,beam_sm[:,0],'k-')
    sm1=plot.plot(radarr,beam_sm[:,1],'k-')
    sm2=plot.plot(radarr,beam_sm[:,2],'k-')
    legmod,=plot.plot(radarr,beam_mod[:,0],'k--',label=r'$B_\mathrm{mod}$')
    plot.plot(radarr,beam_mod[:,0],'k--')
    plot.plot(radarr,beam_mod[:,1],'k--')
    plot.plot(radarr,beam_mod[:,2],'k--')
    plot.yscale('log')
    #plot.xlabel('Radius ["]')
    plot.ylim(1e-4,1.)
    plot.ylabel(r'Normalised response')
    plot.xlim(0.,100.)
    plot.xlabel(r'Angle [arcsec]')
    
    plot.annotate(r'250$\,\mu$m',(37,5.e-3),ha='right',weight='bold',size='x-large',family='sans-serif',va='center')
    plot.annotate(r'350$\,\mu$m',(55,5.e-3),ha='left',weight='bold',size='x-large',family='sans-serif',va='center')
    plot.annotate(r'500$\,\mu$m',(80,5.e-3),ha='left',weight='bold',size='x-large',family='sans-serif',va='center')
    #plot.xticks(xtickv,xtickbl)
    #plot.axvline(x=fwhm[0]/2,color='b',linestyle=':')
    #plot.axvline(x=fwhm[1]/2,color='g',linestyle=':')
    #plot.axvline(x=fwhm[2]/2,color='r',linestyle=':')
    l1=plot.legend((legnep,legmod),(r'Measured',r'Predicted'),loc='upper right',ncol=1,frameon=False)
    #plot.legend((sm0,sm1,sm2),('PSW','PMW','PLW'),loc='lower left',ncol=3)
    plot.gca().lines.remove(legnep)
    plot.gca().lines.remove(legmod)
    plot.gca().add_artist(l1)
    
    
    #plot.subplot(2,1,2)
    #plot.plot(radarr,beam_mod[:,0]-beam_sm[:,0],'b-',label='PSW')
    #plot.plot(radarr,beam_mod[:,1]-beam_sm[:,1],'g-',label='PWM')
    #plot.plot(radarr,beam_mod[:,2]-beam_sm[:,2],'r-',label='PLW')
    #plot.ylim(-0.02,0.02)
    #plot.xlim(0,100)
    #plot.axvline(x=fwhm[0]/2,color='b',linestyle=':')
    #plot.axvline(x=fwhm[1]/2,color='g',linestyle=':')
    #plot.axvline(x=fwhm[2]/2,color='r',linestyle=':')
    #plot.axhline(y=0.,color='k',linestyle='-')
    #plot.ylabel(r'$B_\mathrm{mod} - B_\mathrm{Nep}$')
    #plot.xticks(xtickv,xtickbl)
    #plot.legend(loc='lower right',ncol=3)

    #print 'min/max (PSW):',min((beam_mod[:,0]-beam_sm[:,0])/beam_sm[:,0]),max((beam_mod[:,0]-beam_sm[:,0])/beam_sm[:,0])
    #print 'min/max (PMW):',min((beam_mod[:,1]-beam_sm[:,1])/beam_sm[:,1]),max((beam_mod[:,1]-beam_sm[:,1])/beam_sm[:,1])
    #print 'min/max (PLW):',min((beam_mod[:,2]-beam_sm[:,2])/beam_sm[:,2]),max((beam_mod[:,2]-beam_sm[:,2])/beam_sm[:,2])

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamcomp'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamcomp'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotareath(thetarr,areas,fsuff,fignum,figlab,unit='sr'):
    
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
    plot.plot(thetarr,areas[0,:],'b-',label=r'PSW')
    plot.plot(thetarr,areas[1,:],'g-',label=r'PMW')
    plot.plot(thetarr,areas[2,:],'r-',label=r'PLW')
    plot.axhline(y=area0[0],color='k',linestyle='--',label=r'Current beam area')    
    plot.axhline(y=area0[0],color='b',linestyle='--')
    plot.axhline(y=area0[1],color='g',linestyle='--')
    plot.axhline(y=area0[2],color='r',linestyle='--')
    #plot.title('Beam areas')
    plot.ylabel(ytit)
    plot.xlabel(r'Source FHWM ["]')
    plot.legend(loc='lower right',numpoints=1,ncol=2)
    
    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
        
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamarea_theta'+fsuff+'_'+unit
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamarea_theta'+fsuff+'_'+unit
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
def plotK6_va(alpharr,thetarr,K6,K5,aplot,fsuff,fignum,figlab):

    plot.figure(16)
    plot.clf()
    cols=['b','g','r']
    tit=[r'PSW',r'PMW',r'PLW']
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

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size='x-large',family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/K6_a%.1f%s'%(alpharr[aplot],fsuff)
    else:
        figname='../Docs/Paper/Figs/Fig%s_K6_a%.1f%s'%(fignum,alpharr[aplot],fsuff)
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
def plotKc3_va(alpharr,thetarr,Kc3,Kc2,aplot,fsuff,fignum,figlab):

    plot.figure(17)
    plot.clf()
    #cols=['b','g','r']
    tit=[r'PSW',r'PMW',r'PLW']
    plot.plot(thetarr,Kc3[aplot,:,0]/1.e6,'k-',label=r'250$\,\mu$m')
    plot.plot(thetarr,Kc3[aplot,:,1]/1.e6,'k--',label=r'350$\,\mu$m')
    plot.plot(thetarr,Kc3[aplot,:,2]/1.e6,'k:',label=r'500$\,\mu$m')
    
    plot.xlabel(r'Gaussian Source FWHM [arcsec]')
    plot.ylabel(r'Conversion parameter')
    plot.legend(loc='lower left',ncol=3,frameon=False)
    plot.annotate(r'$\alpha=%.1f$'%(alpharr[aplot]),(9.e2,7000),ha='right',weight='bold',size='x-large',family='sans-serif')
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(1.e1,1.e4)
    plot.xlim(1.,1.e3)
    #plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc3_a%.1f%s'%(alpharr[aplot],fsuff)
    else:
        figname='../Docs/Paper/Figs/Fig%s_Kc3_a%.1f%s'%(fignum,alpharr[aplot],fsuff)
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotKc3_varel(alpharr,thetarr,Kc3rel,Kc2rel,aplot,fsuff,fignum,figlab):

    plot.figure(17)
    plot.clf()
    #cols=['b','g','r']
    tit=[r'PSW',r'PMW',r'PLW']
    plot.plot(thetarr,Kc3rel[aplot,:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(thetarr,Kc3rel[aplot,:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(thetarr,Kc3rel[aplot,:,2],'k:',label=r'500$\,\mu$m')
    
    plot.xlabel(r'Gaussian Source FWHM [arcsec]')
    plot.ylabel(r'Conversion parameter (to surface brightness)')
    plot.legend(loc='lower left',ncol=3,frameon=False)
    plot.annotate(r'$\alpha=%.1f$'%(alpharr[aplot]),(9.e2,600),ha='right',weight='bold',size='x-large',family='sans-serif')
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(1.e-1,1.e3)
    plot.xlim(1.,1.e3)
    #plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc3rel_a%.1f%s'%(alpharr[aplot],fsuff)
    else:
        figname='../Docs/Paper/Figs/Fig%s_Kc3rel_a%.1f%s'%(fignum,alpharr[aplot],fsuff)
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotKc3t_va(alpharr,thetarr,Kc3t,Kc,aplot,fsuff,fignum,fignumz,figlab,figlabz):

    plot.figure(18)
    plot.clf()
    #cols=['b','g','r']  
    plot.plot(thetarr,Kc3t[aplot,:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(thetarr,Kc3t[aplot,:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(thetarr,Kc3t[aplot,:,2],'k:',label=r'500$\,\mu$m')
    
    plot.xlabel(r'Gaussian Source FWHM [arcsec]')
    plot.ylabel(r'Conversion parameter [Jy per Jy output by pipeline]')
    plot.legend(loc='upper left',ncol=1,frameon=False)
    plot.annotate(r'$\alpha=%.1f$'%(alpharr[aplot]),(90,0.2),ha='right',weight='bold',size='x-large',family='sans-serif')
    plot.yscale('linear')
    plot.xscale('log')
    plot.ylim(0,10)
    plot.xlim(1,100)
    #plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))    

    if figlab != '':
        labX=10**(log10(plot.xlim()[0]) + 0.05*(log10(plot.xlim()[1])-log10(plot.xlim()[0])))
        labY=plot.ylim()[0] + 0.05*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc3trel_a%.1f%s'%(alpharr[aplot],fsuff)
    else:
        figname='../Docs/Paper/Figs/Fig%s_Kc3trel_a%.1f%s'%(fignum,alpharr[aplot],fsuff)
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    plot.figure(33)
    plot.clf()
    #cols=['b','g','r']
    plot.plot(thetarr,Kc3t[aplot,:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(thetarr,Kc3t[aplot,:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(thetarr,Kc3t[aplot,:,2],'k:',label=r'500$\,\mu$m')
    
    plot.xlabel(r'Gaussian Source FWHM [arcsec]')
    plot.ylabel(r'Conversion parameter [Jy per Jy output by pipeline]')
    plot.legend(loc='upper left',ncol=1,frameon=False)
    plot.annotate(r'$\alpha=%.1f$'%(alpharr[aplot]),(9,0.891),ha='right',weight='bold',size='x-large',family='sans-serif')
    plot.yscale('linear')
    plot.xscale('log')
    plot.ylim(0.89,0.93)
    plot.xlim(1e-2,10)
    #plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))    

    if figlabz != '':
        labX=10**(log10(plot.xlim()[0]) + 0.05*(log10(plot.xlim()[1])-log10(plot.xlim()[0])))
        labY=plot.ylim()[0] + 0.05*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlabz,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignumz =='-':
        fignamez='../Docs/Paper/Figs/Kc3t-zoom_a%.1f%s'%(alpharr[aplot],fsuff)
    else:
        fignamez='../Docs/Paper/Figs/Fig%s_Kc3t-zoom_a%.1f%s'%(fignumz,alpharr[aplot],fsuff)
    plot.savefig(fignamez+'.png',transparent=False,dpi=300.)
    plot.savefig(fignamez+'.eps',transparent=False)
    

def plotKc3t_varel(alpharr,thetarr,Kc3trel,Kc,aplot,fsuff,fignum,fignumz,figlab,figlabz):

    from numpy import min,max
    
    #print min(Kc3trel[aplot,:,0]),max(Kc3trel[aplot,:,0])
    #print min(Kc3trel[aplot,:,1]),max(Kc3trel[aplot,:,1])
    #print min(Kc3trel[aplot,:,2]),max(Kc3trel[aplot,:,2])

    plot.figure(18)
    plot.clf()
    #cols=['b','g','r']  
    plot.plot(thetarr,Kc3trel[aplot,:,0]*1.e6,'k-',label=r'250$\,\mu$m')
    plot.plot(thetarr,Kc3trel[aplot,:,1]*1.e6,'k--',label=r'350$\,\mu$m')
    plot.plot(thetarr,Kc3trel[aplot,:,2]*1.e6,'k:',label=r'500$\,\mu$m')
    
    plot.xlabel(r'Gaussian Source FWHM [arcsec]')
    plot.ylabel(r'Conversion parameter [Jy per MJy/sr]')
    plot.legend(loc='upper right',ncol=1,frameon=False)
    plot.annotate(r'$\alpha=%.1f$'%(alpharr[aplot]),(90,0.012),ha='right',weight='bold',size='x-large',family='sans-serif')
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(0.01,1.0)
    plot.xlim(1,100)
    #plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))    

    if figlab != '':
        labX=10**(log10(plot.xlim()[0]) + 0.05*(log10(plot.xlim()[1])-log10(plot.xlim()[0])))
        labY=10**(log10(plot.ylim()[0]) + 0.95*(log10(plot.ylim()[1])-log10(plot.ylim()[0])))
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc3trel_a%.1f%s'%(alpharr[aplot],fsuff)
    else:
        figname='../Docs/Paper/Figs/Fig%s_Kc3trel_a%.1f%s'%(fignum,alpharr[aplot],fsuff)
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    plot.figure(33)
    plot.clf()
    #cols=['b','g','r']
    plot.plot(thetarr,Kc3trel[aplot,:,0]*1.e6,'k-',label=r'250$\,\mu$m')
    plot.plot(thetarr,Kc3trel[aplot,:,1]*1.e6,'k--',label=r'350$\,\mu$m')
    plot.plot(thetarr,Kc3trel[aplot,:,2]*1.e6,'k:',label=r'500$\,\mu$m')
    
    plot.xlabel(r'Gaussian Source FWHM [arcsec]')
    plot.ylabel(r'Conversion parameter [Jy per MJy/sr]')
    plot.legend(loc='lower left',ncol=1,frameon=False)
    plot.annotate(r'$\alpha=%.1f$'%(alpharr[aplot]),(9,0.0015),ha='right',weight='bold',size='x-large',family='sans-serif')
    plot.yscale('log')
    plot.xscale('log')
    plot.ylim(0.001,0.1)
    plot.xlim(1e-2,10)
    #plot.title(r'$\alpha=%.1f$'%(alpharr[aplot]))    

    if figlabz != '':
        labX=10**(log10(plot.xlim()[0]) + 0.05*(log10(plot.xlim()[1])-log10(plot.xlim()[0])))
        labY=10**(log10(plot.ylim()[0]) + 0.95*(log10(plot.ylim()[1])-log10(plot.ylim()[0])))
        plot.annotate(figlabz,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignumz =='-':
        fignamez='../Docs/Paper/Figs/Kc3trel-zoom_a%.1f%s'%(alpharr[aplot],fsuff)
    else:
        fignamez='../Docs/Paper/Figs/Fig%s_Kc3trel-zoom_a%.1f%s'%(fignumz,alpharr[aplot],fsuff)
    plot.savefig(fignamez+'.png',transparent=False,dpi=300.)
    plot.savefig(fignamez+'.eps',transparent=False)
    
def plotKc_all(alpharr,Kc,Kce,Kc3,Kc3t,fsuff,fignum,figlab):
   
   plot.figure(19)
   plot.clf()
   tit=[r'PSW',r'PMW',r'PLW']
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
   plot.ylabel(r'Correction factor')
   plot.gca
   plot.legend((legKc,legKce,legKc3t,legpsw,legpmw,legplw),
               (r'$K_\mathrm{c,pt}$',r'$K_\mathrm{c,ext}$',r'$K_\mathrm{c3}^\mathrm{tot}$',r'PSW',r'PMW',r'PLW'),
                loc='lower left',ncol=2)
   plot.gca().lines.remove(legKc)
   plot.gca().lines.remove(legKce)
   plot.gca().lines.remove(legKc3t)
   plot.gca().lines.remove(legpsw)
   plot.gca().lines.remove(legpmw)
   plot.gca().lines.remove(legplw)

   if figlab != '':
       labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
       labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
       plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
   
   if fignum =='-':
       figname='../Docs/Paper/Figs/Kc_all'+fsuff
   else:
       figname='../Docs/Paper/Figs/Fig'+fignum+'_Kc_all'+fsuff
   plot.savefig(figname+'.png',transparent=False,dpi=300.)
   plot.savefig(figname+'.eps',transparent=False)
     
def plotKcKc3t_thsmall(alpharr,Kc,Kc3t,thplot,fsuff,fignum,figlab):

    plot.figure(20)
    plot.clf()

    tit=[r'PSW',r'PMW',r'PLW']
    cols=['b','g','r']
    
    for b in range(3):
        plot.plot(alpharr,1.e3*(Kc[:,b]-Kc3t[:,b])/Kc[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$10^3 \times \ (K_\mathrm{c}-K_\mathrm{c3}^\mathrm{tot})/K_\mathrm{c}$ (Jy per Jy$_\mathrm{pip}$)')
    plot.xlabel(r'Spectral index $\alpha$')
    plot.legend(loc='upper left',ncol=3)
    plot.ylim(-4,1)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/KcKc3t_thsmall'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_KcKc3t_thsmall'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
        
def plotKc2Kc3_thbig(alpharr,Kc3,Kc2,thplot,fsuff,fignum,figlab):

    plot.figure(21)
    plot.clf()
    tit=[r'PSW',r'PMW',r'PLW']
    cols=['b','g','r']
    for b in range(3):
        plot.plot(alpharr,1.e5*(Kc2[:,b]-Kc3[:,b])/Kc2[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$10^5 \times \ (K_\mathrm{c2}-K_\mathrm{c3})/K_\mathrm{c2}$')
    plot.xlabel(r'Spectral index $\alpha$')

    plot.legend(loc='upper right',ncol=3)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/KcKc3t_thbig'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_KcKc3t_thbig'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
def plotKc2Kc3_thbigrel(alpharr,Kc3rel,Kc2rel,thplot,fsuff,fignum,figlab):

    plot.figure(21)
    plot.clf()
    tit=[r'PSW',r'PMW',r'PLW']
    cols=['b','g','r']
    for b in range(3):
        plot.plot(alpharr,1.e5*(Kc2rel[:,b]-Kc3rel[:,b])/Kc2rel[:,b],linestyle='-',color=cols[b],label=tit[b])
    plot.ylabel(r'$10^5 \times \ (K_\mathrm{c2}-K_\mathrm{c3})/K_\mathrm{c2}$')
    plot.xlabel(r'Spectral index $\alpha$')

    plot.legend(loc='upper right',ncol=3)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/KcKc3t_thbigrel'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_KcKc3t_thbigrel'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
def plotKc2(alpharr,Kce,Kc2,area0,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    plot.figure(22)
    plot.clf()
    #plot.subplot(2,1,1)
    
    #lkc2,=plot.plot(alpharr,Kc2[:,0],'k-',label=r'$K_\mathrm{c2}$')
    #lkce,=plot.plot(alpharr,Kce[:,0]/area0[0],'k--',label=r'$K_\mathrm{c,ext}/\Omega_0$')
    plot.plot(alpharr,Kc2[:,0]/1.e6,'k-',label=r'250$\,\mu$m')
    plot.plot(alpharr,Kc2[:,1]/1.e6,'k--',label=r'350$\,\mu$m')
    plot.plot(alpharr,Kc2[:,2]/1.e6,'k:',label=r'500$\,\mu$m')
        
    #plot.plot(alpharr,Kce[:,0]/area0[0],'b--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PSW)')
    #plot.plot(alpharr,Kce[:,1]/area0[1],'g--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PMW)')
    #plot.plot(alpharr,Kce[:,2]/area0[2],'r--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PLW)')
    
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'Conversion factor $K_\mathrm{ColE}(\alpha,\infty)$ [MJy/sr per Jy]')
    plot.ylim(10,110)
    #plot.title('Colour correction factors')        
    plot.legend(ncol=3,loc='upper right',frameon=False)
    #plot.legend((lkc2,lkce),(r'$K_\mathrm{c2}$',r'$K_\mathrm{c,ext}/\Omega_0$'),ncol=2,loc='upper right')
    #plot.subplot(2,1,2)
    #plot.plot(alpharr,(Kce[:,0]/area0[0])/Kc2[:,0],'b-',label='PSW')
    #plot.plot(alpharr,(Kce[:,1]/area0[1])/Kc2[:,1],'g-',label='PMW')
    #plot.plot(alpharr,(Kce[:,2]/area0[2])/Kc2[:,2],'r-',label='PLW')
    #plot.xlabel(r'Spectral index ($\alpha$)')
    #plot.ylabel(r'$(K_\mathrm{c,ext}/\Omega_0)\, / \, K_\mathrm{c2}$')
    #plot.ylim(1.0,1.12)
    #plot.legend(ncol=3,loc='lower right')

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/colcorr2'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorr2'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    ###############################################################

def plotKc2rel(alpharr,Kce,Kc2rel,area0,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    plot.figure(22)
    plot.clf()
    #plot.subplot(2,1,1)
    
    #lkc2,=plot.plot(alpharr,Kc2[:,0],'k-',label=r'$K_\mathrm{c2}$')
    #lkce,=plot.plot(alpharr,Kce[:,0]/area0[0],'k--',label=r'$K_\mathrm{c,ext}/\Omega_0$')
    plot.plot(alpharr,Kc2rel[:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(alpharr,Kc2rel[:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(alpharr,Kc2rel[:,2],'k:',label=r'500$\,\mu$m')
        
    #plot.plot(alpharr,Kce[:,0]/area0[0],'b--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PSW)')
    #plot.plot(alpharr,Kce[:,1]/area0[1],'g--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PMW)')
    #plot.plot(alpharr,Kce[:,2]/area0[2],'r--',label=r'$K_\mathrm{c,ext}/\Omega_0$ (PLW)')
    
    plot.xlabel(r'Spectral index ($\alpha$)')
    plot.ylabel(r'Conversion factor $K_\mathrm{ColE}(\alpha,\infty)$')
    plot.ylim(0.85,1.05)
    #plot.title('Colour correction factors')        
    plot.legend(ncol=3,loc='upper right',frameon=False)
    #plot.legend((lkc2,lkce),(r'$K_\mathrm{c2}$',r'$K_\mathrm{c,ext}/\Omega_0$'),ncol=2,loc='upper right')
    #plot.subplot(2,1,2)
    #plot.plot(alpharr,(Kce[:,0]/area0[0])/Kc2[:,0],'b-',label='PSW')
    #plot.plot(alpharr,(Kce[:,1]/area0[1])/Kc2[:,1],'g-',label='PMW')
    #plot.plot(alpharr,(Kce[:,2]/area0[2])/Kc2[:,2],'r-',label='PLW')
    #plot.xlabel(r'Spectral index ($\alpha$)')
    #plot.ylabel(r'$(K_\mathrm{c,ext}/\Omega_0)\, / \, K_\mathrm{c2}$')
    #plot.ylim(1.0,1.12)
    #plot.legend(ncol=3,loc='lower right')

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/colcorr2rel'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorr2rel'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    ###############################################################

def plotblim(radarr,beam_scl_lim,beam_cmb_lim,beam_fix,brad,fsuff,fignum,figlab):

    import matplotlib.pyplot as plot

    ##############################################################
    ##Plot az-averaged beams
    ##############################################################
    plot.figure(24)
    plot.clf()
    tit=[r'PSW',r'PMW',r'PLW']
    
    plot.plot(radarr,beam_cmb_lim[:,0,0],'k--',label=r'Low-freq. edge')
    plot.plot(radarr,beam_cmb_lim[:,1,0],'k-',label=r'Band centre')
    plot.plot(radarr,beam_cmb_lim[:,2,0],'k:',label=r'High-freq. edge')

    plot.yscale('log')
    plot.ylim(1.e-6,1)
    #plot.xlim(0,350.)
    plot.ylabel(r'Normalised response')
    plot.xlabel(r'Angle [arcsec]')
    plot.legend(ncol=1,loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamlim'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamlim'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)    
    
def plotbcomb(radarr,beam_scl,beam_fix,beam_cmb,bsw,fsuff,fignum,figlab):

    plot.figure(30)
    plot.clf()    
    #plot.plot(radarr,beam_scl[:,0],'k-',lw=1)
    #plot.plot(radarr,beam_scl[:,1],'k--',lw=1)
    #plot.plot(radarr,beam_scl[:,2],'k:',lw=1)
    #plot.plot(radarr,beam_fix[:,0],'b-',lw=1)
    #plot.plot(radarr,beam_fix[:,1],'k--',lw=1)
    #plot.plot(radarr,beam_fix[:,2],'k:',lw=1)
    lpsw,=plot.plot(radarr,beam_cmb[:,0],'k-',label=r'250$\,\mu$m')
    #lcmb,=plot.plot(radarr,beam_cmb[:,0],'k-',label='Combined')
    lpmw,=plot.plot(radarr,beam_cmb[:,1],'k--',label=r'350$\,\mu$m')
    #lfix,=plot.plot(radarr,beam_cmb[:,0],'k--',label='Fixed')
    lplw,=plot.plot(radarr,beam_cmb[:,2],'k:',label=r'500$\,\mu$m')
    #lscl,=plot.plot(radarr,beam_cmb[:,0],'k:',label='Scaled')
    #plot.plot(bsw[0,:],beam_cmb[bsw[0,:],0],'bx')
    #plot.plot(bsw[1,:],beam_cmb[bsw[1,:],1],'gx')
    #plot.plot(bsw[2,:],beam_cmb[bsw[2,:],2],'rx')
    plot.yscale('log')
    plot.ylabel(r'Normalised response')
    plot.xlabel(r'Angle [arcsec]')
    plot.legend(ncol=1,loc='upper right',frameon=False)
    
    #plot.gca().lines.remove(lcmb)
    #plot.gca().lines.remove(lscl)
    #plot.gca().lines.remove(lfix)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum =='-':
        figname='../Docs/Paper/Figs/beamcomb'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beamcomb'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False) 

def plot1(nu1a,et1a,ang1b_1,beam1b_1,ang1b_2,beam1b_2,ang1b_3,beam1b_3,figAntEdge,figAntAng,labAntEdge,labAntAng):
    
    plot.figure(34)
    plot.clf()
    #plot.subplot(1,2,1)
    plot.plot(nu1a,et1a,'k-')
    plot.xlabel(r'Frequency [normalised to band centre]')
    plot.xlim(0.8,1.2)
    plot.xticks(arange(0.8,1.3,0.1))
    plot.ylabel(r'Edge taper [dB]')
    plot.ylim(5,12)
    plot.yticks(arange(5,13,1))

    if labAntEdge != '':
        labX=plot.xlim()[0] + 0.95*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.05*(plot.ylim()[1]-plot.ylim()[0])    
        plot.annotate(labAntEdge,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='right')

    if figAntEdge == '-':
        figname='../Docs/Paper/Figs/freq_edgetaper'
    else:
        figname='../Docs/Paper/Figs/Fig'+figAntEdge+'_freq_edgetaper'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    plot.figure(35)
    plot.clf()
    #plot.subplot(1,2,2)
    plot.plot(ang1b_2,beam1b_2,'k-',label=r'Band centre')
    plot.plot(ang1b_1,beam1b_1,'k--',label=r'Low-freq edge')
    plot.plot(ang1b_3,beam1b_3,'k:',label=r'High-freq edge')
    plot.xlabel(r'Angle [units of $\lambda_0/D$]')
    plot.xlim(0,3)
    plot.ylabel(r'Relative response')
    plot.ylim(0,1)

    plot.legend(loc='upper right', frameon=False)

    if labAntAng != '':
        labX=plot.xlim()[0] + 0.95*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.05*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(labAntAng,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='right')

    if figAntAng == '-':
        figname='../Docs/Paper/Figs/freq_angle_beam-response'
    else:
        figname='../Docs/Paper/Figs/Fig'+figAntAng+'_freq_angle_beam-response'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
def plot2(nu2a,fwhm2a_1,fwhm2a_2,nu2b,area2b_1,area2b_2,figAntFWHM,figAntArea,labAntFWHM,labAntArea):
    
    plot.figure(36)
    plot.clf()
    plot.plot(nu2a,fwhm2a_1,'k-',label=r'Beam FWHM')
    plot.plot(nu2a,fwhm2a_2,'k--',label=r'Power-law fit')
    plot.xlabel(r'Normalised frequency')
    plot.xlim(0.8,1.2)
    plot.ylabel(r'Normalised FWHM')
    plot.ylim(0.8,1.2)
    plot.xticks(arange(0.8,1.3,0.1))
    plot.yticks(arange(0.8,1.3,0.1))

    if labAntFWHM != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.05*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(labAntFWHM,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    plot.legend(loc='upper right',frameon=False)
    
    if figAntFWHM == '-':
        figname='../Docs/Paper/Figs/freq_fwhm'
    else:
        figname='../Docs/Paper/Figs/Fig'+figAntFWHM+'_freq_fwhm'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    plot.figure(37)
    plot.clf()
    plot.plot(nu2b,area2b_1,'k-',label='Beam solid angle')
    plot.plot(nu2b,area2b_2,'k--',label='Power-law fit')
    plot.xlabel(r'Normalised frequency')
    plot.xlim(0.8,1.2)
    plot.ylabel(r'Normalised Beam Solid Angle')
    plot.ylim(0.7,1.5)
    plot.xticks(arange(0.8,1.3,0.1))
    plot.yticks(arange(0.7,1.6,0.1))

    if labAntArea != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.05*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(labAntArea,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    plot.legend(loc='upper right',frameon=False)

    if figAntArea == '-':
        figname='../Docs/Paper/Figs/freq_beam-area'
    else:
        figname='../Docs/Paper/Figs/Fig'+figAntArea+'_freq_beam-area'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)    

def plot3(nu3,ap3,fignum,figlab):
    
    plot.figure(38)
    plot.clf()
    plot.plot(nu3,ap3,'k-')
    plot.xlabel(r'Normalised frequency')
    plot.xlim(0,1.5)
    plot.xticks(arange(0,1.75,0.25))
    #plot.xticks(arange(0,1.1,0.1))
    plot.ylabel(r'Aperture efficiency')
    plot.ylim(0,1)
    #plot.yticks(arange(0.94,1.16,0.02))
    plot.axvline(5./6.,c='k',ls='--')
    plot.axvline(7./6.,c='k',ls='--')
    #plot.grid(c='0',ls=':')
    plot.arrow(1.,0.95,1./6.,0.,length_includes_head=True,head_length=0.05,head_width=0.02,color='k')
    plot.arrow(1.,0.95,-1./6.,0.,length_includes_head=True,head_length=0.05,head_width=0.02,color='k')
    plot.annotate(r'Passband',(1.,0.9),xytext=(1.,0.9),ha='center')
    plot.annotate(r'(R=3)',(1.,0.85),xytext=(1.,0.85),ha='center')

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/freq_ap-eff'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_freq_ap-eff'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plot4(nu4,area4,ap4,fignum,figlab):
    
    plot.figure(39)
    plot.clf()
    plot.plot(nu4,area4,'k-',label=r'Beam Solid Angle')
    plot.plot(nu4,ap4,'k--',label=r'Aperture Efficiency')
    plot.xlabel(r'Normalised frequency')
    plot.xlim(0.8,1.2)
    plot.xticks(arange(0.8,1.3,0.1))
    plot.ylabel(r'Normalised quantity')
    plot.ylim(0.7,1.5)
    plot.yticks(arange(0.7,1.6,0.1))
    
    plot.legend(loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/freq_abs_beam-area_ap-eff'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_freq_abs_beam-area_ap-eff'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotKct(temparr1,Kct1,beta1,temparr2,Kct2,beta2,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    ###############################################################
    from scipy import where
    
    plot.figure(39)
    plot.clf()
       
    plot.plot(temparr1,Kct1[:,0],'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct1[:,1],'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct1[:,2],'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr2,Kct2[:,0],'ko',markevery=10,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct2[:,1],'ks',mfc='None',markevery=10,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct2[:,2],'ko',mfc='None',markevery=10,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)

    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{ColP}(T,\beta)$')
    plot.ylim(0.85,1.10)
    plot.xlim(0,40)
    #plot.title('Colour correction factors')
    plot.legend(ncol=1,loc='upper right',frameon=False)
    
    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/colcorrt'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorrt'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

###    plot.figure(40)
###    plot.clf()
###
###    plot.plot(temparr1,Kct1[:,0],'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
###    plot.plot(temparr1,Kct1[:,1],'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
###    plot.plot(temparr1,Kct1[:,2],'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
###    plot.plot(temparr2,Kct2[:,0],ls='-',c='0.7',lw=4,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
###    plot.plot(temparr2,Kct2[:,1],ls='--',c='0.7',lw=4,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
###    plot.plot(temparr2,Kct2[:,2],ls=':',c='0.7',lw=4,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)
###
###    plot.xlabel(r'Temperature (K)')
###    plot.ylabel(r'Colour Correction factor $K_{ColP}(T,\beta)$')
###    plot.ylim(0.86,1.06)
###    plot.xlim(0,40)
###    #plot.title('Colour correction factors')
###    plot.legend(ncol=2,loc='upper right',frameon=False)
###    
###    if fignum == '-':
###        figname='../Docs/Paper/Figs/colcorrt_grey'+fsuff
###    else:
###        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorrt_grey'+fsuff
###    plot.savefig(figname+'.png',transparent=False,dpi=300.)
###    plot.savefig(figname+'.eps',transparent=False)
###    
 
def plotKct2(temparr1,Kct21,beta1,temparr2,Kct22,beta2,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors over temperatures
    ###############################################################

    plot.figure(41)
    plot.clf()
    
    plot.plot(temparr1,Kct21[:,0]/1.e6,'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21[:,1]/1.e6,'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21[:,2]/1.e6,'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr2,Kct22[:,0]/1.e6,'ko',markevery=10,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22[:,1]/1.e6,'ks',mfc='None',markevery=10,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22[:,2]/1.e6,'ko',mfc='None',markevery=10,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)

    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{ColE}(T,\beta)$ [MJy/sr per Jy]')
    plot.ylim(0.,140.)
    plot.xlim(0,40)
    #plot.title('Colour correction factors')
    plot.legend(ncol=2,loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/colcorr2t'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorr2t'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

###    ###Plot grey version
###
###    plot.figure(42)
###    plot.clf()
###    
###    plot.plot(temparr1,Kct21[:,0]/1.e6,'k-',label=r'250$\,\mu$m ($\beta=%.1f$'%beta1)
###    plot.plot(temparr1,Kct21[:,1]/1.e6,'k--',label=r'350$\,\mu$m ($\beta=%.1f$'%beta1)
###    plot.plot(temparr1,Kct21[:,2]/1.e6,'k:',label=r'500$\,\mu$m ($\beta=%.1f$'%beta1)
###    plot.plot(temparr2,Kct22[:,0]/1.e6,ls='-',c='0.7',lw=4,label=r'250$\,\mu$m ($\beta=%.1f)$'%beta2)
###    plot.plot(temparr2,Kct22[:,1]/1.e6,ls='--',c='0.7',lw=4,label=r'350$\,\mu$m ($\beta=%.1f)$'%beta2)
###    plot.plot(temparr2,Kct22[:,2]/1.e6,ls=':',c='0.7',lw=4,label=r'500$\,\mu$m ($\beta=%.1f)$'%beta2)
### 
###    plot.xlabel(r'Temperature (K)')
###    plot.ylabel(r'Colour Correction factor $K_\mathrm{ColE}(T,\beta)$ [MJy/sr per Jy]')
###    plot.ylim(0.,130.)
###    plot.xlim(0,40)
###    #plot.title('Colour correction factors')
###    plot.legend(ncol=2,loc='upper right',frameon=False)
###
###    if figlab != '':
###        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
###        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
###        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
###    
###    if fignum == '-':
###        figname='../Docs/Paper/Figs/colcorr2t_grey'+fsuff
###    else:
###        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorr2t_grey'+fsuff
###    plot.savefig(figname+'.png',transparent=False,dpi=300.)
###    plot.savefig(figname+'.eps',transparent=False)

def plotKct2rel(temparr1,Kct21rel,beta1,temparr2,Kct22rel,beta2,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors over temperatures
    ###############################################################

    plot.figure(41)
    plot.clf()
    
    plot.plot(temparr1,Kct21rel[:,0],'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21rel[:,1],'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21rel[:,2],'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr2,Kct22rel[:,0],'ko',markevery=10,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22rel[:,1],'ks',mfc='None',markevery=10,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22rel[:,2],'ko',mfc='None',markevery=10,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)

    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{ColE}(T,\beta)$')
    plot.ylim(0.85,1.15)
    plot.xlim(0,40)
    #plot.title('Colour correction factors')
    plot.legend(ncol=2,loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/colcorr2trel'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorr2trel'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotKcr(radlist,alphaplot,Kc2rad,bsw,fsuff,fignum,figlab):
    ############################################################
    ## plot colour corrections for range of beam radii
    ############################################################

   
    plot.figure(43)
    plot.clf()

    plot.axvline(bsw[0],c=('0.7'),ls='-')
    plot.axvline(bsw[1],c=('0.7'),ls='--')
    plot.axvline(bsw[2],c=('0.7'),ls=':')
    
    plot.plot(radlist,Kc2rad[:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(radlist,Kc2rad[:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(radlist,Kc2rad[:,2],'k:',label=r'500$\,\mu$m')
#    plot.plot(radlist,Kc2rad[:,0],'kx')
#    plot.plot(radlist,Kc2rad[:,1],'ko')
#    plot.plot(radlist,Kc2rad[:,2],'ks')

    plot.xlabel(r'Beam radius (arcmin)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{ColE}(\alpha,-1)$ [MJy/sr per Jy]')
    plot.legend(ncol=1,loc='upper right',frameon=False)

    plot.xlim(0,1000)
    plot.ylim(0,0.2)
    plot.yticks(arange(0,0.22,0.02))

    #plot.annotate(r'$\alpha=%.1f$'%alphaplot,(0.95*max(radlist),0.95*plot.ylim()[1]),ha='right')
    
    plot.annotate(r'Edges of scaled beams',(bsw[2]+12,0.7*plot.ylim()[1]),color='0.7',ha='center',rotation='vertical')
    plot.annotate(r'250$\,\mu$m',(bsw[0]+12,0.95*plot.ylim()[1]),color='0.7',ha='center',rotation='vertical')
    plot.annotate(r'350$\,\mu$m',(bsw[1]+12,0.95*plot.ylim()[1]),color='0.7',ha='center',rotation='vertical')
    plot.annotate(r'500$\,\mu$m',(bsw[2]+12,0.95*plot.ylim()[1]),color='0.7',ha='center',rotation='vertical')

    plot.title(r'$\alpha=%.1f$'%alphaplot)    

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/colcorr2_R'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorr2_R'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    

def plotuniant(resList,nu0List,alphaList,KcUniAnt,resPlot,nu0Plot,fignumR,fignumNu0,labnumR,labnumNu0):
    
    styles=['-','--',':','-.']

    ###plot colcorr for various Res    
    plot.figure(44)
    plot.clf()
    nRes=len(resPlot)
    for r in range(nRes):
        plot.plot(alphaList,KcUniAnt[resPlot[r],:],c='k',ls=styles[r],label=r'$R=%d$'%(resList[resPlot[r]]))
    
    plot.xlim(-4,4)
    plot.ylim(0.85,1.05)
    plot.xlabel(r'Spectral index $(\alpha)$')
    plot.ylabel(r'Colour correction factor $K_\mathrm{ColP}(\alpha,-1,\nu_0)$')

    plot.legend(loc='lower right',frameon=False)
    
    if labnumR != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(labnumR,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignumR == '-':
        figname='../Docs/Paper/Figs/colcorrUniAntvarR'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignumR+'_colcorrUniAntvarR'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    ###plot colcorr for various nu0
    plot.figure(45)
    plot.clf()
    nNu0=len(nu0Plot)
    for n in range(nNu0):
        if nu0List[nu0Plot[n]]==0:
            label=r'$\nu_0$: centre'
        elif nu0List[nu0Plot[n]] > 0:
            label=r'$\nu_0$: +%d\%%'%nu0List[nu0Plot[n]]
        else:
            label=r'$\nu_0$: %d\%%'%nu0List[nu0Plot[n]]
        plot.plot(alphaList,KcUniAnt[nu0Plot[n],:],c='k',ls=styles[n],label=label)
    
    plot.xlim(-4,4)
    plot.ylim(0.85,1.05)
    plot.xlabel(r'Spectral index $(\alpha)$')
    plot.ylabel(r'Colour correction factor $K_\mathrm{ColP}(\alpha,-1,\nu_0)$')

    plot.legend(loc='lower left',frameon=False)

    if labnumNu0 != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(labnumNu0,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignumNu0 == '-':
        figname='../Docs/Paper/Figs/colcorrUniAntvarnu0'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignumNu0+'_colcorrUniAntvarnu0'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotuniabs(resList,nu0List,alphaList,KcUniAbs,resPlot,nu0Plot,fignumR,fignumNu0,labnumR,labnumNu0):
    
    styles=['-','--',':','-.']

    ###plot colcorr for various Res    
    plot.figure(46)
    plot.clf()
    nRes=len(resPlot)
    for r in range(nRes):
        plot.plot(alphaList,KcUniAbs[resPlot[r],:],c='k',ls=styles[r],label=r'$R=%d$'%(resList[resPlot[r]]))
    
    plot.xlim(-4,4)
    plot.ylim(0.85,1.05)
    plot.xlabel(r'Spectral index $(\alpha)$')
    plot.ylabel(r'Colour correction factor $K_\mathrm{ColP}(\alpha,-1,\nu_0)$')

    plot.legend(loc='lower left',frameon=False)

    if labnumR != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(labnumR,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignumR == '-':
        figname='../Docs/Paper/Figs/colcorrUniAbsvarR'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignumR+'_colcorrUniAbsvarR'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    ###plot colcorr for various nu0
    plot.figure(47)
    plot.clf()
    nNu0=len(nu0Plot)
    for n in range(nNu0):
        if nu0List[nu0Plot[n]]==0:
            label=r'$\nu_0$: centre'
        elif nu0List[nu0Plot[n]] > 0:
            label=r'$\nu_0$: +%d\%%'%nu0List[nu0Plot[n]]
        else:
            label=r'$\nu_0$: %d\%%'%nu0List[nu0Plot[n]]
        plot.plot(alphaList,KcUniAbs[nu0Plot[n],:],c='k',ls=styles[n],label=label)
    
    plot.xlim(-4,4)
    plot.ylim(0.85,1.05)
    plot.xlabel(r'Spectral index $(\alpha)$')
    plot.ylabel(r'Colour correction factor $K_\mathrm{ColP}(\alpha,-1,\nu_0)$')

    plot.legend(loc='lower left',frameon=False)

    if labnumNu0 != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(labnumNu0,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignumNu0 == '-':
        figname='../Docs/Paper/Figs/colcorrUniAbsvarnu0'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignumNu0+'_colcorrUniAbsvarnu0'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)


def plotuniantext(resList,nu0List,alphaList,KcUniAntExt,resPlot,nu0Plot,fignum,figlab):

    ###plot colcorr for various extended source
    plot.figure(48)
    plot.clf()
    nRes=len(resPlot)
    for r in range(nRes):
        plot.plot(alphaList,KcUniAntExt[resPlot[r],:],c='k',ls='-',label=r'$R=%d$'%(resList[resPlot[r]]))
    
    plot.xlim(-4,4)
    plot.ylim(0.85,1.05)
    plot.xlabel(r'Spectral index $(\alpha)$')
    plot.ylabel(r'Colour correction factor $K_\mathrm{ColE}(\alpha,\infty,-1,\nu_0)$')

    #plot.legend(loc='lower left',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum == '-':
        figname='../Docs/Paper/Figs/colcorrUniAntExt'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorrUniAntExt'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotuniabsext(resList,nu0List,alphaList,KcUniAbsExt,resPlot,nu0Plot,fignum,figlab):

    ###plot colcorr for various extended source
    plot.figure(49)
    plot.clf()
    nRes=len(resPlot)
    for r in range(nRes):
        plot.plot(alphaList,KcUniAbsExt[resPlot[r],:],c='k',ls='-',label=r'$R=%d$'%(resList[resPlot[r]]))
    
    plot.xlim(-4,4)
    plot.ylim(0.85,1.05)
    plot.xlabel(r'Spectral index $(\alpha)$')
    plot.ylabel(r'Colour correction factor $K_\mathrm{ColE}(\alpha,\infty,-1,\nu_0)$')

    #plot.legend(loc='lower left',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')
    
    if fignum == '-':
        figname='../Docs/Paper/Figs/colcorrUniAbsExt'
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_colcorrUniAbsExt'
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

def plotKcKc2(alpharr,Kc,Kc2,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    plot.figure(50)
    plot.clf()
    
    plot.plot(alpharr,Kc[:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(alpharr,Kc[:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(alpharr,Kc[:,2],'k:',label=r'500$\,\mu$m')
        
    plot.xlabel(r'Source spectral index ($\alpha_S$)')
    plot.ylabel(r'Colour Correction Factor')
    #plot.ylim(10,110)
    plot.legend(ncol=3,loc='upper right',frameon=False)
    plot.title(r'Colour Correction Factor vs.~Source Spectral index (Point Source)')
    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kc'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    plot.figure(51)
    plot.clf()
    
    plot.plot(alpharr,Kc2[:,0]/1.e6,'k-',label=r'250$\,\mu$m')
    plot.plot(alpharr,Kc2[:,1]/1.e6,'k--',label=r'350$\,\mu$m')
    plot.plot(alpharr,Kc2[:,2]/1.e6,'k:',label=r'500$\,\mu$m')
        
    plot.xlabel(r'Source spectral index ($\alpha_S$)')
    plot.ylabel(r'Colour Correction Factor [MJy/sr per Jy]')
    plot.ylim(10,110)
    plot.legend(ncol=3,loc='upper right',frameon=False)
    plot.title(r'Colour Correction Factor vs.~Source Spectral index (Extended Source)')
    
    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc2'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kc2'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    ###############################################################

def plotKcKc2rel(alpharr,Kc,Kc2rel,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors
    plot.figure(50)
    plot.clf()
    
    plot.plot(alpharr,Kc[:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(alpharr,Kc[:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(alpharr,Kc[:,2],'k:',label=r'500$\,\mu$m')
        
    plot.xlabel(r'Source spectral index ($\alpha_S$)')
    plot.ylabel(r'Colour Correction Factor')
    #plot.ylim(10,110)
    plot.legend(ncol=3,loc='upper right',frameon=False)
    plot.title(r'Colour Correction Factor vs.~Source Spectral index (Point Source)')
    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kc'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    plot.figure(51)
    plot.clf()
    
    plot.plot(alpharr,Kc2rel[:,0],'k-',label=r'250$\,\mu$m')
    plot.plot(alpharr,Kc2rel[:,1],'k--',label=r'350$\,\mu$m')
    plot.plot(alpharr,Kc2rel[:,2],'k:',label=r'500$\,\mu$m')
        
    plot.xlabel(r'Source spectral index ($\alpha_S$)')
    plot.ylabel(r'Colour Correction Factor')
    plot.ylim(0.75,1.05)
    plot.legend(ncol=3,loc='upper right',frameon=False)
    plot.title(r'Colour Correction Factor vs.~Source Spectral index (Extended Source)')
    
    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum =='-':
        figname='../Docs/Paper/Figs/Kc2rel'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kc2rel'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    ###############################################################

def plotKctKct2(temparr1,Kct1,Kct21,beta1,temparr2,Kct2,Kct22,beta2,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors over temperatures
    ###############################################################

    plot.figure(52)
    plot.clf()
    
    plot.plot(temparr1,Kct1[:,0],'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct1[:,1],'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct1[:,2],'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr2,Kct2[:,0],'ko',markevery=10,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct2[:,1],'ks',mfc='None',markevery=10,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct2[:,2],'ko',mfc='None',markevery=10,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)

    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{CP}(T,\beta)$')
    plot.ylim(0.85,1.10)
    plot.xlim(0,40)
    plot.title('Colour correction factors vs Source Temperature (Point Source)')
    plot.legend(ncol=2,loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/Kct'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kct'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    plot.figure(53)
    plot.clf()
    
    plot.plot(temparr1,Kct21[:,0]/1.e6,'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21[:,1]/1.e6,'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21[:,2]/1.e6,'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr2,Kct22[:,0]/1.e6,'ko',markevery=10,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22[:,1]/1.e6,'ks',mfc='None',markevery=10,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22[:,2]/1.e6,'ko',mfc='None',markevery=10,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)

    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{CE}(T,\beta)$ [MJy/sr per Jy]')
    plot.ylim(0.,140.)
    plot.xlim(0,40)
    plot.title('Colour correction factors vs Source Temperature (Extended Source)')
    plot.legend(ncol=2,loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/Kct2'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kct2'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)


    ###############################################################

def plotKctKct2rel(temparr1,Kct1,Kct21rel,beta1,temparr2,Kct2,Kct22rel,beta2,fsuff,fignum,figlab):
    ##############################################################
    ##plot colour correction factors over temperatures
    ###############################################################

    plot.figure(52)
    plot.clf()
    
    plot.plot(temparr1,Kct1[:,0],'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct1[:,1],'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct1[:,2],'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr2,Kct2[:,0],'ko',markevery=10,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct2[:,1],'ks',mfc='None',markevery=10,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct2[:,2],'ko',mfc='None',markevery=10,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)

    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{CP}(T,\beta)$')
    plot.ylim(0.85,1.10)
    plot.xlim(0,40)
    plot.title('Colour correction factors vs Source Temperature (Point Source)')
    plot.legend(ncol=2,loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/Kct'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kct'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)
    
    plot.figure(53)
    plot.clf()
    
    plot.plot(temparr1,Kct21rel[:,0],'k-',label=r'250$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21rel[:,1],'k--',label=r'350$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr1,Kct21rel[:,2],'k:',label=r'500$\,\mu$m ($\beta=%.1f$)'%beta1)
    plot.plot(temparr2,Kct22rel[:,0],'ko',markevery=10,label=r'250$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22rel[:,1],'ks',mfc='None',markevery=10,label=r'350$\,\mu$m ($\beta=%.1f$)'%beta2)
    plot.plot(temparr2,Kct22rel[:,2],'ko',mfc='None',markevery=10,label=r'500$\,\mu$m ($\beta=%.1f$)'%beta2)

    plot.xlabel(r'Temperature (K)')
    plot.ylabel(r'Colour Correction factor $K_\mathrm{CE}(T,\beta)$')
    plot.ylim(0.85,1.15)
    plot.xlim(0,40)
    plot.title('Colour correction factors vs Source Temperature (Extended Source)')
    plot.legend(ncol=2,loc='upper right',frameon=False)

    if figlab != '':
        labX=plot.xlim()[0] + 0.05*(plot.xlim()[1]-plot.xlim()[0])
        labY=plot.ylim()[0] + 0.95*(plot.ylim()[1]-plot.ylim()[0])
        plot.annotate(figlab,(labX,labY),xycoords='data',weight='bold',size=24,family='sans-serif',ha='center',va='center')

    if fignum == '-':
        figname='../Docs/Paper/Figs/Kct2rel'+fsuff
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_Kct2rel'+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)


    ###############################################################

def plotprettybeam(fignum,figlab,fitsIn=None,band='psw'):
    import aplpy as ap
    
    if fitsIn == None:
        if band == 'psw':
            fitsIn='../Inputs/spire_beams_measured/psw_beam_1arcsec.fits'
        elif band == 'pmw':
            '../Inputs/spire_beams_measured/pmw_beam_1arcsec.fits'
        elif band == 'plw':
            '../Inputs/spire_beams_measured/plw_beam_1arcsec.fits'
        else:
            print 'Unknown band: '+band
        
    fig=ap.FITSFigure(fitsIn)
    fig.show_grayscale(stretch='log',vmin=1e-5,vmax=1,invert=True)

    fig.recenter(0,0,width=10/60.,height=10/60.)
    
    fig.hide_tick_labels()
    fig.hide_axis_labels()
    if figlab!='':    
        fig.add_label(0.05,0.95,figlab,relative='True',weight='bold',size=48)

    if fignum == '-':
        figname='../Docs/Paper/Figs/beammap_'+band
    else:
        figname='../Docs/Paper/Figs/Fig'+fignum+'_beammap_'+band
    fig.save(figname+'.png')
    fig.save(figname+'.eps')
    
if __name__ == "__main__":
    maincode()