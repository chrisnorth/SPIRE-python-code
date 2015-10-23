##compare results for point sources and extended sources

from csv import reader
import argparse
import string
import sys

import matplotlib.pyplot as plot

from numpy import array,zeros,size,shape,zeros_like,max,min,NaN
from scipy import where

def maincode():
    parser=argparse.ArgumentParser(description='Herschel-SPIRE Extended Emission Calibration')
    parser.add_argument('--brad',action='store',default=600.,type=float,dest='brad',help='Radius to integrate beam out to (arcsec). Default=350')
    parser.add_argument('--ind',action='store',default=0.85,type=float,dest='ind',help='Power law index with which beam FWHM changes with frequency (FWHM \\propto \\nu^{ind}). Default=0.65')
    parser.add_argument('--bzlim',action='store',default=None,type=float,dest='bzlim',help='Zero-limit for beamfiles, below which the values are replaced by BZVAL. Default=None (i.e. don\'t set limit)')
    parser.add_argument('--bzval',action='store',default=0.,type=float,dest='bzval',help='Value with which to replace beam values below the zero-limit (BZLIM). Default=0')
    parser.add_argument('--rsrf', action='store',default='M',dest='rsrftype',help='RSRF type to use [M=Measured | T=Top-hat]')
    parser.add_argument('--apeff', action='store',default='R',dest='apftype',help='Aperture Efficiency type to use [R=Real | U=Uniform]')
    parser.add_argument('--nepmod', action='store',default='ESA4',dest='nepmod',help='Model to use for Neptune spectrum ["ESA2"|"ESA4"]')
    
    args=parser.parse_args()
    brad=args.brad
    ind=args.ind
    rsrftype=args.rsrftype
    apftype=args.apftype
    nepmod=args.nepmod
    bzlim=args.bzlim
    bzval=args.bzval

    band=range(3)

    fwhm0=array([17.6,23.9, 35.2])
            
    fsuff='_theta_final_newBprofBS_Br%d_Ind%.2f'%(int(brad),ind)
    if nepmod != '':
        fsuff=fsuff+'_%s'%(nepmod)
    if bzlim != None:
        fsuff=fsuff+'_Bl%1g_Bv%1g'%(bzlim,bzval)
    if rsrftype == 'T':
        fsuff=fsuff+'_noRSRF'
    if apftype == 'U':
        fsuff=fsuff+'_noApEff'
        
    print 'File_suffix= '+fsuff
    
    ###Read in Kc###
    
    file='../Outputs/Kc'+fsuff+'.csv'
    kc_input=reader(open(file))
    alist=[]
    kc_psw=[]
    kc_pmw=[]
    kc_plw=[]
    Kpip=zeros(3)
    for row in kc_input:
        if string.find(row[0],'#') < 0:
            alist.append(float(row[0]))
            kc_psw.append(float(row[1]))
            kc_pmw.append(float(row[2]))
            kc_plw.append(float(row[3]))

    alpharr=array(alist)
    nalph=size(alpharr)
    print kc_psw
    Kc=zeros((nalph,3))
    Kc[:,0]=kc_psw
    Kc[:,1]=kc_pmw
    Kc[:,2]=kc_plw

    ###Read in Kc2###
        
    file='../Outputs/Kc2'+fsuff+'.csv'
    kc2_input=reader(open(file))
    kc2_psw=[]
    kc2_pmw=[]
    kc2_plw=[]
    for row in kc2_input:
        if string.find(row[0],'#') < 0:
            kc2_psw.append(float(row[1]))
            kc2_pmw.append(float(row[2]))
            kc2_plw.append(float(row[3]))
            
    nalph=size(alpharr)
    Kc2=zeros((nalph,3))
    Kc2[:,0]=kc2_psw
    Kc2[:,1]=kc2_pmw
    Kc2[:,2]=kc2_plw

    ###Read in Kc2b###
        
    file='../Outputs/Kc2b'+fsuff+'.csv'
    kc2b_input=reader(open(file))
    kc2b_psw=[]
    kc2b_pmw=[]
    kc2b_plw=[]
    for row in kc2b_input:
        if string.find(row[0],'#') < 0:
            kc2b_psw.append(float(row[1]))
            kc2b_pmw.append(float(row[2]))
            kc2b_plw.append(float(row[3]))
            
    nalph=size(alpharr)
    Kc2b=zeros((nalph,3))
    Kc2b[:,0]=kc2b_psw
    Kc2b[:,1]=kc2b_pmw
    Kc2b[:,2]=kc2b_plw

    ###Read in Kc3###
    
    file='../Outputs/Kc3_PSW'+fsuff+'.csv'
    kc3psw_input=reader(open(file))
    thetalist=[]    
    kc3psw_in=[]
    for row in kc3psw_input:
        if string.find(row[0],'#') < 0:
            thetalist.append(float(row[0]))
            kc3x=[]
            for a in range(nalph):
                kc3x.append(float(row[a+1]))
            kc3psw_in.append(array(kc3x))

    thetarr=array(thetalist)

    file='../Outputs/Kc3_PMW'+fsuff+'.csv'
    kc3pmw_input=reader(open(file))
    kc3pmw_in=[]
    for row in kc3pmw_input:
        if string.find(row[0],'#') < 0:
            kc3x=[]
            for a in range(nalph):
                kc3x.append(float(row[a+1]))
            kc3pmw_in.append(array(kc3x))

    file='../Outputs/Kc3_PLW'+fsuff+'.csv'
    kc3plw_input=reader(open(file))
    kc3plw_in=[]
    for row in kc3plw_input:
        if string.find(row[0],'#') < 0:
            kc3x=[]
            for a in range(nalph):
                kc3x.append(float(row[a+1]))
            kc3plw_in.append(array(kc3x))

    nth=size(thetarr)
    Kc3=zeros((nalph,nth,3))
    for a in range(nalph):
        for t in range(nth):
            Kc3[a,t,0]=kc3psw_in[t][a]
            Kc3[a,t,1]=kc3pmw_in[t][a]
            Kc3[a,t,2]=kc3plw_in[t][a]

    ###Read in Kc3t###
        
    file='../Outputs/Kc3t_PSW'+fsuff+'.csv'
    kc3tpsw_input=reader(open(file))
    thetalist=[]    
    kc3tpsw_in=[]
    for row in kc3tpsw_input:
        if string.find(row[0],'#') < 0:
            kc3tx=[]
            for a in range(nalph):
                kc3tx.append(float(row[a+1]))
            kc3tpsw_in.append(array(kc3tx))

    file='../Outputs/Kc3t_PMW'+fsuff+'.csv'
    kc3tpmw_input=reader(open(file))
    kc3tpmw_in=[]
    for row in kc3tpmw_input:
        if string.find(row[0],'#') < 0:
            kc3tx=[]
            for a in range(nalph):
                kc3tx.append(float(row[a+1]))
            kc3tpmw_in.append(array(kc3tx))

    file='../Outputs/Kc3t_PLW'+fsuff+'.csv'
    kc3tplw_input=reader(open(file))
    kc3tplw_in=[]
    for row in kc3tplw_input:
        if string.find(row[0],'#') < 0:
            kc3tx=[]
            for a in range(nalph):
                kc3tx.append(float(row[a+1]))
            kc3tplw_in.append(array(kc3tx))

    Kc3t=zeros((nalph,nth,3))
    for a in range(nalph):
        for t in range(nth):
            Kc3t[a,t,0]=kc3tpsw_in[t][a]
            Kc3t[a,t,1]=kc3tpmw_in[t][a]
            Kc3t[a,t,2]=kc3tplw_in[t][a]
    
    Kc3_Kc2=zeros_like(Kc3)
    Kc3t_Kc=zeros_like(Kc3t)
    Kc3_Kc2_diff=zeros_like(Kc3)
    Kc3t_Kc_diff=zeros_like(Kc3t)
    Kmax=zeros_like(Kc3t)
    for b in band:
        for a in range(nalph):
            Kc3_Kc2[a,:,b]=Kc3[a,:,b]/Kc2[a,b]
            Kc3t_Kc[a,:,b]=Kc3t[a,:,b]/Kc[a,b]
            Kc3_Kc2_diff[a,:,b]=100.*(Kc3[a,:,b]-Kc2[a,b])/Kc2[a,b]
            Kc3t_Kc_diff[a,:,b]=100.*(Kc3t[a,:,b]-Kc[a,b])/Kc[a,b]
            for t in range(nth):
                Kmax[a,t,b]=min([Kc3_Kc2[a,t,b],Kc3t_Kc[a,t,b]])
            Kmax[a,:,b]=where(thetarr < 0.1*fwhm0[b],NaN,Kmax[a,:,b])
            Kmax[a,:,b]=where(thetarr > 10*fwhm0[b],NaN,Kmax[a,:,b])
            
    
    i1=14
    a1=alpharr[i1]
    i2=6
    a2=alpharr[i2]
    print a1,a2
    
    print Kc3[i1,:,0]
    print Kc3[i1,:,1]
    print Kc3[i1,:,2]
    plot.figure(1)
    plot.clf()
    plot.plot(thetarr,Kc3[i1,:,0],'b-')
    plot.plot(thetarr,Kc3[i1,:,1],'g-')
    plot.plot(thetarr,Kc3[i1,:,2],'r-')
    plot.axhline(Kc2[i1,0],c='b',ls=':')
    plot.axhline(Kc2[i1,1],c='g',ls=':')
    plot.axhline(Kc2[i1,2],c='r',ls=':')
    plot.xscale('log')
    plot.yscale('log')
    plot.xlim(0.1,1000)
    plot.ylabel(r'$K_{c3}$ (Partially extended, surface brightness')
    plot.xlabel(r'Source size, $\theta$')
    plot.draw()
    
    plot.figure(2)
    plot.clf()
    plot.plot(thetarr,Kc3t[i1,:,0],'b-')
    plot.plot(thetarr,Kc3t[i1,:,1],'g-')
    plot.plot(thetarr,Kc3t[i1,:,2],'r-')
    plot.axhline(Kc[i1,0],c='b',ls=':')
    plot.axhline(Kc[i1,1],c='g',ls=':')
    plot.axhline(Kc[i1,2],c='r',ls=':')
    plot.xscale('log')
    plot.yscale('log')
    plot.xlim(0.1,1000)
    plot.ylabel(r'$K_{c3}^{tot}$ (Partially extended, flux density')
    plot.xlabel(r'Source size, $\theta$')
    plot.draw()
    
    thetarel=zeros((nth,3))
    for b in band:
        thetarel[:,b]=thetarr[:]/fwhm0[b]
        
    plot.figure(3)
    plot.clf()
    plot.plot(thetarel[:,b],Kc3_Kc2_diff[i1,:,0],'b-')
    plot.plot(thetarel[:,b],Kc3_Kc2_diff[i1,:,1],'g-')
    plot.plot(thetarel[:,b],Kc3_Kc2_diff[i1,:,2],'r-')
    plot.plot(thetarel[:,b],Kc3t_Kc_diff[i1,:,0],'b--')
    plot.plot(thetarel[:,b],Kc3t_Kc_diff[i1,:,1],'g--')
    plot.plot(thetarel[:,b],Kc3t_Kc_diff[i1,:,2],'r--')    
    plot.axhline(10.,c='k',ls=':')         
    plot.axhline(1.,c='k',ls=':')
    plot.xscale('log')
    plot.yscale('log')
    plot.xlim(0.01,100)
    plot.ylim(0.1,1000.)
    plot.ylabel(r'% change by assuming partially resolved source')
    plot.xlabel(r'Source FWHM/Beam FWHM')
    plot.annotate('Total source flux',(0.07,30),rotation=55)
    plot.annotate('Peak surface brightness',(2,45),rotation=-47)
    plot.draw()
    plot.show()
    
    plot.figure(4)
    plot.clf()
    plot.plot(thetarel[:,b],Kc3_Kc2[i1,:,0],'b-')
    plot.plot(thetarel[:,b],Kc3_Kc2[i1,:,1],'g-')
    plot.plot(thetarel[:,b],Kc3_Kc2[i1,:,2],'r-')
    plot.plot(thetarel[:,b],Kc3t_Kc[i1,:,0],'b--')
    plot.plot(thetarel[:,b],Kc3t_Kc[i1,:,1],'g--')
    plot.plot(thetarel[:,b],Kc3t_Kc[i1,:,2],'r--')    
    plot.axhline(1.1,c='k',ls=':')   
    plot.axhline(1.,c='k',ls=':')
    plot.xscale('log')
    plot.yscale('log')
    plot.xlim(0.01,100)
    plot.ylim(0.1,1000.)
    plot.ylabel(r'Partially resolve / Nominal (Point or Fully extended)')
    plot.xlabel(r'Source FWHM/Beam FWHM')
    plot.annotate('Total source flux',(3,500),rotation=55)
    plot.annotate('Peak surface brightness',(0.05,500),rotation=-55)
    plot.draw()
    plot.show()
    
if __name__=="__main__":
    maincode()