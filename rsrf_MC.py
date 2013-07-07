
### Vary RSRF by set limits and run monte-carlo

from numpy import array,arange,zeros,zeros_like,mean,std,min,max,median
from numpy.random import normal
from scipy import where
from scipy.interpolate import interp1d
from scipy.io.idl import readsav
from scipy.integrate import trapz
from csv import reader

from rsrf import getrsrf,getapf,bandedge,readsvg
from calc_conv import calc_k4,calc_kc

import matplotlib as mpl
import matplotlib.pyplot as plot

import sys

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['legend.fontsize']='x-large'
mpl.rcParams['lines.linewidth']=2
mpl.rcParams['axes.linewidth']=2
mpl.rcParams['axes.labelsize']='x-large'
mpl.rcParams['axes.titlesize']='xx-large'
mpl.rcParams['xtick.labelsize']='x-large'
mpl.rcParams['ytick.labelsize']='x-large'
mpl.rcParams['font.family']='serif'
mpl.rcParams['text.usetex']=True
mpl.rcParams['font.weight']='bold'
mpl.rcParams['font.serif']=['Times', 'Palatino', 'New Century Schoolbook', 'Bookman', 'Computer Modern Roman']
mpl.rcParams['font.sans-serif']=['Helvetica', 'Avant Garde', 'Computer Modern Sans serif']
mpl.rcParams['font.cursive']=['Zapf Chancery']
mpl.rcParams['font.monospace']=['Courier', 'Computer Modern Typewriter']

#def maincode():
#
c=299792458. #speed of light
#set band centres
wlc=(250.e-6,350.e-6,500.e-6)
wlc=array(wlc)
nuc=c/wlc
print nuc/1.e12
band=arange(3)
bandStr=['PSW','PMW','PLW']
errStr='Meas'
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
print 'Getting RSRF...'
file='../Inputs/SPIRE-Phot-RSRF.csv'
rsrf_input=reader(open(file))
nnuin=len(list(rsrf_input))

rsrfin=zeros((nnuin,3))
wlin=zeros(nnuin)
nuin=zeros(nnuin)
r=0
rsrf_input=reader(open(file))
for row in rsrf_input:
    wlin[r]=float(row[0])*1.e-6
    rsrfin[r,0]=float(row[1])
    rsrfin[r,1]=float(row[2])
    rsrfin[r,2]=float(row[3])
    r=r+1

nuin=c/(array(wlin))

ilimin=zeros((3,2))
nulimin=zeros((3,2))
ilimin_2=zeros((3,2))
nulimin_2=zeros((3,2))
for b in band:
    (ilimin[b,:],nulimin[b,:])=bandedge(nuin,rsrfin[:,b])
    (ilimin_2[b,:],nulimin_2[b,:])=bandedge(nuin,rsrfin[:,b],fact=2.)
    print 'Band %d input limits: [%.2f:%.2f] GHz ([%.2f:%.2f] um)' % (b,nulimin[b,0]/1.e9,nulimin[b,1]/1.e9,c/nulimin[b,1]*1.e6,c/nulimin[b,0]*1.e6)

rsrfarr0=zeros((nnu,3))
ilim=zeros((3,2))
nulim=zeros((3,2))
ilim_2=zeros((3,2))
nulim_2=zeros((3,2))
for b in band:
    intf=interp1d(nuin,rsrfin[:,b],bounds_error=False,fill_value=0)
    rsrf=intf(nuarr)
    (ilim[b,:],nulim[b,:])=bandedge(nuarr,rsrf)
    (ilim_2[b,:],nulim_2[b,:])=bandedge(nuarr,rsrf,fact=2.)
    print 'Band %d limits: [%.2f:%.2f] GHz ([%.2f:%.2f] um)' % (b,nulim[b,0]/1.e9,nulim[b,1]/1.e9,c/nulim[b,1]*1.e6,c/nulim[b,0]*1.e6)
    #print 'Band %d limits: [%.2f:%.2f] um' % (b,c/nulim[b,1]*1.e6,c/nulim[b,0]*1.e6)
    rsrfarr0[:,b]=rsrf

#########################
###  Set RSRF errors  ###
#########################
relerr=zeros((nnuin,3))
rsrfinsd=zeros((nnuin,3))
if errStr=='Uni':
    relErr=0.1 #constant 10% error
    print 'Using uniform error of %.1f%%...'%(relErr*100.)
    #rsrfinerr=zeros((nnuin,3,2))
    rsrfinsd=zeros((nnuin,3))
    for b in band:
        rsrfinsd[ilimin_2[b,0]:ilimin_2[b,1],b] = relErr * rsrfin[ilimin_2[b,0]:ilimin_2[b,1],b]
elif errStr=='Meas':
    print 'Using measured error...'
    rsrfinsd=zeros((nnuin,3))
    specin=zeros((nnuin,3))
    specinnofilt=zeros((nnuin,3))
    sdin=zeros((nnuin,3))
    for b in band:
        ##set filenames
        filemeas='../RSRF/Locke/PHOTrsrf_%s_coadded-meas.dat'%bandStr[b]
        filesd='../RSRF/Locke/PHOTrsrf_%s_coadded-sd.dat'%bandStr[b]
        
        ##read in measured spectrum and sd
        wn,meanspec,detspec = readsvg(filemeas,verbose=False,dofilt=True)
        wnnofilt,meanspecnofilt,detspecnofilt = readsvg(filemeas,verbose=False,dofilt=False)
        #wnsd,meansd,detsd = readsvg(filesd,wnin=wn,verbose=True,dofilt=True)
        freqsd=wn*c*1.e2
        
        specint=interp1d(freqsd,meanspec,bounds_error=False,fill_value=0.)
        specin[:,b]=specint(nuin)
        specintnofilt=interp1d(freqsd,meanspecnofilt,bounds_error=False,fill_value=0.)
        specinnofilt[:,b]=specintnofilt(nuin)
        
        ##calculate relative error
        nc=len(wn)
        newsd=zeros(nc)
        for coord in range(nc):
            newsd[coord]=std(detspec[coord,:])

        sdint=interp1d(freqsd,newsd,bounds_error=False,fill_value=0.)
        sdin[:,b]=sdint(nuin)
        
        relerrin=abs(newsd/meanspec)
        maxerr=5.
        relerrin=where(relerrin>maxerr,maxerr,relerrin)
        
        ##interpolate to input rsrf grid
        relerrint=interp1d(freqsd,relerrin,bounds_error=False,fill_value=0.)
        relerr[:,b]=relerrint(nuin)
        rsrfinsd[:,b]=relerr[:,b] * rsrfin[:,b]
        print 'Error range: %.3f,%.3f'%(min(relerrin),max(relerrin))
else:
    print 'Unknown error type. Must be ["Meas"|"Uni"]'
    sys.exit()
###############################
##  Get aperture efficiency  ##
###############################

print 'Getting Aperture Efficiency...'

apfarr=zeros((nnu,3))
for b in band:
    apfarr[:,b]=getapf("R",b,nuarr)
    
##############################
##  Calculate KMonP, KColP, Kcal  ##
##############################

alphapip=-1.
alpha_nep=array((1.26,1.39,1.45))
alpharr=arange(-4,5.5,0.5)
nalph=alpharr.size

print 'Calculating Kcal, KMonP,KColP for nominal version...'
KMonP0=zeros((nalph,3))
#KMonPpip0=zeros(3)
KColP0=zeros((nalph,3))
Kcal0=zeros(3)
alpha_nep=array((1.26,1.39,1.45))

nuc0=zeros(3)
nulim0=zeros((3,2))
rsrflim0=zeros(3)
rsrfmed0=zeros(3)
rsrfR0=zeros(3)

rsrfwid0=zeros_like(rsrfarr0)
Kcalwid0=zeros_like(Kcal0)
KMonPwid0=zeros_like(KMonP0)
KColPwid0=zeros_like(KColP0)
for b in band:
    nuc0[b]=trapz(nuarr*rsrfarr0[:,b])/trapz(rsrfarr0[:,b])
    ilim_med,nulim_med=bandedge(nuarr,rsrfarr0[:,b],fact=1.,method='median')
    nulim0[b,:]=[nulim_med[0],nulim_med[1]]
    rsrflim0[b]=nulim0[b,1]-nulim0[b,0]
    rsrfR0[b]=nuc0[b]/rsrflim0[b]
    rsrfwid0[ilim_med[0]:ilim_med[1],b]=1.
    rsrfmed0[b]=median(rsrfarr0[ilim_med[0]:ilim_med[1],b])
    #KMonPpip0[b] = calc_k4(alphapip,rsrfarr0[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
    Kcal0[b] = calc_k4(alpha_nep[b],rsrfarr0[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
    Kcalwid0[b] = calc_k4(alpha_nep[b],rsrfwid0[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
    for a in arange(nalph):
        KMonP0[a,b] = calc_k4(alpharr[a],rsrfarr0[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
        KColP0[a,b] = calc_kc(alpharr[a],rsrfarr0[:,b],apfarr[:,b],nuarr,dnu,nuc[b],alphapip)
        KMonPwid0[a,b] = calc_k4(alpharr[a],rsrfwid0[:,b],apfarr[:,b],nuarr,dnu,nuc[b])
        KColPwid0[a,b] = calc_kc(alpharr[a],rsrfwid0[:,b],apfarr[:,b],nuarr,dnu,nuc[b],alphapip)

wlc0=c/nuc0
wllim0=c/nulim0
print 'Band centres: [%.2f, %.2f,%.2f] GHz ([%.2f,%.2f,%.2f] um)'%(nuc0[0]/1.e9,nuc0[1]/1.e9,nuc0[2]/1.e9,wlc0[0]*1.e6,wlc0[1]*1.e6,wlc0[2]*1.e6)
print 'Low-freq edge: [%.2f, %.2f,%.2f] GHz ([%.2f,%.2f,%.2f] um)'%(nulim0[0,0]/1.e9,nulim0[1,0]/1.e9,nulim0[2,0]/1.e9,wllim0[0,0]*1.e6,wllim0[1,0]*1.e6,wllim0[2,0]*1.e6)
print 'High-freq edge: [%.2f, %.2f,%.2f] GHz ([%.2f,%.2f,%.2f] um)'%(nulim0[0,1]/1.e9,nulim0[1,1]/1.e9,nulim0[2,1]/1.e9,wllim0[0,1]*1.e6,wllim0[1,1]*1.e6,wllim0[2,1]*1.e6)
print 'Band widths: [%.2f, %.2f,%.2f] GHz ([%.2f,%.2f,%.2f] um)'%(rsrflim0[0]/1.e9,rsrflim0[1]/1.e9,rsrflim0[2]/1.e9,1.e6*(wllim0[0,1]-wllim0[0,0]),1.e6*(wllim0[1,1]-wllim0[1,0]),1.e6*(wllim0[2,1]-wllim0[2,0]))
print 'Band resolutions: [%.2f, %.2f, %.2f]'%(rsrfR0[0],rsrfR0[1],rsrfR0[2])
##############################
##  Do random realisations  ##
##############################

nmc=1000
doplot = nmc < 1000
doplot=False
print 'Calculating KMonP and KColP for %d random realisations...'%nmc
rsrfinMC=[]
KpipMC=zeros((3,nmc))
KMonPMC=zeros((nalph,3,nmc))
KColPMC=zeros((nalph,3,nmc))
KcalMC=zeros((3,nmc))

rsrflimMC=zeros((3,nmc))
rsrfRMC=zeros((3,nmc))
nulimMC=zeros((3,nmc,2))
nucMC=zeros((3,nmc))

KpipwidMC=zeros_like(KpipMC)
KMonPwidMC=zeros_like(KMonPMC)
KColPwidMC=zeros_like(KColPMC)
KcalwidMC=zeros_like(KcalMC)
for m in range(nmc):
    rsrfinRand=zeros_like(rsrfin)
    rsrfRand=zeros_like(rsrfarr0)
    for b in band:
        ## Make random RSRF
        rsrfinRand[:,b] = normal(0.,1.,nnuin) * rsrfinsd[:,b] + rsrfin[:,b]
        #rsrfinRand[:,b] = where(rsrfinRand[:,b]<0,rsrf)
        intf=interp1d(nuin,rsrfinRand[:,b],bounds_error=False,fill_value=0)
        rsrfarr=intf(nuarr)
        
        ##calculate band edges
        ilim_med,nulim_med=bandedge(nuarr,rsrfarr,fact=1.,method='median',lim=rsrfmed0[b])
        nulimMC[b,m,:]=nulim_med
        nucMC[b,m]=trapz(nuarr*rsrfarr)/trapz(rsrfarr)
        rsrflimMC[b,m]=nulimMC[b,m,1]-nulimMC[b,m,0]
        rsrfRMC[b,m]=nucMC[b,m]/rsrflimMC[b,m]
        
        rsrfwid=zeros_like(rsrfarr)
        rsrfwid[ilim_med[0]:ilim_med[1]]=1.
        ##calculate K parameters
        KcalMC[b,m] = calc_k4(alpha_nep[b],rsrfarr,apfarr[:,b],nuarr,dnu,nuc[b])/Kcal0[b]
        KcalwidMC[b,m] = calc_k4(alpha_nep[b],rsrfwid,apfarr[:,b],nuarr,dnu,nuc[b])/Kcal0[b]
        KpipMC[b,m] = KcalMC[b,m] * calc_k4(-1,rsrfarr,apfarr[:,b],nuarr,dnu,nuc[b])
        for a in arange(nalph):
            KMonPMC[a,b,m] = KcalMC[b,m] * calc_k4(alpharr[a],rsrfarr,apfarr[:,b],nuarr,dnu,nuc[b])
            KColPMC[a,b,m] = KcalMC[b,m] * calc_kc(alpharr[a],rsrfarr,apfarr[:,b],nuarr,dnu,nuc[b],alphapip)
            KMonPwidMC[a,b,m] = KcalwidMC[b,m] * calc_k4(alpharr[a],rsrfwid,apfarr[:,b],nuarr,dnu,nuc[b])
            KColPwidMC[a,b,m] = KcalwidMC[b,m] * calc_kc(alpharr[a],rsrfwid,apfarr[:,b],nuarr,dnu,nuc[b],alphapip)
            
    rsrfinMC.append(rsrfinRand)

############################
##  Calculate statistics  ##
############################
KMonPmean=zeros((nalph,3))
KMonPsd=zeros((nalph,3))
KColPmean=zeros((nalph,3))
KColPsd=zeros((nalph,3))
KMonPmeanrelmax=zeros(3)
KColPmeanrelmax=zeros(3)
KMonPsdmax=zeros(3)
KColPsdmax=zeros(3)
Kcalmean=zeros(3)
Kcalsd=zeros(3)

rsrflomeanrel=zeros(3)
rsrfhimeanrel=zeros(3)
rsrfwidmeanrel=zeros(3)
rsrflomean=zeros(3)
rsrfhimean=zeros(3)
rsrfwidmean=zeros(3)
rsrflosd=zeros(3)
rsrfhisd=zeros(3)
rsrfwidsd=zeros(3)
nucmean=zeros(3)
nucsd=zeros(3)
rsrfRmean=zeros(3)
rsrfRsd=zeros(3)

KMonPwidmean=zeros((nalph,3))
KMonPwidsd=zeros((nalph,3))
KColPwidmean=zeros((nalph,3))
KColPwidsd=zeros((nalph,3))
KMonPwidmeanrelmax=zeros(3)
KColPwidmeanrelmax=zeros(3)
KMonPwidsdmax=zeros(3)
KColPwidsdmax=zeros(3)
Kcalwidmean=zeros(3)
Kcalwidsd=zeros(3)

KMonPmeanrelrun=zeros((3,nmc))
KColPmeanrelrun=zeros((3,nmc))
KMonPsdrun=zeros((3,nmc))
KColPsdrun=zeros((3,nmc))

for b in band:
    for a in range(nalph):
        KMonPmean[a,b]=mean(KMonPMC[a,b,:])
        KMonPsd[a,b]=std(KMonPMC[a,b,:])
        KColPmean[a,b]=mean(KColPMC[a,b,:])
        KColPsd[a,b]=std(KColPMC[a,b,:])

        KMonPwidmean[a,b]=mean(KMonPwidMC[a,b,:])
        KMonPwidsd[a,b]=std(KMonPwidMC[a,b,:])
        KColPwidmean[a,b]=mean(KColPwidMC[a,b,:])
        KColPwidsd[a,b]=std(KColPwidMC[a,b,:])

    #rsrflomean[b]=mean(nulimMC[b,:,0]/1.e12)
    #rsrfhimean[b]=mean(nulimMC[b,:,1]/1.e12)
    #rsrfwidmean[b]=mean(nulimMC[b,:,1]-nulimMC[b,:,0])/1.e12
    rsrflomean[b]=100.*mean(nulimMC[b,:,0]-nulim0[b,0])/nulim0[b,0]
    rsrfhimean[b]=100.*mean(nulimMC[b,:,1]-nulim0[b,1])/(nulim0[b,1])
    rsrfwidmean[b]=100.*mean((nulimMC[b,:,1]-nulimMC[b,:,0])-(nulim0[b,1]-nulim0[b,0]))/(nulim0[b,1]-nulim0[b,0])
    
    rsrflosd[b]=100.*std(nulimMC[b,:,0])/nulim0[b,0]
    rsrfhisd[b]=100.*std(nulimMC[b,:,1])/nulim0[b,1]
    rsrfwidsd[b]=100.*std(nulimMC[b,:,1]-nulimMC[b,:,0])/(nulim0[b,1]-nulim0[b,0])

    nucmean[b]=100.*(mean(nucMC[b,:])-nuc0[b])/nuc0[b]
    nucsd[b]=100.*std(nucMC[b,:])/nuc0[b]
    rsrfRmean[b]=100.*(mean(rsrfRMC[b,:])-rsrfR0[b])/rsrfR0[b]
    rsrfRsd[b]=100.*std(rsrfRMC[b,:])/rsrfR0[b]

    KMonPmeanrelmax[b]=100.*max((KMonPmean[:,b]-KMonP0[:,b])/KMonP0[:,b])
    KColPmeanrelmax[b]=100.*max((KColPmean[:,b]-KColP0[:,b])/KColP0[:,b])
    KMonPsdmax[b]=100.*max(KMonPsd[:,b]/KMonP0[:,b])
    KColPsdmax[b]=100.*max(KColPsd[:,b]/KColP0[:,b])
    Kcalmean[b] = 100.*mean(KcalMC[b,:]-1.)
    Kcalsd[b] = 100.*std(KcalMC[b,:])

    KMonPwidmeanrelmax[b]=100.*max((KMonPwidmean[:,b]-KMonPwid0[:,b])/KMonPwid0[:,b])
    KColPwidmeanrelmax[b]=100.*max((KColPwidmean[:,b]-KColPwid0[:,b])/KColPwid0[:,b])
    KMonPwidsdmax[b]=100.*max(KMonPwidsd[:,b]/KMonPwid0[:,b])
    KColPwidsdmax[b]=100.*max(KColPwidsd[:,b]/KColPwid0[:,b])
    Kcalwidmean[b] = 100.*mean(KcalwidMC[b,:]-1.)
    Kcalwidsd[b] = 100.*std(KcalwidMC[b,:])
    
    for m in range(nmc):
        KMonPmeanrelrun[b,m]=(mean(KMonPMC[nalph-1,b,0:m]) - KMonP0[nalph-1,b])/KMonP0[nalph-1,b]
        KColPmeanrelrun[b,m]=(mean(KColPMC[nalph-1,b,0:m]) - KColP0[nalph-1,b])/KColP0[nalph-1,b]
        KMonPsdrun[b,m]=std(KMonPMC[nalph-1,b,0:m])/KMonP0[nalph-1,b]
        KColPsdrun[b,m]=std(KColPMC[nalph-1,b,0:m])/KColP0[nalph-1,b]

print '---'
print 'Kcal mean offset (%%): [%.3g , %.3g , %.3g]'%(Kcalmean[0],Kcalmean[1],Kcalmean[2])
print 'Kcal standard deviation (%%): [%.3g , %.3g , %.3g]'%(Kcalsd[0],Kcalsd[1],Kcalsd[2])
print 'Kcal.KMonP max mean offset (%%): [%.3g , %.3g , %.3g]'%(KMonPmeanrelmax[0],KMonPmeanrelmax[1],KMonPmeanrelmax[2])
print 'Kcal.KColP max mean offset (%%): [%.3g , %.3g , %.3g]'%(KColPmeanrelmax[0],KColPmeanrelmax[1],KColPmeanrelmax[2])
print 'Kcal.KMonP max standard deviation (%%): [%.3g , %.3g , %.3g]'%(KMonPsdmax[0],KMonPsdmax[1],KMonPsdmax[2])
print 'Kcal.KColP max standard deviation (%%): [%.3g , %.3g , %.3g]'%(KColPsdmax[0],KColPsdmax[1],KColPsdmax[2])
print '---'
print 'RSRF low-edge mean offset (%%): [%.3g , %.3g , %.3g]'%(rsrflomean[0],rsrflomean[1],rsrflomean[2])
print 'RSRF low-edge standard deviation (%%): [%.3g , %.3g , %.3g]'%(rsrflosd[0],rsrflosd[1],rsrflosd[2])
print 'RSRF high-edge mean offset (%%): [%.3g , %.3g , %.3g]'%(rsrfhimean[0],rsrfhimean[1],rsrfhimean[2])
print 'RSRF high-edge standard deviation (%%): [%.3g , %.3g , %.3g]'%(rsrfhisd[0],rsrfhisd[1],rsrfhisd[2])
print 'RSRF width mean offset (%%): [%.3g , %.3g , %.3g]'%(rsrfwidmean[0],rsrfwidmean[1],rsrfwidmean[2])
print 'RSRF width standard deviation (%%): [%.3g , %.3g , %.3g]'%(rsrfwidsd[0],rsrfwidsd[1],rsrfwidsd[2])
print 'RSRF centre mean offset (%%): [%.3g , %.3g , %.3g]'%(nucmean[0],nucmean[1],nucmean[2])
print 'RSRF centre standard deviation(%%): [%.3g , %.3g , %.3g]'%(nucsd[0],nucsd[1],nucsd[2])
print 'RSRF relative width mean offset (%%): [%.3g , %.3g , %.3g]'%(rsrfRmean[0],rsrfRmean[1],rsrfRmean[2])
print 'RSRF relative width standard deviation (%%): [%.3g , %.3g , %.3g]'%(rsrfRsd[0],rsrfRsd[1],rsrfRsd[2])
print '---'
print 'Kcal mean offset wid (%%): [%.3g , %.3g , %.3g]'%(Kcalwidmean[0],Kcalwidmean[1],Kcalwidmean[2])
print 'Kcal standard deviation wid (%%): [%.3g , %.3g , %.3g]'%(Kcalwidsd[0],Kcalwidsd[1],Kcalwidsd[2])
print 'Kcal.KMonP max mean offset wid (%%): [%.3g , %.3g , %.3g]'%(KMonPwidmeanrelmax[0],KMonPwidmeanrelmax[1],KMonPwidmeanrelmax[2])
print 'Kcal.KColP max mean offset wid (%%): [%.3g , %.3g , %.3g]'%(KColPwidmeanrelmax[0],KColPwidmeanrelmax[1],KColPwidmeanrelmax[2])
print 'Kcal.KMonP max standard deviation wid (%%): [%.3g , %.3g , %.3g]'%(KMonPwidsdmax[0],KMonPwidsdmax[1],KMonPwidsdmax[2])
print 'Kcal.KColP max standard deviation wid (%%): [%.3g , %.3g , %.3g]'%(KColPwidsdmax[0],KColPwidsdmax[1],KColPwidsdmax[2])

####################
##  Plot figures  ##
####################

if doplot:
    print 'Plotting figures...'
    npK=5
    if npK > nmc:
        npK = nmc
    npR=1
    if npR > nmc:
        npR = nmc
    
    plot.figure(1,figsize=(12,6))
    plot.clf()
    cols=['b','g','r']
    for b in band:
        for m in range(nmc):
            plot.plot(nuin/1.e12,rsrfinMC[m][:,b],c=(0.7,0.7,0.7,0.01),lw=1)
        for p in range(npR):
            plot.plot(nuin/1.e12,rsrfinMC[p][:,b],c='k',lw=1)
        plot.plot(nuin/1.e12,rsrfin[:,b],c=cols[b],ls='-')
        #plot.plot(nuin/1.e12,rsrfin[:,b],'kx')
        plot.plot(nuin/1.e12,rsrfin[:,b]-rsrfinsd[:,b],c=cols[b],ls='--')
        plot.plot(nuin/1.e12,rsrfin[:,b]+rsrfinsd[:,b],c=cols[b],ls='--')
        plot.annotate(bandStr[b],(nuc[b]/1.e12,0.1),color=cols[b],weight='bold',size=24,family='sans-serif',ha='center',va='center')
    plot.xlim(nulimin_2[2,0]/1.e12,nulimin_2[0,1]/1.e12)
    plot.ylim(0,1)
    plot.xlabel('Frequency, THz')
    plot.ylabel('Spectral response function')
    plot.draw()
    plot.savefig('../RSRF/Figs/RSRF_%sErr_MC%d.eps'%(errStr,nmc),transparent=False)
    plot.savefig('../RSRF/Figs/RSRF_%sErr_MC%d.png'%(errStr,nmc),transparent=True)
    
    plot.figure(2,figsize=(12,6))
    plot.clf()
    for b in band:
        plot.subplot(1,3,b+1)
        for m in range(nmc):
            plot.plot(alpharr,KMonPMC[:,b,m],'0.7',lw=1)
        plot.plot(alpharr,KMonP0[:,b],'k-')
        #plot.axhline(1.,c='k',ls=':')
        plot.xlabel(r'Spectral Index, $\alpha$')
        plot.ylim(0.75,1.05)
        if b==0:
            plot.ylabel(r'$K_\mathrm{Kcal}(\alpha) K_\mathrm{MonP}(\alpha)$')
            plot.yticks(arange(0.75,1.1,0.05))
        else:
            plot.yticks(arange(0.75,1.1,0.05),('','','','','','',''))
        plot.title('%s'%(bandStr[b]))
    plot.draw()
    plot.savefig('../RSRF/Figs/KMonP_%sErr_MC%d.eps'%(errStr,nmc),transparent=False)
    plot.savefig('../RSRF/Figs/KMonP_%sErr_MC%d.png'%(errStr,nmc),transparent=True)
        
    plot.figure(3,figsize=(12,6))
    plot.clf()
    for b in band:
        plot.subplot(1,3,b+1)
        for m in range(nmc):
            plot.plot(alpharr,KColPMC[:,b,m],'0.7',lw=1)
        plot.plot(alpharr,KColP0[:,b],'k-')
        #plot.axhline(1.,c='k',ls=':')
        plot.xlabel(r'Spectral Index, $\alpha$')
        plot.ylim(0.75,1.05)
        if b==0:
            plot.ylabel(r'$K_\mathrm{Cal}(\alpha) K_\mathrm{ColP}(\alpha)$')
            plot.yticks(arange(0.75,1.1,0.05))
        else:
            plot.yticks(arange(0.75,1.1,0.05),('','','','','','',''))
        plot.title('%s'%(bandStr[b]))
    plot.draw()
    plot.savefig('../RSRF/Figs/KColP_%sErr_MC%d.eps'%(errStr,nmc),transparent=False)
    plot.savefig('../RSRF/Figs/KColP_%sErr_MC%d.png'%(errStr,nmc),transparent=True)
    
    plot.figure(4,figsize=(12,6))
    plot.clf()
    for b in band:
        plot.subplot(1,3,b+1)
        for m in range(nmc):
            plot.plot(alpharr,(KMonPMC[:,b,m]-KMonP0[:,b])/KMonP0[:,b],'0.7',lw=1.)
        for p in range(npK):
            plot.plot(alpharr,(KMonPMC[:,b,p]-KMonP0[:,b])/KMonP0[:,b],'b-',lw=1.)
        #plot.axhline(0.,c='k')
        plot.plot(alpharr,(KMonPmean[:,b]-KMonP0[:,b])/KMonP0[:,b],'k-')
        plot.plot(alpharr,((KMonPmean[:,b]-KMonPsd[:,b])-KMonP0[:,b])/KMonP0[:,b],'k--')
        plot.plot(alpharr,((KMonPmean[:,b]+KMonPsd[:,b])-KMonP0[:,b])/KMonP0[:,b],'k--')
        plot.ylim(-0.02,0.02)
        plot.xlabel(r'Spectral Index, $\alpha$')
        if b==0:
            plot.ylabel(r'$(K_\mathrm{Cal}K_\mathrm{MonP} - K_\mathrm{MonP}^\mathrm{true}) / K_\mathrm{MonP}^\mathrm{true}$')
            plot.yticks(arange(-0.02,0.025,0.005))
        else:
            plot.yticks(arange(-0.02,0.025,0.005),('','','','','','','','',''))
        plot.title('%s'%(bandStr[b]))
    plot.draw()
    plot.savefig('../RSRF/Figs/KMonPrel_%sErr_MC%d.eps'%(errStr,nmc),transparent=False)
    plot.savefig('../RSRF/Figs/KMonPrel_%sErr_MC%d.png'%(errStr,nmc),transparent=True)
        
    plot.figure(5,figsize=(12,6))
    plot.clf()
    for b in band:
        plot.subplot(1,3,b+1)
        for m in range(nmc):
            plot.plot(alpharr,(KColPMC[:,b,m]-KColP0[:,b])/KColP0[:,b],'0.7',lw=1)
        for p in range(npK):
            plot.plot(alpharr,(KColPMC[:,b,p]-KColP0[:,b])/KColP0[:,b],'b-',lw=1)
        #plot.axhline(0.,c='k')
        plot.plot(alpharr,(KColPmean[:,b]-KColP0[:,b])/KColP0[:,b],'k-')
        plot.plot(alpharr,((KColPmean[:,b]-KColPsd[:,b])-KColP0[:,b])/KColP0[:,b],'k--')
        plot.plot(alpharr,((KColPmean[:,b]+KColPsd[:,b])-KColP0[:,b])/KColP0[:,b],'k--')
        plot.ylim(-0.02,0.02)
        plot.xlabel(r'Spectral Index, $\alpha$')
        if b==0:
            plot.ylabel(r'$(K_\mathrm{Cal}K_\mathrm{ColP} - K_\mathrm{ColP}^\mathrm{true}) / K_\mathrm{ColP}^\mathrm{true}$')
            plot.yticks(arange(-0.02,0.025,0.005))
        else:
            plot.yticks(arange(-0.02,0.025,0.005),('','','','','','','','',''))
        plot.title('%s'%(bandStr[b]))
    plot.draw()
    plot.savefig('../RSRF/Figs/KColPrel_%sErr_MC%d.eps'%(errStr,nmc),transparent=False)
    plot.savefig('../RSRF/Figs/KColPrel_%sErr_MC%d.png'%(errStr,nmc),transparent=True)
    
    plot.figure(6,figsize=(12,6))
    plot.clf()
    #cols=['b','g','r']
    cols=['k','k','k']
    for b in band:
        plot.subplot(1,3,b+1)
        #plot.plot(arange(nmc)+1,KMonPmeanrelrun[b,:],c=cols[b],ls='-')
        plot.plot(arange(nmc)+1,KMonPsdrun[b,:],c=cols[b],ls='-')
        #plot.plot(arange(nmc)+1,KColPmeanrelrun[b,:],c=cols[b],ls='--')
        plot.plot(arange(nmc)+1,KColPsdrun[b,:],c=cols[b],ls='--')
        plot.xscale('log')
        if b==1:
            plot.xlabel('Number of realisations')
        plot.ylim(0,0.02)
        if b==0:
            plot.ylabel('Standard Deviation')
            plot.yticks(arange(0,0.022,0.002))
        else:
            plot.yticks(arange(0,0.022,0.002),('','','','','','','','','','',''))
        plot.title(bandStr[b])
    plot.draw()
    plot.savefig('../RSRF/Figs/Krun_%sErr_MC%d.eps'%(errStr,nmc),transparent=False)
    plot.savefig('../RSRF/Figs/Krun_%sErr_MC%d.png'%(errStr,nmc),transparent=True)
    
    if errStr=='Meas':
        plot.figure(7,figsize=(12,6))
        plot.clf()
        cols=['b','g','r']
        for b in band:
            plot.plot(nuin/1.e12,specin[:,b],c=cols[b],ls='-')
            plot.plot(nuin/1.e12,specinnofilt[:,b],c=cols[b],ls=':')
            plot.plot(nuin/1.e12,sdin[:,b],c=cols[b],ls='--')
            plot.annotate(bandStr[b],(nuc[b]/1.e12,-0.15),color=cols[b],weight='bold',size=24,family='sans-serif',ha='center',va='center')
        plot.xlabel('Frequency, THz')
        plot.ylabel('Measured Spectrum (arb. units)')
        
        plot.draw()
        plot.savefig('../RSRF/Figs/MeasSpec.eps',transparent=False)
        plot.savefig('../RSRF/Figs/MeasSpec.png',transparent=False)
        
        plot.figure(8,figsize=(12,6))
        plot.clf()
        for b in band:
            plot.plot(nuin/1.e12,relerr[:,b],c=cols[b],ls='-')
            plot.annotate(bandStr[b],(nuc[b]/1.e12,0.015),color=cols[b],weight='bold',size=24,family='sans-serif',ha='center',va='center')
            plot.axvline(nulim[b,0]/1.e12,c=cols[b],ls=':')
            plot.axvline(nulim[b,1]/1.e12,c=cols[b],ls=':')
        plot.xlabel('Frequency, THz')
        plot.ylabel('Relative error')
        plot.yscale('log')
        plot.draw()
        plot.savefig('../RSRF/Figs/RelErr.eps',transparent=False)
        plot.savefig('../RSRF/Figs/RelErr.png',transparent=False)

plot.figure(8,figsize=(12,6))
plot.clf()
plot.subplot(2,3,1)
cols=['b','g','r']
for b in band:
    plot.hist((nulimMC[b,:,0]-nulim0[b,0])/nulim0[b,0],bins=10,align='mid',histtype='step',color=cols[b])
plot.xlabel(r'$(\nu-\nu_0)/\nu_0$')
plot.axvline(0,c='k',ls=':')
plot.title('Low edge')

plot.subplot(2,3,2)
cols=['b','g','r']
for b in band:
    plot.hist((nulimMC[b,:,1]-nulim0[b,1])/nulim0[b,1],bins=10,align='mid',histtype='step',color=cols[b])
plot.xlabel(r'$(\nu-\nu_0)/\nu_0$')
plot.axvline(0,c='k',ls=':')
plot.title('High edge')

plot.subplot(2,3,3)
cols=['b','g','r']
for b in band:
    plot.hist((rsrflimMC[b,:]-rsrflim0[b])/rsrflim0[b],bins=10,align='mid',histtype='step',color=cols[b])
plot.xlabel(r'$(\Delta\nu-\Delta\nu_0)/\Delta\nu_0$')
plot.axvline(0,c='k',ls=':')
plot.title('Band width')

plot.subplot(2,3,4)
cols=['b','g','r']
for b in band:
    plot.hist((nucMC[b,:]-nuc0[b])/nuc0[b],bins=10,align='mid',histtype='step',color=cols[b])
plot.xlabel(r'$(\nu_c-\nu_{c,0})/\nu_{c,0}$')
plot.axvline(0,c='k',ls=':')
plot.title('Band centre')

plot.subplot(2,3,5)
cols=['b','g','r']
for b in band:
    plot.hist((rsrfRMC[b,:]-rsrfR0[b])/rsrfR0[b],bins=10,align='mid',histtype='step',color=cols[b])
plot.xlabel(r'$(R-R_0)/R_0$')
plot.axvline(0,c='k',ls=':')
plot.title('Relative Bandwidth')

plot.savefig('../RSRF/Figs/Bandlim_%sErr_MC%d.png'%(errStr,nmc),transparent=True)
plot.savefig('../RSRF/Figs/Bandlim_%sErr_MC%d.eps'%(errStr,nmc),transparent=True)
plot.draw()

plot.figure(9,figsize=(12,6))
plot.clf()
cols=('b','g','r')
for b in band:
    plot.plot(nuarr/1.e12,rsrfarr0[:,b],c=cols[b],ls='--')
    plot.plot([nulim0[b,0]/1.e12,nulim0[b,0]/1.e12],[0,rsrfmed0[b]],c=cols[b],ls='-')
    plot.plot([nulim0[b,1]/1.e12,nulim0[b,1]/1.e12],[0,rsrfmed0[b]],c=cols[b],ls='-')
    plot.plot(nulim0[b,:]/1.e12,[rsrfmed0[b],rsrfmed0[b]],c=cols[b],ls='-')
    plot.annotate(bandStr[b],(nuc[b]/1.e12,0.1),color=cols[b],weight='bold',size=24,family='sans-serif',ha='center',va='center')
    #plot.plot([nulim0[b,0]*(1-rsrflosd/100.)/1.e12,nulim0[b,0]*(1.-rsrflosd/100.)/1.e12],[0,rsrfmed0[b]],c=cols[b],ls=':')
    #plot.plot([nulim0[b,0]*(1.+rsrflosd/100.)/1.e12,nulim0[b,0]*(1.+rsrflosd/100.)/1.e12],[0,rsrfmed0[b]],c=cols[b],ls=':')
    #plot.plot([nulim0[b,1]*(1-rsrfhisd/100.)/1.e12,nulim0[b,1]*(1.-rsrfhisd/100.)/1.e12],[0,rsrfmed0[b]],c=cols[b],ls=':')
    #plot.plot([nulim0[b,1]*(1.+rsrfhisd/100.)/1.e12,nulim0[b,1]*(1.+rsrfhisd/100.)/1.e12],[0,rsrfmed0[b]],c=cols[b],ls=':')
plot.ylim(0,1)
plot.ylabel('Spectral Response Function')
plot.xlabel('Frequency (THz)')
plot.savefig('../RSRF/Figs/RSRF_Tophat.png',transparent=True)
plot.savefig('../RSRF/Figs/RSRF_TopHat.eps',transparent=True)

plot.draw()        


plot.show()

#if __name__=="__main__":
#    maincode()