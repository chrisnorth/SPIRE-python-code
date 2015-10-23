from numpy import *
import scipy
from beam import *
from rsrf import *
import scipy.integrate as integrate

#Beamtype options:
#  G=Gaussian
#  E=Elliptical Gaussian
#  M=Measured beams
#  T=Theoretical (modelled) beams)
beamtype='T'

#Set beam grid options
bgrid=1. #grid size in arcsec
bwid=701 #needs to be odd to ensure that there is a (0,0) element
brad=350
bzlim=1.e-8
bzval=0.
npx=bwid/bgrid
beams=zeros((npx,npx,3))
for b in band:
    beamx,beamy,beams[:,:,b]=getbeam(beamtype,b,bgrid,bwid,bzlim=bzlim,bzval=bzval)

c=299792458. #speed of light
#set band centres
wlc=(250.e-6,350.e-6,500.e-6)
wlc=array(wlc)
nuc=c/wlc
band=arange(3)

#Set up nu array
numin=300.e9 #300GHz
numax=1800.e9 #1790 GHz
dnu=1.e9 #1GHz step
nnu=1+(numax - numin)/dnu
nul=range(int(nnu))
nuarr=array(nul)*dnu + numin
print nnu
nu0=nuc

#get RSRF
rsrftype='M'
rsrfarr=zeros((nnu,3))
for b in band:
    rsrfarr[:,b]=getrsrf(rsrftype,b,nuarr)

#Aperture Efficiency options:
apftype='R'
apfarr=zeros((nnu,3))
for b in band:
    apfarr[:,b]=getapf(apftype,b,nuarr)

#Spectral index of neptune
alph_nep=array((1.26,1.39,1.45))

beammeas=zeros((npx,npx,3))
areameas=zeros(3)

for b in band:
    print 'Band %d'%b
    intbeam=scipy.interpolate.RectBivariateSpline(beamx[:,0],beamy[0,:],beams[:,:,b])
    denom_arg = rsrfarr[:,b] * apfarr[:,b] * nuarr**alph_nep[b]
#    for n in arange(nnu):
#        print n,'of',nnu
#        beamx_nu=beamx[:,0]*nuarr[n]/nu0[b]
#        beamy_nu=beamy[:,0]*nuarr[n]/nu0[b]
#        beam_xynu[:,:,n]=intbeam(beamx_nu,beamy_nu)
#
#    for x in arange(npx):
#        for y in arange(npx):
#            print x,y,'of',npx
#            num_arg=rsrfarr[:,b] * apfarr[:,b] * nuarr**alph_nep[b] * beam_xynu[x,y,:]
#            beammeas[x,y,b]=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)
    
    for x in arange(npx):
        print 100.*x/npx
        for y in arange(npx):
            #calculate beam at that position over frequency range
            x_nu=beamx[x,y]*nuarr/nu0[b]
            y_nu=beamy[x,y]*nuarr/nu0[b]
            beam_xy=intbeam.ev(x_nu,y_nu)

            #calculate measured beam            
            num_arg = rsrfarr[:,b] * apfarr[:,b] * nuarr**alph_nep[b] * beam_xy
            beammeas[x,y,b]=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)

    #calculate beam area
    areameas[b]=sum(beammeas[:,:,b])

print 'Areas: [%.2f, %.2f, %.2f]' % (areameas[0],areameas[1],areameas[2])
##Will calculate this empirically eventually
#a_bm=array([1.0184,1.0164,1.0204])
#nu0_bm=nuc*a_bm
nuL_bm=array([1033.,741.,491.])*1.e9
nuU_bm=array([1418.,1008.,732.])*1.e9

#measured areas from
area_mjg=[435.8,776.1,1623.9]
brad=350.
area=zeros((3,2))
for b in band:
    area[b,0]=beamarea(beamx,beamy,beams[:,:,b],brad=brad)
    area[b,1]=sum(beams[:,:,b])
print 'Beam areas: [%.2f, %.2f, %.2f]' % (area[0,0],area[1,0],area[2,0])
print 'Beam areas: [%.2f, %.2f, %.2f]' % (area[0,1],area[1,1],area[2,1])

##Set frequency arrays for beam area calculations
dnu_beam=10.e9 #freq spacing
nnu_bm=zeros((3))
for b in band:
    nnu_bm[b]=int((nuU_bm[b]-nuL_bm[b])/dnu_beam)

print 'Computing area (PSW)'
area_bm0=zeros(nnu_bm[0])
nu_bm0=arange(nnu_bm[0])*dnu_beam + nuL_bm[0]
for n in arange(nnu_bm[0]):
    regrid_bm=beam_regrid(beamx,beamy,beams[:,:,0],nu0_bm[0],nu_bm0[n])
    area_bm0[n]=beamarea(beamx,beamy,regrid_bm,brad=brad)

a0=abs(area_bm0-area_mjg[0])
nu0=nu_bm0[where(a0 == min(a0))]
afact
print 'PSW Frequency: %.2f a=%.4f' %(nu0_0[0]/1.e9,nu0_1[0]/1.e9)

print 'Computing area (PMW)'
area_bm1=zeros((nnu_bm[1],2))
nu_bm1=arange(nnu_bm[1])*dnu_beam + nuL_bm[1]
for n in arange(nnu_bm[1]):
    regrid_bm=beam_regrid(beamx,beamy,beams[:,:,1],nu0_bm[1],nu_bm1[n])
    area_bm1[n,0]=beamarea(beamx,beamy,regrid_bm,brad=brad)
    area_bm1[n,1]=sum(regrid_bm)
a1_0=abs(area_bm1[:,0]-area_mjg[1])
a1_1=abs(area_bm1[:,1]-area_mjg[1])
nu1_0=nu_bm1[where(a1_0 == min(a1_0))]
nu1_1=nu_bm1[where(a1_1 == min(a1_1))]
print 'Frequency: %.2f / %.2f GHz' %(nu1_0[0]/1.e9,nu1_1[0]/1.e9)    
    
print 'Computing area (PLW)'
area_bm2=zeros((nnu_bm[2],2))
nu_bm2=arange(nnu_bm[2])*dnu_beam + nuL_bm[2]
for n in arange(nnu_bm[2]):
    regrid_bm=beam_regrid(beamx,beamy,beams[:,:,2],nu0_bm[2],nu_bm2[n])
    area_bm2[n,0]=beamarea(beamx,beamy,regrid_bm,brad=brad)
    area_bm2[n,1]=sum(regrid_bm)
a2_0=abs(area_bm2[:,0]-area_mjg[2])
a2_1=abs(area_bm2[:,1]-area_mjg[2])
nu2_0=nu_bm2[where(a2_0 == min(a2_0))]
nu2_1=nu_bm2[where(a2_1 == min(a2_1))]
print 'Frequency: %.2f / %.2f GHz' %(nu2_0[0]/1.e9,nu2_1[0]/1.e9)