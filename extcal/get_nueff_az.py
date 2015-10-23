from numpy import *
import scipy
from beam import *
from rsrf import *
import scipy
import scipy.integrate as integrate
import scipy.interpolate as interp


c=299792458. #speed of light
#set band centres
wlc=(250.e-6,350.e-6,500.e-6)
wlc=array(wlc)
nuc=c/wlc
band=arange(3)

#Beamtype options:
#  G=Gaussian
#  E=Elliptical Gaussian
#  M=Measured beams
#  T=Theoretical (modelled) beams)
beamtype='M'

#Set beam grid options
bgrid=1. #grid size in arcsec
bwid=2000 #full width of map in arcsec
#bzlim=1.e-8
bzlim=None
bzval=0.
npx=bwid/bgrid
beams=zeros((npx,npx,3))
for b in band:
    beamx,beamy,beams[:,:,b]=getbeam(beamtype,b,bgrid,bwid,bzlim=bzlim,bzval=bzval)

#Set up nu array
numin=300.e9 #300GHz
numax=1800.e9 #1790 GHz
dnu=1.e9 #1GHz step
nnu=1+(numax - numin)/dnu
nul=range(int(nnu))
nuarr=array(nul)*dnu + numin
print nnu

a_arr=array([[1.,1.01],[1.,1.02],[1.,1.03]])
nu0=nuc*a_arr[:,0]

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
alpha_nep=array((1.26,1.39,1.45))

radarr=arange(1000)
areameas=zeros(3)
for b in band:
    print 'Band %d'%b
    areameas[b]=measbeam(beamx,beamy,beams[:,:,b],radarr,nuarr,rsrfarr[:,b],nu0[b],alpha_nep[b])

#nrad=radarr.size
#beam_sm=zeros((nrad,3))
#beammeas=zeros((nrad,3))
#beam_sm_nu=zeros(nnu)
#for b in band:
#    print 'Band',b
#    beam_sm[:,b]=beam_azsm(beamx,beamy,beams[:,:,b],radarr)
#    beam_sm_int=interp.interp1d(radarr,beam_sm[:,b])
#    denom_arg = rsrfarr[:,b] * nuarr**alph_nep[b]
#    for r in arange(nrad):
#        r_nu=radarr[r]*nuarr/nu0[b]
#        inr=scipy.where(r_nu < max(radarr))
#        outr=scipy.where(r_nu >= max(radarr))
#        beam_sm_nu[inr]=beam_sm_int(r_nu[inr])
#        beam_sm_nu[outr]=0.
#        num_arg = rsrfarr[:,b] * nuarr**alph_nep[b] * beam_sm_nu
#        beammeas[r,b]=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)
#    #calculate beam area
#    areameas[b]=integrate.trapz(beammeas[:,b]*2.*pi*radarr,radarr)


print 'Meas areas: [%.2f, %.2f, %.2f]' % (areameas[0],areameas[1],areameas[2])
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
    area[b,0]=beamarea(beamx,beamy,beams[:,:,b],brad)
    area[b,1]=sum(beams[:,:,b])

print 'Beam Areas (350"): [%.2f, %.2f, %.2f]' % (area[0,0],area[1,0],area[2,0])
print 'Beam Areas (sum): [%.2f, %.2f, %.2f]' % (area[0,1],area[1,1],area[2,1])