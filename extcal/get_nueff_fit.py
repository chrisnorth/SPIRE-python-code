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

brad=1000.

area=zeros(3)
for b in band:
    area[b]=beamarea(beamx,beamy,beams[:,:,b],brad=brad)

print 'Beam Areas (brad): [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])

#Set up nu array
numin=300.e9 #300GHz
numax=1800.e9 #1790 GHz
dnu=1.e9 #1GHz step
nnu=1+(numax - numin)/dnu
nul=range(int(nnu))
nuarr=array(nul)*dnu + numin
print nnu

a_init=array([[1.,1.01],[1.,1.02],[1.,1.03]])

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

aprec=0.0001 #relative precision required on area
maxit=10.
radarr=arange(brad)
a_fin=zeros(3)
ameas_fin=zeros(3)
arel_fin=zeros(3)
for b in band:
    print 'Band %d'%b
    a_arr=a_init[b,:]
    ameas=zeros(2)
    beam_sm=beam_azsm(beamx,beamy,beams[:,:,b],radarr)
    ameas[0]=measbeam(radarr,beam_sm,nuarr,rsrfarr[:,b],nuc[b]*a_arr[0],alpha_nep[b])
    ameas[1]=measbeam(radarr,beam_sm,nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b])
    adiff=ameas-area[b]
    arel=(ameas-area[b])/area[b]
    print 'a=%.4f, A=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[0],ameas[0],adiff[0],arel[0])
    print 'a=%.4f, A=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
    it=1
    done=False
    while done==False:
        grad=(adiff[1]-adiff[0])/(a_arr[1]-a_arr[0])
        a_arrnew=a_arr[1] - adiff[1]/grad
        a_arr[0]=a_arr[1]
        ameas[0]=ameas[1]
        a_arr[1]=a_arrnew
        ameas[1]=measbeam(radarr,beam_sm,nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b])
        adiff=ameas-area[b]        
        arel=(ameas-area[b])/area[b]
        print 'a=%.4f, Area=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
        it=it+1
        if abs(arel[1]) < aprec:
            done=True
        elif it >= maxit:
            done=True
    a_fin[b]=a_arr[1]
    ameas_fin[b]=ameas[1]
    arel_fin[b]=arel[1]
    print 'Final a=%.4f, Area=%.2f (Rel. error=%.6f)' % (a_fin[b],ameas_fin[b],arel_fin[b])
    print 'Nu(eff)=%.2f GHz'%(nuc[b]*a_fin[b]/1.e9)

print 'Meas areas: [%.2f, %.2f, %.2f]' % (ameas_fin[0],ameas_fin[1],ameas_fin[2])
print 'True areas: [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])
##Will calculate this empirically eventually
#a_bm=array([1.0184,1.0164,1.0204])
#nu0_bm=nuc*a_bm
nuL_bm=array([1033.,741.,491.])*1.e9
nuU_bm=array([1418.,1008.,732.])*1.e9

#measured areas from
area_mjg=[435.8,776.1,1623.9]