##Calculate band centre and width of bands from RSRF

from numpy import array,zeros,arange,mean,median
from rsrf import getrsrf,bandedge
import matplotlib.pyplot as plot

c=299792458. #speed of light
#set band centres
wlc=(250.e-6,350.e-6,500.e-6)
wlc=array(wlc)
nuc=c/wlc
#print nuc/1.e12

band=arange(3)
bandStr=['PSW','PMW','PLW']

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
rsrf=zeros((nnu,3))
rsrftype='M'
nulim_max=zeros((3,2))
ilim_max=zeros((3,2))
nulim_med=zeros((3,2))
ilim_med=zeros((3,2))
rsrfmean=zeros(3)
rsrfmed=zeros(3)
for b in band:
    rsrf[:,b]=getrsrf(rsrftype,b,nuarr)
    ilim_max[b,:],nulim_max[b,:]=bandedge(nuarr,rsrf[:,b],fact=1.)
    ilim_med[b,:],nulim_med[b,:]=bandedge(nuarr,rsrf[:,b],fact=1.,method='median')
    rsrfmean[b]=mean(rsrf[ilim_max[b,0]:ilim_max[b,1],b])
    rsrfmed[b]=median(rsrf[ilim_max[b,0]:ilim_max[b,1],b])
    print 'band %d (max): '%(b),ilim_max[b,:],nulim_max[b,:]/1.e12,rsrfmean[b]
    print 'band %d (med): '%(b),ilim_med[b,:],nulim_med[b,:]/1.e12,rsrfmed[b]
plot.figure(1,figsize=(12,6))
plot.clf()
for b in band:
    plot.subplot(1,3,b+1)
    plot.hist(rsrf[ilim_max[b,0]:ilim_max[b,1],b],bins=50,align='mid',normed=True,histtype='stepfilled',color='b')
    plot.axvline(rsrfmean[b],c='r')
    plot.axvline(rsrfmed[b],c='g')
plot.draw()

plot.figure(2,figsize=(12,6))
plot.clf()
cols=('b','g','r')
for b in band:
    plot.plot(nuarr/1.e12,rsrf[:,b],c=cols[b],ls='-')
    plot.plot([nulim_med[b,0]/1.e12,nulim_med[b,0]/1.e12],[0,rsrfmed[b]],c=cols[b],ls=':')
    plot.plot([nulim_med[b,1]/1.e12,nulim_med[b,1]/1.e12],[0,rsrfmed[b]],c=cols[b],ls=':')
    plot.plot(nulim_med[b,:]/1.e12,[rsrfmed[b],rsrfmed[b]],c=cols[b],ls=':')

    plot.plot([nulim_max[b,0]/1.e12,nulim_max[b,0]/1.e12],[0,rsrfmean[b]],c=cols[b],ls='--')
    plot.plot([nulim_max[b,1]/1.e12,nulim_max[b,1]/1.e12],[0,rsrfmean[b]],c=cols[b],ls='--')
    plot.plot(nulim_max[b,:]/1.e12,[rsrfmean[b],rsrfmean[b]],c=cols[b],ls='--')

plot.draw()
plot.show()