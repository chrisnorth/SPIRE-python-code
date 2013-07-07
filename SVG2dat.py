###Read in spectral data from SVG files and calculate relative errors

import matplotlib.pyplot as plot
from rsrf import readsvg
from numpy import mean,std,zeros,sum,sqrt,exp,min
from scipy import where
from scipy.io.idl import readsav
from scipy.interpolate import interp1d
import sys

bandStr='PSW'
filemeas='../RSRF/Locke/PHOTrsrf_%s_coadded-meas.dat'%bandStr
filesd='../RSRF/Locke/PHOTrsrf_%s_coadded-sd.dat'%bandStr
filesav='../RSRF/Locke/%s_allIFGMS_SVRvars.sav'%bandStr

##Read in IDL SAV file
savvar=readsav(filesav,verbose=False)
#sys.exit()

##Read in SVG dat files
wn,meanspec,detspec = readsvg(filemeas,verbose=False,dofilt=True)
wnnofilt,meanspecnofilt,detspecnofilt = readsvg(filemeas,verbose=False,dofilt=False)
wnsd,meansd,detsd = readsvg(filesd,wnin=wn,verbose=False,dofilt=True)

c=299792458. #speed of light
freq=wn*c*1.e2

nspec=len(detspec[0,:])
nsd=len(detsd[0,:])
print nspec,' spec ; ',nsd,' sd'
nc=len(wn)
newsd=zeros(nc)
quadsd=zeros(nc)
newmean=zeros(nc)
weights=1./detsd**2
#weightspec=detspec/detsd**2
for c in range(nc):
    newmean[c]=sum(weights[c,:]*detspec[c,:])/sum(weights[c,:])
    newsd[c]=std(detspec[c,:])
    quadsd[c]=sqrt(sum(weights[c,:]**2 * detsd[c,:]**2))/sum(weights[c,:])
print sum(weights[c,:])




###interpolate IDL files
atmint=interp1d(savvar.wn,savvar.atm)
#atm=1./exp(atmint(wn))
atm=atmint(wn)

hbbint=interp1d(savvar.wn,savvar.hbb)
hbb=hbbint(wn)
hbb=hbb/mean(hbb)

respint=interp1d(savvar.wn,savvar.resp)
resp=respint(wn)
resp=resp/mean(resp)

savspecint=interp1d(savvar.wn,savvar.spec)
savspec=savspecint(wn)

savrsrfint=interp1d(savvar.wn,savvar.rsrf)
savrsrf=savrsrfint(wn)

#cal=savrsrf/savspec
#calspec=meanspec*cal
#calsd=newsd*cal

#meancal=(meanspec-min(meanspec))*(wn/20.)**2/(hbb*resp*atm)
#sdcal=newsd*(wn/20.)**2/(hbb*resp*atm)
#relerr=sdcal/meancal

meancal = meanspec
sdcal = newsd

maxerr=5.
relerr = abs(sdcal/meancal)
#relerr = where(relerr < 0.,-relerr,relerr)
relerr = where(relerr > maxerr,maxerr,relerr)

plot.figure(7)
plot.clf()
print 'Spec:'
for d in range(nspec):
    #print len(detspec)
    plot.plot(wn,detspec[:,d],c='r',lw=0.5)
plot.plot(wn,detspec[:,0],c='r',lw=0.5,label='Individual detector')
plot.plot(wn,meanspec,c='k',lw=2.,label='Plotted mean')
plot.plot(wn,newmean,c='b',ls='--',lw=1.,label='Weighted mean')
plot.plot(wn,savspec,c='g',ls='--',lw=2.,label='IDL spectrum')
plot.legend(loc='upper right')
plot.draw()
plot.savefig('../RSRF/Locke/spec%s_mean.eps'%bandStr)
plot.savefig('../RSRF/Locke/spec%s_mean.png'%bandStr)

plot.figure(8)
plot.clf()
print 'SD:'
for d in range(nspec):
    #print len(detsd[:,d])
    plot.plot(wnsd,detsd[:,d],c='r',lw=0.5)
plot.plot(wn,detsd[:,d],c='r',lw=0.5,label='Individual detector')
plot.plot(wn,meansd,c='k',lw=2.,label='Plotted mean')
plot.plot(wn,newsd,c='b',ls='--',lw=2.,label='Simple SD')
#plot.plot(wnspec,quadsd,c='g',ls='--',lw=2.,label='Quadrature SD')

plot.legend(loc='upper left')
plot.draw()
plot.savefig('../RSRF/Locke/spec%s_sd.eps'%bandStr)
plot.savefig('../RSRF/Locke/spec%s_sd.png'%bandStr)


plot.figure(9)
plot.clf()
plot.plot(wn,meancal,c='k',ls='-',label='Measured spectrum')
#plot.plot(wn,calspec,c='k',ls='--',label='Calibrated spectrum')
plot.plot(wn,sdcal,c='k',ls='--',label='Calibrated SD')
plot.plot(wn,relerr,c='k',ls=':',label='Relative error')
plot.plot(wn,savrsrf,c='r',ls=':',label='IDL RSRF')
plot.plot(wn,savspec,c='r',ls=':',label='IDL spec')
#plot.legend(loc='upper left')
plot.yscale('log')
plot.draw()

plot.figure(10)
plot.clf()
plot.plot(wn,savrsrf,c='k',ls='-',label='RSRF')
plot.plot(wn,savrsrf*(1.-relerr),c='k',ls='--',label='Low')
plot.plot(wn,savrsrf*(1.+relerr),c='k',ls='--',label='High')
#plot.plot(wn,relerr,'r:',label='Relative Error')
plot.xlabel('Wavenumber')
plot.ylabel('RSRF')
plot.legend(loc='upper left')
plot.draw()

plot.figure(11)
plot.clf()
plot.plot(wn,meanspec,'k-')
plot.plot(wn,meanspecnofilt,'k:')
plot.plot(wn,meanspec*(1-relerr),'k--')
plot.plot(wn,meanspec*(1+relerr),'k--')
plot.yscale('log')
plot.draw()

plot.show()