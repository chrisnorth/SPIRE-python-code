
from numpy import zeros,zeros_like,sqrt,arctan2,pi,isnan,ones,arange,min,max,size
from numpy import argsort,sum,array

from scipy import where
from scipy.special import j1,j0
from scipy.signal import convolve2d
from scipy.interpolate import interp1d,interp2d
from scipy.integrate import trapz

from csv import reader
import string
import matplotlib.pyplot as plot
import sys

from beam import beam_azsm
from airy import airyconv,airyconv_MG

############################################
## Convolve airy ring with square pixel
############################################

##Units are lambda/D (l/D)
## Max radius is 3 l/D

#steps: 0.01 l/D

maxRad=10.
dd=0.01 ##pixel size in array

pixCtr = maxRad / dd

##make position/radius Arrays
npix = int(2.*maxRad/dd + 1)
xArr = zeros((npix,npix))
yArr = zeros_like(xArr)

for x in range(npix):
    xArr[x,:] = -maxRad + x*dd
    yArr[:,x] = -maxRad + x*dd

radArr = sqrt(xArr**2 + yArr**2)
azArr = arctan2(yArr,xArr)

##remove zero
radArr = where(radArr == 0.,1.e-12,radArr)

###radius is the angle in units of lambda/d

##Make Airy function
bessArr = j1(pi*radArr)
airyArr = ((2./pi) * (bessArr / radArr))**2
#sys.exit()

################################################
### Calculate aperture efficiency
################################################

radListTot=arange(0.,300.,0.01)
radListTot=where(radListTot == 0.,1.e-12,radListTot)
nRadTot=len(radListTot)
airyList= ((2./pi) * (j1(pi*radListTot)/radListTot))**2
apEffList=1. - (j0(pi*radListTot))**2 - (j1(pi*radListTot))**2.
airyTot=sum(airyList*2.*pi*radListTot*dd)

radList=arange(0.,maxRad,dd)
radList=where(radList == 0.,1.e-12,radList)
nRad=len(radList)

###Checking integration
apEffList2=zeros(nRad)
for r in range(nRad):
    apEffList2[r]=trapz(airyList[0:r]*2.*pi*radList[0:r]/airyTot,radList[0:r])

#sys.exit()
plot.figure(1)
plot.clf()
plot.plot(radListTot,airyList,'k-',label='Airy function')
plot.plot(radListTot,apEffList,'k--',label='Aperture Efficiency')
plot.xlim(0,3)
plot.xlabel('Angle')
plot.ylabel('Airy Fn / Ap Eff')
plot.savefig('../Plots/airyFn_apEff.png',transparent=False,dpi=300)
plot.savefig('../Plots/airyFn_apEff.eps',transparent=False)

################################################
### Read in MathCad values
################################################
fileMathcad='../Inputs/airy_conv_Mathcad.dat'
fMcad=reader(open(fileMathcad))
widList=[]
fwhmMathcad=[]
for row in fMcad:
    if string.find(row[0],'#') < 0:
        print row
        widList.append(float(row[2]))
        fwhmMathcad.append(float(row[3]))

widList=array(widList)
fwhmMathcad=array(fwhmMathcad)
print widList
print fwhmMathcad
##Define freq, wavelength and pixel width lists
#widList=array([0.3,0.4,0.5,0.6,0.7,1.,1.5,2.])
nuList=2.*widList
wlList=1./nuList

################################################
### Do convolution with kernel
################################################

### Test version
#widList=[1.0]

print widList
nWid=len(widList)

convArrCN=zeros((nRad,nWid))
apEffArrCN=zeros(nWid)
fwhmArrCN=zeros(nWid)
beamAreaCN=zeros(nWid)

convArrMG=zeros((3,nRad,nWid))
apEffArrMG=zeros((3,nWid))
fwhmArrMG=zeros((3,nWid))
beamAreaMG=zeros((3,nWid))

#convArr2=zeros((nRad,nWid))
#apEffArr2=zeros(nWid)
#fwhmArr2=zeros(nWid)
#beamArea2=zeros(nWid)

extent=min(xArr),max(xArr),min(yArr),max(yArr)

for w in range(nWid):
    wid=widList[w]

    ###############################################
    ## Full convolution
    ###############################################
    print 'Convolving with %.2f lambda/D pixel...'%(wid)

    convAzCN=airyconv(xArr,yArr,airyArr,wid,radList)

    apEffArrCN[w]=max(convAzCN)
    convAzCN=convAzCN/max(convAzCN)
    
    convArrCN[:,w]=convAzCN
    convIntCN=interp1d(radList,convAzCN,kind='cubic')
    convSortCN=argsort(convArrCN[:,w])
    convInvCN=interp1d(convArrCN[convSortCN,w],radList[convSortCN])
    fwhmArrCN[w]=2.*convInvCN(0.5) ## in units of l/D

    beamAreaCN[w]=trapz(convAzCN*2.*pi*radList*wlList[w]**2.,radList) ##in units of (l0/D)^2
    print 'FHWM (CN)= %.4f'%(fwhmArrCN[w])
    print 'ApEff (CN) =%.4f'%(apEffArrCN[w])
    print 'Beam(FWHM) (CN)=%.4f'%(convIntCN(fwhmArrCN[w]/2.))
    print 'Beam Area (CN)=%.4f'%(beamAreaCN[w])

    ################################################
    ##Matt's approximate method
    ###############################################
    print "Doing Matt's approximate method"
    (convMG1,convMG2,convAzMG)=airyconv_MG(wid,radList)    
    convMGall=array([convMG1,convMG2,convAzMG])

    for c in range(3):
        apEffArrMG[c,w]=max(convMGall[c,:])
        convMGall[c,:]=convMGall[c,:]/max(convMGall[c,:])
        convArrMG[c,:,w]=convMGall[c,:]
        convIntMG=interp1d(radList,convMGall[c,:])
        convSortMG=argsort(convArrMG[c,:,w])
        convInvMG=interp1d(convArrMG[c,convSortMG,w],radList[convSortMG])
        fwhmArrMG[c,w]=2.*convInvMG(0.5) ## in units of l/D

        beamAreaMG[c,w]=trapz(convMGall[c,:]*2.*pi*radList*wlList[w]**2.,radList) ##in units of (l0/D)^2
        print 'FHWM (MG %d)= %.4f'%(c,fwhmArrMG[c,w])
        print 'ApEff (MG %d)=%.4f'%(c,apEffArrMG[c,w])
        print 'Beam(FWHM) (MG %d)=%.4f'%(c,convIntMG(fwhmArrMG[c,w]/2.))
        print 'Beam Area (MG %d)=%.4f'%(c,beamAreaMG[c,w])
    
    ######################################
    ## In-line convolution for testing purposes only
    ######################################
#    npKer=wid / dd
#    print 'wid:',wid
#    print 'dd:',dd
#    kernel=ones((npKer,npKer))
#    kernel = kernel * (dd**2.) / airyTot
#    print 'airyTot:',airyTot
#    convSame=convolve2d(airyArr,kernel,mode="same")
#    
#    convAz2=beam_azsm(xArr,yArr,convSame,radList,nsamp=360)
#    
#    apEffArr2[w]=max(convAz2)
#    convArr2[:,w]=convAz2/max(convAz2)
#
#    convSort2=argsort(convArr2[:,w])
#    convInv2=interp1d(convArr2[convSort,w],radList[convSort2])
#    fwhmArr2[w]=2.*convInv2(0.5)
#
#    beamArea2[w]=trapz(convAz2*2.*pi*radList,radList)

#    print 'FHWM 2= %.4f'%(fwhmArr2[w])
#    print 'ApEff 2=%.4f'%(apEffArr2[w])
#
#    #print 'Beam(FWHM)=%.4f'%(convInt(fwhmArr[w]/2.))
#    
#    print 'Beam Area 2=%.4f'%(beamArea2[w])

    #########################################
    ## Plotting for testing putposes only
    #########################################

#    kernelMask = zeros_like(airyArr)
#    kernelMask=where(abs(xArr) <= wid/2.,1.,0.)
#    kernelMask=where(abs(yArr) <= wid/2.,kernelMask,0.)    
#    apEffArr2[w]=sum(airyArr*kernelMask)* (dd**2.) / airyTot
    
#    plot.figure(2)
#    plot.clf()
#    plot.subplot(2,2,1)
#    plot.imshow(airyArr,extent=extent,vmin=0,vmax=1)
#    plot.title('Airy function')
#    
#    plot.subplot(2,2,2)
#    plot.imshow(kernelMask,extent=extent,vmin=0,vmax=1)
#    plot.title('Kernel mask')
#    
#    plot.subplot(2,2,3)
#    plot.imshow(convSame,extent=extent,vmin=0,vmax=1)
#    plot.title('Convolution')
#    plot.colorbar()
#    
#    plot.figure(3)
#    plot.clf()
#    plot.plot(xArr,airyArr[:,pixCtr],'k-')
#    plot.plot(xArr,convSame[:,pixCtr],'k--')


###########################################
##Set up output files
###########################################
fileOut='../Inputs/airy_conv_CN_MG.dat'
fOut=open(fileOut,'w')
line='#nu_norm, wl_norm, pix width, FWHM (CN), FWHM (MG), Ap Eff (CN), ApEff (MG), Beam Area (CN), Beam Area (MG)\n'
fOut.write(line)
for w in range(nWid):
    line='%.3f, %.3f, %.3f, %.5f, %.5f, %.5f, %.5f, %.5f, %.5f\n'% \
        (nuList[w],wlList[w],widList[w],\
        fwhmArrCN[w],fwhmArrMG[2,w],apEffArrCN[w],apEffArrMG[2,w],\
        beamAreaCN[w],beamAreaMG[2,w])
    fOut.write(line)
fOut.close()

fileBeam='../Inputs/airy_beams_MG.dat'
fBeam=open(fileBeam,'w')
line='#Rad'
for w in range(nWid):
    line='%s,%.3f'%(line,widList[w])
line=line+'\n'
fBeam.write(line)

for r in range(nRad):
    line='%.3f'%(radList[r])
    for w in range(nWid):
        line='%s,%.5g'%(line,convArrCN[r,w])
    line=line+'\n'
    fBeam.write(line)
fBeam.close()

#######################################
## Plot results of convolution
#######################################

plot.figure(4)
plot.clf()
styles=('-','--',':','-.')
colors=('k','r','0.9','g','0.8','b','0.7','c','0.6','m','0.5','y','0.3')
lCN,=plot.plot(radList,convArrCN[:,0],color=colors[0],linestyle='-',label=r'CN (conv)'%(widList[0]))
lMG2,=plot.plot(radList,convArrMG[2,:,0],color=colors[0],linestyle='--',label=r'MG (avg))'%(widList[0]))
#lMG0,=plot.plot(radList,convArrMG[0,:,0],color=colors[w],linestyle='-.',label=r'MG (axis)'%(widList[0]))
#lMG1,=plot.plot(radList,convArrMG[1,:,0],color=colors[w],linestyle=':',label=r'MG (diag)'%(widList[0]))
for w in range(nWid):
    plot.plot(radList,convArrCN[:,w],color=colors[w],linestyle='-',label=r'%.2f $\lambda/D$'%(widList[w]))
    plot.plot(radList,convArrMG[2,:,w],color=colors[w],linestyle='-.',label='_nolegend_')
    #plot.plot(radList,convArrMG[0,:,w],color=colors[w],linestyle='--',label='_nolegend_')
    #plot.plot(radList,convArrMG[1,:,w],color=colors[w],linestyle=':',label='_nolegend_')
plot.xlim(0,3)
plot.xlabel('Angular radius')
plot.ylabel('Beam response')
plot.legend(loc='upper right')
plot.savefig('../Plots/airy_conv_CN_MG.png',transparent=False,dpi=300)
plot.savefig('../Plots/airy_conv_CN_MG.eps',transparent=False)

plot.figure(5)
plot.clf()
styles=('-','--',':','-.')
colors=('k','r','0.9','g','0.8','b','0.7','c','0.6','m','0.5','y','0.3')
for w in range(nWid):
    plot.plot(radList,(convArrMG[2,:,w]-convArrCN[:,w])/convArrCN[:,w],color=colors[w],linestyle='-',label=r'%.2f $\lambda/D$'%(widList[w]))
    #plot.plot(radList,convArrMG[0,:,w],color=colors[w],linestyle='--',label='_nolegend_')
    #plot.plot(radList,convArrMG[1,:,w],color=colors[w],linestyle=':',label='_nolegend_')
plot.xlim(0,3)
plot.xlabel('Angular radius')
plot.ylabel('Relative Beam response (MG-CN)/CN')
plot.savefig('../Plots/airy_conv_rel_CN_MG.png',transparent=False,dpi=300)
plot.savefig('../Plots/airy_conv_rel_CN_MG.eps',transparent=False)

plot.figure(6)
plot.clf()
plot.plot(radListTot*2.,apEffList,'k-',label='Square side')
plot.plot(radListTot*sqrt(2),apEffList,'k--',label='Diagonal')
plot.plot(widList,apEffArrCN,'ko',label='Square pixel (CN)')
plot.plot(widList,apEffArrMG[2,:],'ks',label='Square pixel (MG)')
plot.xlim(0,3)
plot.xlabel('Pixel Width')
plot.ylabel('Aperture Efficiency')
plot.legend(loc='lower right')
plot.savefig('../Plots/ApEff_CN_MG.png',transparent=False,dpi=300)
plot.savefig('../Plots/ApEff_CN_MG.eps',transparent=False)

plot.figure(7)
plot.clf()
plot.plot(widList,(apEffArrMG[2,:]-apEffArrCN)/apEffArrCN,'k-')
plot.xlim(0,1.5)
plot.xlabel('Pixel Width')
plot.ylabel('Relative Aperture Efficiency (MG-CN)/CN')
plot.savefig('../Plots/ApEff_rel_CN_MG.png',transparent=False,dpi=300)
plot.savefig('../Plots/ApEff_rel_CN_MG.eps',transparent=False)

plot.figure(8)
plot.clf()
plot.plot(widList,fwhmArrCN,'k-',label='FWHM (CN)')
plot.plot(widList,fwhmArrMG[2,:],'k--',label='FWHM (MG)')
plot.plot(widList,fwhmMathcad,'ko',label='FWHM (MathCad)')
plot.xlabel('Pixel width')
plot.ylabel('FWHM')
plot.legend(loc='upper left')
plot.savefig('../Plots/FWHM_CN_MG_MathCad.png',transparent=False,dpi=300)
plot.savefig('../Plots/FWHM_CN_MG_MathCad.eps',transparent=False)

plot.figure(9)
plot.clf()
plot.plot(widList,(fwhmArrMG[2,:]-fwhmArrCN)/fwhmArrCN,'k-',label='FWHM (MG-CN)/CN')
plot.plot(widList,(fwhmMathcad-fwhmArrCN)/fwhmArrCN,'k--',label='FWHM (MathCad-CN)/CN')
plot.xlabel('Pixel width')
plot.ylabel('Relative FWHM (MG-CN)/CN')
plot.legend(loc='upper left')
plot.savefig('../Plots/FWHM_rel_CN_MG_MathCad.png',transparent=False,dpi=300)
plot.savefig('../Plots/FWHM_rel_CN_MG_MathCad.eps',transparent=False)

plot.figure(10)
plot.clf()
plot.plot(widList,beamAreaCN,'k-',label='Beam Area (CN)')
plot.plot(widList,beamAreaMG[2,:],'k--',label='Beam Area (MG)')
plot.xlabel('Pixel width')
plot.ylabel('Beam Area')
plot.legend(loc='upper right')
plot.savefig('../Plots/beamarea_CN_MG.png',transparent=False,dpi=300)
plot.savefig('../Plots/beamarea_CN_MG.eps',transparent=False)

plot.figure(11)
plot.clf()
plot.plot(widList,(beamAreaMG[2,:]-beamAreaCN)/beamAreaCN,'k-')
plot.xlabel('Pixel width')
plot.ylabel('Relative Beam Area (MG-CN)/CN')
plot.savefig('../Plots/beamarea_rel_CN_MG.png',transparent=False,dpi=300)
plot.savefig('../Plots/beamarea_rel_CN_MG.eps',transparent=False)

##############################################
#####   Do interpolations
##############################################
###
#####Normalised aperture efficiency
###apEffIntCN=interp1d(nuList,apEffArrCN,kind='cubic')
###apEffCN0=apEffIntCN(1.0)
###apEffNormCN=interp1d(nuList,apEffArrCN/apEffCN0,kind='cubic')
###
#####Normalised FWHM
###fwhmIntCN=interp1d(nuList,fwhmArrCN,kind='cubic')
###fwhmCN0=fwhmIntCN(1.0)
###fhwmNormCN=interp1d(nuList,fwhmArrCN/fwhmCN0,kind='cubic')
###
#####Normalised beam area
###areaIntCN=interp1d(nuList,beamAreaCN,kind='cubic')
###areaCN0=areaIntCN(1.0)
###areaNormCN=interp1d(nuList,beamAreaCN/areaCN0,kind='cubic')
###
#####normalised wavelength list
###wlNorm=interp1d(nuList,wlList,kind='cubic')
###
#####Normalised frequency list
###nuNorm=arange(0.5,1.51,0.01)
###
###########################################
##### Set up bandpass variables
###########################################
#####Spectral resolution
###Res=3
###
###nu0=1.0
###nuL=1. - 1./(2*Res) ##lower freq
###nuU=1. + 1./(2*Res) ##upper freq
###
###nuB=array([nuL,nu0,nuU])
###wlB=1./nuB
###widB=0.5*nuB
###
###print 'Convolving with %.3f lambda/D pixel (lower edge of band)'%(widB[0])
###beamL=airyconv(xArr,yArr,airyArr,widB[0],radList)
###beamL=beamL/max(beamL)
###fwhmL=fwhmInt(nuB[0])
###print 'FWHM=%.3f'%(fwhmL)
###beamLInt=interp1d(radList,beamL,kind='cubic')
###print 'Beam(FWHM)=%.3f'%(beamLInt(fwhmL))
###
###print 'Convolving with %.3f lambda/D pixel (centre of band)'%(widB[1])
###beamC=airyconv(xArr,yArr,airyArr,widB[1],radList)
###beamC=beamC/max(beamC)
###fwhmC=fwhmInt(nuB[1])
###print 'FWHM=%.3f'%(fwhmC)
###beamCInt=interp1d(radList,beamC,kind='cubic')
###print 'Beam(FWHM)=%.3f'%(beamCInt(fwhmC))
###
###print 'Convolving with %.3f lambda/D pixel (upper edge of band)'%(widB[2])
###beamU=airyconv(xArr,yArr,airyArr,widB[2],radList)
###beamU=beamU/max(beamU)
###fwhmU=fwhmInt(nuB[2])
###print 'FWHM=%.3f'%(fwhmU)
###beamUInt=interp1d(radList,beamU,kind='cubic')
###print 'Beam(FWHM)=%.3f'%(beamUInt(fwhmU))
###
###plot.figure(8)
###plot.clf()
###plot.plot(radList,beamC,'k-',label='Centre')
###plot.plot(radList/nuL,beamL,'k--',label='Lower edge')
###plot.plot(radList/nuU,beamU,'k:',label='Upper edge')
###plot.xlabel(r'Angle ($\lambda/D$)')
###plot.ylabel('Beam Response')
###plot.legend(loc='upper right')
###
######Plot beam profile at upper/lower edge
###plot.figure(9)
###plot.clf()
###for n in range(nWid):
###    plot.plot(radList,beamInt(radList,nuList[n]),linestyle='-')
###    plot.plot(radList,convArr[:,n],marker='o')
###plot.xlim(0,3)
###plot.xlabel('Radius')
###plot.ylabel('Beam profile')
###
###plot.figure(10)
###plot.clf()
###plot.plot(nuNorm,areaNorm(nuNorm),apEffNorm(nuNorm))

plot.show()