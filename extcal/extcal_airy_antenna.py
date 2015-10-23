
from numpy import pi,arange,zeros,zeros_like,sqrt,arctan2,array
from scipy.special import j1
from scipy import where

import matplotlib.pyplot as plot

from airy import airyconv,getairyconv,calcAiryTot,makeAiryArr
from calc_conv import calc_k4uni,calc_kcuni,calc_k5uni,calc_kc2uni
from beam import getbeamareauni
from rsrf import getapeffuni

fileConvAbs='../Inputs/airy_conv_Final.dat'
#fileBeamAbs='../Inputs/airy_beams.dat'

################################
##Get input factors
################################
print 'Reading in inputs'
(fwhmIntAbs,apEffIntAbs,beamAreaIntAbs)=getairyconv(fileConvAbs,norm=True)
(fwhmAbs,apEffAbs,beamAreaAbs)=getairyconv(fileConvAbs,norm=False)

##Get antenna-coupled factors
fileApEffAnt='../Inputs/freq_ap-eff.csv'
apEffIntAnt=getapeffuni(fileApEffAnt)

fileBeamAreaAnt='../Inputs/freq_beam-area.csv'
beamAreaIntAnt=getbeamareauni(fileBeamAreaAnt)

####################################
##Calculate Airy functions
####################################
print 'Calculating Airy function'
maxRad=10.
dd=0.01 ##pixel size in array

(xArr,yArr,airyArr)=makeAiryArr(maxRad,dd)
airyTot=calcAiryTot()

radList=arange(0.,maxRad,dd)
radList=where(radList == 0.,1.e-12,radList)
nRad=len(radList)

########################################
## Set up bandpass variables
########################################
##Spectral resolution
Res=3.

nu0=1. #band centre
nuL=1. - 1./(2.*Res) #band lower edge
nuU=1. + 1./(2.*Res) #band upper edge

##lower,centre,upper
nuB=array([nuL,nu0,nuU])
wlB=1./nuB
widB=0.5*nuB

dNu=0.001

####print 'Convolving with %.3f lambda/D pixel (lower edge of band)'%(widB[0])
####beamL=airyconv(xArr,yArr,airyArr,widB[0],radList)
####beamL=beamL/max(beamL)
####
####print 'Convolving with %.3f lambda/D pixel (centre of band)'%(widB[1])
####beamC=airyconv(xArr,yArr,airyArr,widB[1],radList)
####beamC=beamL/max(beamC)
####
####print 'Convolving with %.3f lambda/D pixel (upper edge of band)'%(widB[2])
####beamU=airyconv(xArr,yArr,airyArr,widB[2],radList)
####beamU=beamU/max(beamU)
####beamLInt=beamIntAbs(radList,widB[0])
####beamCInt=beamIntAbs(radList,widB[1])
####beamUInt=beamIntAbs(radList,widB[2])
####
######################################################
###### Plot to check interpolation works
######################################################
####
####plot.figure(1)
####plot.clf()
####plot.plot(radList,beamC,'ko',label='Centre')
####plot.plot(radList,beamU,'ks',label='Upper')
####plot.plot(radList,beamL,'k^',label='Lower')
####
####plot.plot(radList,beamCInt,'k-',label='Centre (Int)')
####plot.plot(radList,beamUInt,'k-',label='Upper (Int)')
####plot.plot(radList,beamLInt,'k-',label='Lower (Int)')
####

print 'Calculating K factors'
##Set alpha array
alphaArr=arange(-4,5.5,0.5)
alpha0=-1.
nAlph=len(alphaArr)

resList=[3.,3.,3.,5.,10.,2.,3.]
nu0List=[0.,-3.,+3.,0.,0.,0.,0.]

nRes=len(resList)

k4arrAbs=zeros((nAlph,nRes))
kcarrAbs=zeros((nAlph,nRes))
k5arrAbs=zeros((nAlph,nRes))
kc2arrAbs=zeros((nAlph,nRes))

k4arrAnt=zeros((nAlph,nRes))
kcarrAnt=zeros((nAlph,nRes))
k5arrAnt=zeros((nAlph,nRes))
kc2arrAnt=zeros((nAlph,nRes))

for r in range(nRes):
    res=resList[r]
    nuL=1. - 1./(2.*res) #band lower edge
    nuU=1. + 1./(2.*res) #band upper edge
    nu0=1.0 + (nu0List[r]/100.)*(nuU-nuL)
    nuB=[nuL,nu0,nuU]

    print 'Res, nu0 = %d, %.2f'%(res,nu0)
    #####################################
    ##Set normalised frequency/wl/pixel width lists
    ######################################
    nuNorm=arange(nuB[0],nuB[2]+dNu,dNu)
    wlNorm=1./nuNorm
    widNorm=0.5*nuNorm

    fwhmNormAbs=fwhmIntAbs(nuNorm)
    apEffNormAbs=apEffIntAbs(nuNorm)
    beamAreaNormAbs=beamAreaIntAbs(nuNorm)

    print 'FWHM(nuL,nu0,nuU)=%.3f, %.3f, %.3f'%(fwhmAbs(nuL),fwhmAbs(nu0),fwhmAbs(nuU))
    print 'FWHM/nu(nuL,nu0,nuU)=%.3f, %.3f, %.3f'%(fwhmAbs(nuL)/nuL,fwhmAbs(nu0)/nu0,fwhmAbs(nuU)/nuU)
    print 'norm FWHM/nu(nuL,nu0,nuU)=%.3f, %.3f, %.3f'%(fwhmIntAbs(nuL)/nuL,fwhmIntAbs(nu0)/nu0,fwhmIntAbs(nuU)/nuU)
    print 'Beam Area (nuL,nu0,nuU)=%.3f, %.3f, %.3f'%(beamAreaAbs(nuL),beamAreaAbs(nu0),beamAreaAbs(nuU))
    print 'norm Beam Area (nuL,nu0,nuU)=%.3f, %.3f, %.3f'%(beamAreaIntAbs(nuL),beamAreaIntAbs(nu0),beamAreaIntAbs(nuU))
    apEffAnt=apEffIntAnt(nuNorm)
    beamAreaNormAnt=nuNorm**-1.75
    #beamAreaNormAnt=beamAreaIntAnt(nuNorm)

    ###################################################
    ## Calculate K params
    ###################################################
    for a in range(nAlph):    
        k4arrAbs[a,r]=calc_k4uni(alphaArr[a],nuNorm,apEffNormAbs,nu0)
        kcarrAbs[a,r]=calc_kcuni(alphaArr[a],nuNorm,apEffNormAbs,nu0,alpha0)
        k5arrAbs[a,r]=calc_k5uni(alphaArr[a],nuNorm,apEffNormAbs,beamAreaNormAbs,nu0)
        kc2arrAbs[a,r]=calc_kc2uni(alphaArr[a],nuNorm,apEffNormAbs,beamAreaNormAbs,nu0,alpha0)
        
        k4arrAnt[a,r]=calc_k4uni(alphaArr[a],nuNorm,apEffAnt,nu0)
        kcarrAnt[a,r]=calc_kcuni(alphaArr[a],nuNorm,apEffAnt,nu0,alpha0)
        k5arrAnt[a,r]=calc_k5uni(alphaArr[a],nuNorm,apEffAnt,beamAreaNormAnt,nu0)
        kc2arrAnt[a,r]=calc_kc2uni(alphaArr[a],nuNorm,apEffAnt,beamAreaNormAnt,nu0,alpha0)

#########################################
## Output to files
#########################################
filek4Abs='../Outputs/K4_Uniform_Absorber.csv'    
filekcAbs='../Outputs/Kc_Uniform_Absorber.csv'
filek5Abs='../Outputs/K5_Uniform_Absorber.csv'
filekc2Abs='../Outputs/Kc2_Uniform_Absorber.csv'

filek4Ant='../Outputs/K4_Uniform_Antenna.csv'
filekcAnt='../Outputs/Kc_Uniform_Antenna.csv'
filek5Ant='../Outputs/K5_Uniform_Antenna.csv'
filekc2Ant='../Outputs/Kc2_Uniform_Antenna.csv'

fk4Abs=open(filek4Abs,'w')
fkcAbs=open(filekcAbs,'w')
fk5Abs=open(filek5Abs,'w')
fkc2Abs=open(filekc2Abs,'w')
fk4Ant=open(filek4Ant,'w')
fkcAnt=open(filekcAnt,'w')
fk5Ant=open(filek5Ant,'w')
fkc2Ant=open(filekc2Ant,'w')

line1='#res'
line2='#alpha/nu0'
for r in range(nRes):
    line1='%s,%d'%(line1,resList[r])
    line2='%s,%d'%(line2,nu0List[r])
line1=line1+'\n'
line2=line2+'\n'

fk4Abs.write(line1)
fkcAbs.write(line1)
fk5Abs.write(line1)
fkc2Abs.write(line1)
fk4Ant.write(line1)
fkcAnt.write(line1)
fk5Ant.write(line1)
fkc2Ant.write(line1)

fk4Abs.write(line2)
fkcAbs.write(line2)
fk5Abs.write(line2)
fkc2Abs.write(line2)
fk4Ant.write(line2)
fkcAnt.write(line2)
fk5Ant.write(line2)
fkc2Ant.write(line2)

for a in range(nAlph):
    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,k4arrAbs[a,r])
    line=line+'\n'
    fk4Abs.write(line)
    
    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,kcarrAbs[a,r])
    line=line+'\n'
    fkcAbs.write(line)

    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,k5arrAbs[a,r])
    line=line+'\n'
    fk5Abs.write(line)

    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,kc2arrAbs[a,r])
    line=line+'\n'
    fkc2Abs.write(line)

    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,k4arrAnt[a,r])
    line=line+'\n'
    fk4Ant.write(line)
    
    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,kcarrAnt[a,r])
    line=line+'\n'
    fkcAnt.write(line)

    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,k5arrAnt[a,r])
    line=line+'\n'
    fk5Ant.write(line)

    line='%.1f'%(alphaArr[a])
    for r in range(nRes):
        line='%s,%.9f'%(line,kc2arrAnt[a,r])
    line=line+'\n'
    fkc2Ant.write(line)

fk4Abs.close()
fkcAbs.close()
fk5Abs.close()
fkc2Abs.close()
fk4Ant.close()
fkcAnt.close()
fk5Ant.close()
fkc2Ant.close()

##################################################
## Make new plots
##################################################

##Plot normalised FWHM
plot.figure(2)
plot.clf()
plot.plot(nuNorm,fwhmNormAbs*wlNorm,'k-')
plot.xlim(0.6,1.4)
plot.xlabel('Normalised frequency')
plot.ylabel('Normalised FHWM')

plot.figure(3)
plot.clf()
plot.plot(nuNorm,beamAreaNormAbs*apEffNormAbs,'k-')
plot.xlabel('Normalised frequency')
plot.ylabel('Beam Solid Angle')
plot.xlim(0.75,1.25)
plot.ylim(0.99,1.002)

##Plot normalised quantities
plot.figure(4)
plot.clf()
plot.plot(nuNorm,apEffNormAbs,'k--',label='Aperture Efficiency')
plot.plot(nuNorm,beamAreaNormAbs,'k-',label='Beam Area')
plot.plot(nuNorm,apEffNormAbs*beamAreaNormAbs,'k:',label='ApEff * Beam Area')
plot.xlabel('Normalised frequency')
plot.ylabel('Normalised quantity')
plot.title('Figure 4 - Beam Area & ApEff (Absorber)')
plot.xlim(0.7,1.3)
plot.ylim(0.5,1.5)
plot.legend(loc='lower center')

plot.figure(5)
plot.clf()
styles=('k-','k--','k:','k-.')
cols=(0,3,4,5)
for r in range(4):
    col=cols[r]
    plot.plot(alphaArr,kcarrAnt[:,col],styles[r],label='R=%d'%(resList[col]))
plot.legend(loc='lower right')
plot.xlabel('Spectral Index')
plot.ylabel('Colour Correction Factor')
plot.title('Figure 5 - KColP, var R (Antenna)')

plot.figure(6)
plot.clf()
styles=('k-','k--','k:')
cols=(0,1,2)
for r in range(3):
    col=cols[r]
    plot.plot(alphaArr,kcarrAnt[:,col],styles[r],label=r'$\nu0=%.2f$'%(nu0List[col]))
plot.legend(loc='lower right')
plot.xlabel('Spectral Index')
plot.ylabel('Colour Correction Factor')
plot.title('Figure 6 - KColP, var nu0 (Antenna)')

plot.figure(7)
plot.clf()
plot.plot(alphaArr,kc2arrAnt[:,0])
plot.xlabel('Spectral Index')
plot.ylabel('Colour Correction Factor')
plot.title('Figure 7 - KExt (Antenna)')

plot.figure(8)
plot.clf()
styles=('k-','k--','k:','k-.')
cols=(0,3,4,5)
for r in range(4):
    col=cols[r]
    plot.plot(alphaArr,kcarrAbs[:,col],styles[r],label='R=%d'%(resList[col]))
plot.legend(loc='lower right')
plot.xlabel('Spectral Index')
plot.ylabel('Colour Correction Factor')
plot.xlim(-4,4)
plot.ylim(0.85,1.05)
#plot.yticks((0.85,0.89,0.93,0.97,1.01,1.05))
plot.title('Figure 8 - KColP, var R (Absorber)')

plot.figure(9)
plot.clf()
styles=('k-','k--','k:')
cols=(0,1,2)
for r in range(3):
    col=cols[r]
    plot.plot(alphaArr,kcarrAbs[:,col],styles[r],label=r'$\nu0=%.2f$'%(nu0List[col]))
plot.legend(loc='lower right')
plot.xlabel('Spectral Index')
plot.ylabel('Colour Correction Factor')
plot.xlim(-4,4)
plot.ylim(0.8,1.05)
plot.title('Figure 9 - KColP, var nu0 (Absorber)')

plot.figure(10)
plot.clf()
plot.plot(alphaArr,kc2arrAbs[:,0],'k-')
plot.xlabel('Spectral Index')
plot.ylabel('Colour Correction Factor')
plot.xlim(-4,4)
plot.ylim(0.8,1.05)
plot.title('Figure 10 - KExt (Absorber)')

plot.show()