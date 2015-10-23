# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 13:59:12 2012

@author: chris
"""

from csv import reader
import string
from numpy import array,zeros,arange,min,max
from scipy.interpolate import interp1d

import matplotlib.pyplot as plot


nuNorm=[]
wlNorm=[]
widNorm=[]
fwhmList=[]
apEffList=[]
beamAreaList=[]

fileConv='../Inputs/airy_conv.dat' ##hard-coded
fConv=reader(open(fileConv,'r'))
for row in fConv:
    if string.find(row[0],'#') < 0:
        nuNorm.append(float(row[0]))
        wlNorm.append(float(row[1]))
        widNorm.append(float(row[2]))
        fwhmList.append(float(row[3]))
        apEffList.append(float(row[4]))
        beamAreaList.append(float(row[5]))

nuNorm=array(nuNorm)
wlNorm=array(wlNorm)
widNorm=array(widNorm)
fwhmList=array(fwhmList)
apEffList=array(apEffList)
beamAreaList=array(beamAreaList)

fwhmInt=interp1d(nuNorm,fwhmList,kind='cubic')
apEffInt=interp1d(nuNorm,apEffList,kind='cubic')
beamAreaInt=interp1d(nuNorm,beamAreaList,kind='cubic')

fwhm0=fwhmInt(1.)
apEff0=apEffInt(1.)
beamArea0=beamAreaInt(1.)
fwhmNorm=interp1d(nuNorm,fwhmList/fwhm0,kind='cubic')
apEffNorm=interp1d(nuNorm,apEffList/apEff0,kind='cubic')
beamAreaNorm=interp1d(nuNorm,beamAreaList/beamArea0,kind='cubic')

nuArr=arange(min(nuNorm),max(nuNorm),0.01)

plot.figure(20)
plot.clf()
plot.plot(nuNorm,fwhmList,'ko')
plot.plot(nuArr,fwhmInt(nuArr),'k-')
plot.title('FWHM')

plot.figure(21)
plot.clf()
plot.plot(nuNorm,apEffList,'ko')
plot.plot(nuArr,apEffInt(nuArr),'k-')
plot.title('Ap Eff')

plot.figure(22)
plot.clf()
plot.plot(nuNorm,beamAreaList,'ko')
plot.plot(nuArr,beamAreaInt(nuArr),'k-')
plot.title('Ap Eff')
plot.show()