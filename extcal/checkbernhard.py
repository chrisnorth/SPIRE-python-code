import pyfits
from numpy import arange,zeros_like,sqrt,floor,isfinite,array
from scipy import where
from beam import beam2prof,beam2prof2,beamarea_az,modbeam_area
from csv import reader
import string
import matplotlib.pyplot as plot

###########################################################
#### Checks to do:
### Plot incremental radius with completeness
### Tolerancing on BG level
### Tolerancing on peak value
### Plot bg maps
###########################################################

###########################################################
### Set normalisation constants                         ###
###########################################################

bgLev=[-8.e-6,-2.e-5,-2.e-5]
peakNorm=[150.49599008,96.70399888,58.608041759]

areas0=[450,795,1665]
areas1=[462,825,1690]
areaLim=[[445,465],[785,835],[1650,1700]]
bandStr=['PSW','PMW','PLW']

###########################################################
### Read in files                                       ###
###########################################################

fitsIn=['../Inputs/beamProfs_Bernhard_cln_bgsub_PSW.dat', \
        '../Inputs/beamProfs_Bernhard_cln_bgsub_PMW.dat', \
        '../Inputs/beamProfs_Bernhard_cln_bgsub_PLW.dat']

#profsIn=[]
#profsInNorm=[]
profsSubBg=[]
#areasIn=[]
areasSubBg=[]
nPixIn=[]
nValIn=[]
radList=[]
for b in range(3):
    fIn=reader(open(fitsIn[b],'r'))
    radList=[]
    #profIn=[]
    #profInNorm=[]
    profSubBg=[]
    #areaIn=[]
    areaBgSub=[]
    nPix=[]
    nVal=[]
    for row in fIn:
        if string.find(row[0],'#') < 0:
            radList.append(float(row[0]))
            #profIn.append(float(row[1]))
            profSubBg.append(float(row[1]))
            #profInNorm.append(float(row[2]))
            #areaIn.append(float(row[2]))
            areaBgSub.append(float(row[2]))
            nPix.append(float(row[3]))
            nVal.append(float(row[4]))
    #profsIn.append(array(profIn))
    profsSubBg.append(array(profSubBg))
    #profInNorm=(array(profIn)/peakNorm[b])
    #profsInNorm.append(array(profInNorm))
    areasSubBg.append(array(areaBgSub))
    nPixIn.append(array(nPix))
    nValIn.append(array(nVal))

radList=array(radList)

###########################################
### Subtract background                 ###
### Calcuate areas                      ###
###########################################
#profsInNorm=[]
#areasNorm=[]
#areasSubBg=[]
#for b in range(3):
    #profsSubBg.append(profsInNorm[b]-bgLev[b])
    #areasNorm.append(modbeam_area(radList,profsInNorm[b],radList))
    #areasSubBg.append(modbeam_area(radList,profsSubBg[b],radList))
    
for b in range(3):
    #######################################################
    ## Plot beam profile
    #######################################################    
    plot.figure(10*b+1)
    plot.clf()
    plot.plot(radList,profsSubBg[b],c='k')
    plot.xlim(0,1000)
    plot.yscale('log')
    plot.ylim(1e-8,1)
    plot.xlabel('radius [arcsec]')
    plot.ylabel('Beam Profile (%s)'%(bandStr[b]))
    plot.title(bandStr[b])
    plot.savefig('../Plots/beamprof_Bernhard_cln_bgsub_%s.png'%bandStr[b])
    plot.savefig('../Plots/beamprof_Bernhard_cln_bgsub_%s.eps'%bandStr[b])

    #######################################################
    ## Plot beam profile (zoom)
    #######################################################    
    plot.figure(10*b+2)
    plot.clf()
    plot.axhline(y=0.,c='k',ls='-')
    #plot.axhline(y=bgLev[b],c='k',ls=':')
    plot.plot(radList,profsSubBg[b],c='k')
    plot.xlim(0,1000)
    plot.ylim(-1e-4,1e-4)    
    plot.xlabel('radius [arcsec]')
    plot.ylabel('Beam Profile (%s)'%(bandStr[b]))
    plot.title(bandStr[b])
    plot.savefig('../Plots/beamprof-zoom_Bernhard_cln_bgsub_%s.png'%bandStr[b])
    plot.savefig('../Plots/beamprof-zoom_Bernhard_cln_bgsub_%s.eps'%bandStr[b])
    
    #######################################################
    ## Plot beam area
    #######################################################    
    plot.figure(10*b+3)
    plot.clf()
    #plot.plot(radList,areasNorm[b],c='b')
    plot.plot(radList,areasSubBg[b],c='k')
    plot.xlim(0,1000)
    plot.xlabel('radius [arcsec]')
    plot.ylabel('Beam Area (%s)'%(bandStr[b]))
    plot.title(bandStr[b])
    plot.savefig('../Plots/beamarea_Bernhard_cln_bgsub_%s.png'%bandStr[b])
    plot.savefig('../Plots/beamarea_Bernhard_cln_bgsub_%s.eps'%bandStr[b])
    
    #######################################################
    ## Plot coverage
    #######################################################    
    plot.figure(10*b+4)
    plot.clf()
    plot.plot(radList,nValIn[b]/nPixIn[b])
    plot.xlim(0,1000)
    plot.xlabel('radius [arcsec]')
    plot.ylabel('Map completeness (%s)'%(bandStr[b]))
    plot.title(bandStr[b])
    plot.savefig('../Plots/coverage_Bernhard_cln_bgsub_%s.png'%bandStr[b])
    plot.savefig('../Plots/coverage_Bernhard_cln_bgsub_%s.eps'%bandStr[b])
    
    #######################################################
    ## Plot beam area (zoom) & coverage
    #######################################################    
    plot.figure(10*b+5)
    plot.clf()
    plot.subplot(2,1,1)
    plot.axhline(y=areas0[b],c='k',ls=':')
    plot.axhline(y=areas1[b],c='k',ls=':')
    #plot.plot(radList,areasIn[b])
    #plot.plot(radList,areasNorm[b])
    plot.plot(radList,areasSubBg[b])
    plot.xlim(0,1000)
    plot.ylim(areaLim[b])
    #plot.xlabel('radius [arcsec]')
    plot.ylabel('Beam Area (%s)'%(bandStr[b]))
    plot.title(bandStr[b])
    
    plot.subplot(2,1,2)
    plot.plot(radList,nValIn[b]/nPixIn[b])
    plot.xlim(0,1000)
    plot.ylim(0,1)
    plot.yscale('linear')
    plot.xlabel('radius [arcsec]')
    plot.ylabel('Map completeness (%s)'%(bandStr[b]))
    #plot.title(bandStr[b])
    plot.savefig('../Plots/beamarea-zoom_cov_Bernhard_cln_bgsub_%s.png'%bandStr[b])
    plot.savefig('../Plots/beamarea-zoom_cov_Bernhard_cln_bgsub_%s.eps'%bandStr[b])
    
plot.show()