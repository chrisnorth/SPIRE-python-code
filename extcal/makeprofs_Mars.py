import pyfits
import sys
from numpy import arange,zeros_like,sqrt,floor,isfinite,isnan,min,max
from scipy import where
from beam import beam2prof,beam2prof2,beamarea_az,modbeam_area

###############################################################
### Read in Bernhard's Beam Maps
### Produce beam profiles
###############################################################

######################################
### Set filenames                  ###
######################################

map=4
if map ==1:
    fitsIn=['../Inputs/Mars_5000532c/Mars_5000532c_PSW-cutTrim.fits', \
        '../Inputs/Mars_5000532c/Mars_5000532c_PMW-cutTrim.fits', \
        '../Inputs/Mars_5000532c/Mars_5000532c_PLW-cutTrim.fits']
    outFiles=['../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_PSW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_PMW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_PLW.dat']
elif map ==2:
    fitsIn=['../Inputs/Mars_5000532c/Mars_5000532c_PSW-cutTrim.fits', \
        '../Inputs/Mars_5000532c/Mars_5000532c_PMW-cutTrim.fits', \
        '../Inputs/Mars_5000532c/Mars_5000532c_PLW-cutTrim.fits']
    outFiles=['../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_PSW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_PMW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_PLW.dat']
elif map==3:
    fitsIn=['../Inputs/Mars_5000532c/1arcsec/Mars_mapPSW_5000532c_1arcsec-cutTrim.fits', \
        '../Inputs/Mars_5000532c/1arcsec/Mars_mapPMW_5000532c_1arcsec-cutTrim.fits', \
        '../Inputs/Mars_5000532c/1arcsec/Mars_mapPLW_5000532c_1arcsec-cutTrim.fits']
    outFiles=['../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_1arcsec_PSW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_1arcsec_PMW.dat', \
        '../Inputs/Mars_5000532c/beamProfs_Mars_5000532c_1arcsec_PLW.dat']
elif map==4:
    fitsIn=['../Inputs/Mars_50004e4b/1arcsec/Mars_mapPSW_50004e4b_1arcsec-cutTrim.fits', \
        '../Inputs/Mars_50004e4b/1arcsec/Mars_mapPMW_50004e4b_1arcsec-cutTrim.fits', \
        '../Inputs/Mars_50004e4b/1arcsec/Mars_mapPLW_50004e4b_1arcsec-cutTrim.fits']
    outFiles=['../Inputs/Mars_50004e4b/beamProfs_Mars_50004e4b_1arcsec_PSW.dat', \
        '../Inputs/Mars_50004e4b/beamProfs_Mars_50004e4b_1arcsec_PMW.dat', \
        '../Inputs/Mars_50004e4b/beamProfs_Mars_50004e4b_1arcsec_PLW.dat']

bands=['PSW','PMW','PLW']

##################################################
### Set background level and central positions ###
##################################################

########################################
### Read in maps & calculate radii   ###
########################################
    
mapsIn=[]
masksIn=[]
xArrs=[]
yArrs=[]
radArrs=[]
radIntArrs=[]
dx0=[]
cPxs=[[0,0],[0,0],[0,0]]

for b in range(3):
    print 'Reading maps for band %d'%(b+1)
    mapsIn.append(pyfits.getdata(fitsIn[b]))

    masksIn.append(where(isfinite(mapsIn[b]),1.,0.))
    masksIn[b]=where(isnan(mapsIn[b]),1.,masksIn[b])
    mapsIn[b]=where(isfinite(mapsIn[b]),mapsIn[b],0.)
 
    beamx=zeros_like(mapsIn[b])
    beamy=zeros_like(mapsIn[b])

    hdr=pyfits.getheader(fitsIn[b],0)
    nx0=hdr.get('NAXIS1')
    ny0=hdr.get('NAXIS2')
    dx0.append(abs(hdr.get('CDELT1'))*3600.) #in arcsec
    cPxs[b][0]=hdr.get('CRPIX1')
    cPxs[b][1]=hdr.get('CRPIX2')

    ##################################################
    ###   Calculate x,y positiona
    ##################################################
        
    for x in range(nx0):
        beamx[:,x]=(x - cPxs[b][0]) * dx0[b]
        beamy[:,x]=(arange(ny0) - cPxs[b][1]) * dx0[b]
    print cPxs[b]
    print min(beamx),max(beamx),min(beamy),max(beamy)
    #sys.exit()
    beamrad = sqrt(beamx**2 + beamy**2)
    beamradint=floor(beamrad)
    xArrs.append(beamx)
    yArrs.append(beamy)
    radArrs.append(beamrad)
    radIntArrs.append(beamradint)

#############################################
### Make beam profiles                    ###
#############################################

maxRad=700.

radList=[]
profsIn=[]
areasIn=[]
nPixIn=[]
nValIn=[]
weightsIn=[]
profsClean=[]
areasClean=[]
nPixClean=[]
nValClean=[]
#weightsClean=[]
#profsMasksIn=[]
#profsMasksClean=[]
profsInNorm=[]
#areasNorm=[]
areasClean=[]

for b in range(3):
    print 'Calculating beam profile for band %d (Input)'%(b+1)
    radList.append(arange(maxRad/dx0[b])*dx0[b])
    nRad=len(radList[b])
    (prof,area,nPix,nVal)=beam2prof2(radArrs[b],mapsIn[b],radList[b],mapMask=masksIn[b],retAll=True,verbose=True)
    profsClean.append(prof)
    areasClean.append(area*dx0[b]**2)
    nPixClean.append(nPix)
    nValClean.append(nVal)
    #weightsIn.append(nPix/nVal)
    #print 'Calculating beam profile for band %d (Clean)'%(b+1)
    #(prof,area,nPix,nVal)=beam2prof2(radArrs[b],mapsClean[b],radList,mapMask=masksClean[b],retAll=True)
    #profsClean.append(prof)
    #areasClean.append(area)
    #nPixClean.append(nPix)
    #nValClean.append(nVal)
    #weightsClean.append(nPix/nVal)
#    print 'Calculating beam profile for band %d (Input Mask)'%(b+1)
#    profsMasksIn.append(beam2prof2(radArrs[b],masksIn[b],radList,mapMask=None,retAll=False))
#    print 'Calculating beam profile for band %d (Clean Mask)'%(b+1)
#    profsMasksClean.append(beam2prof2(radArrs[b],masksClean[b],radList,mapMask=None,retAll=False))

    ############################################
    ### Normalise beam profiles              ###
    ############################################

    #profsInNorm.append((profsIn[b] - bgLev[b])/peakNorm[b])    

    #areasNorm.append(modbeam_area(radList,profsInNorm[b],radList))
    #areasClean.append(modbeam_area(radList,profsClean[b],radList))
    
    oFile=open(outFiles[b],'w')
    line='#rad,Norm,Area,nPix,nVal\n'
    oFile.write(line)
    for r in range(nRad):
        line='%.1f, %.9g, %.9g, %d, %d\n'%\
            (radList[b][r],profsClean[b][r],areasClean[b][r],nPixClean[b][r],nValClean[b][r])
        oFile.write(line)
    oFile.close()

import matplotlib.pyplot as plot
plot.figure(1)
plot.clf()
plot.plot(radList[0],profsClean[0],'b-')
plot.plot(radList[1],profsClean[1],'g-')
plot.plot(radList[2],profsClean[2],'r-')
plot.yscale('log')
plot.draw()
plot.show()
