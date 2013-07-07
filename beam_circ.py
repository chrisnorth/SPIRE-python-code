#Test effect of circularising beam

from numpy import arange,zeros,sqrt,log,exp,sum,array,concatenate
import sys
import matplotlib.pyplot as plot
from scipy.interpolate import interp1d
from scipy import where
import pyfits
import aplpy

from beam import getbeam,getnewbeamprofile,comb_beam,beamarea,beamarea_az,beam_azsm

bandstr=['psw','pmw','plw']
beamtype='M'
bpix=1.
brad=500.
bwid=brad*2.
radArr=arange(0,brad+bpix,bpix)
bradInt=2.*brad
radArrInt=arange(0,bradInt+bpix,bpix)
nRad=len(radArr)

##Read in beam maps anf circularised profiles
beamMaps=[]
beamMapAreas=[]

beamProfs=[]
beamProfsInt=[]
beamProfAreas=[]
beamProfErr=[]

beamProfs2=[]
beamProfs2Int=[]
beamProf2Areas=[]
beamProf2Err=[]
beamHdr=[]
for b in range(3):
    (xArr,yArr,beamMapBand)=getbeam(beamtype,b,bpix,bwid,bzlim=None,verbose=True)
    rArr=sqrt(xArr**2 + yArr**2)
    beamMapBand[where(rArr > brad)] = 0.
    beamMaps.append(beamMapBand)    
    
    beamMapAreas.append(beamarea(xArr,yArr,beamMapBand,brad=brad))
    
    beamSclBand,beamFixBand=getnewbeamprofile(b,radArr,bzlim=None,verbose=True)
    beamCombBand=comb_beam(beamSclBand,beamFixBand)
    beamProfs.append(array(beamCombBand))
    beamCombBandInt=concatenate((beamCombBand,zeros(bradInt-brad)))
    beamProfsInt.append(beamCombBandInt)
    beamProfAreas.append(beamarea_az(radArr,beamCombBand,brad=brad))
    beamProfErr.append((beamProfAreas[b]-beamMapAreas[b])/beamMapAreas[b])
    
    beamProf2Band=beam_azsm(xArr,yArr,beamMapBand,radArr)
    
    beamProfs2.append(beamProf2Band)
    beamProf2BandInt=concatenate((beamProf2Band,zeros(bradInt-brad)))
    beamProfs2Int.append(beamProf2BandInt)
    beamProf2Areas.append(beamarea_az(radArr,beamProf2Band,brad=brad))
    beamProf2Err.append((beamProf2Areas[b]-beamMapAreas[b])/beamMapAreas[b])
        
print 'Beam areas (map):  [%.3f, %.3f, %.3f]'%(beamMapAreas[0],beamMapAreas[1],beamMapAreas[2])
print 'Beam areas (prof): [%.3f, %.3f, %.3f]'%(beamProfAreas[0],beamProfAreas[1],beamProfAreas[2])
print 'Beam areas (prof2): [%.3f, %.3f, %.3f]'%(beamProf2Areas[0],beamProf2Areas[1],beamProf2Areas[2])

print 'Relative error (prof-map): [%.3f, %.3f, %.3f]'%(100.*beamProfErr[0],100.*beamProfErr[1],100.*beamProfErr[2])
print 'Relative error (prof2-map): [%.3f, %.3f, %.3f]'%(100.*beamProf2Err[0],100.*beamProf2Err[1],100.*beamProf2Err[2])

#plot.figure(1)
#plot.clf()
#plot.plot(radArr,abs(beamProfs[0]-beamProfs2[0]),'b-')
#plot.plot(radArr,abs(beamProfs[1]-beamProfs2[1]),'g-')
#plot.plot(radArr,abs(beamProfs[2]-beamProfs2[2]),'r-')
##plot.plot(radArr,beamProfs[0],'g-')
##plot.plot(radArr,beamProfs[1],'g-')
##plot.plot(radArr,beamProfs[2],'r-')
#plot.plot(radArr,beamProfs2[0],'b--')
#plot.plot(radArr,beamProfs2[1],'g--')
#plot.plot(radArr,beamProfs2[2],'r--')
#plot.yscale('log')

#plot.show()

nXmap=len(xArr[:,0])
nYmap=len(yArr[0,:])
beamMapCirc=[]
beamDiffCirc=[]
beamDiffRel=[]
beamMapCirc2=[]
beamDiffCirc2=[]
beamDiffRel2=[]

fileCirc=[]
fileCirc2=[]
fileDiff=[]
fileDiff2=[]
fileDiffRel=[]
fileDiffRel2=[]
fileBeam=[]

for b in range(3):
    #beamMapCircBand=zeros((nXmap,nYmap))
    #beamMapCirc2Band=zeros((nXmap,nYmap))
    beamProfInterp=interp1d(radArrInt,beamProfsInt[b])
    beamProf2Interp=interp1d(radArrInt,beamProfs2Int[b])
    #for x in range(nXmap):
    #    beamMapCircBand[x,:]=interp1d()
    beamMapCirc.append(beamProfInterp(rArr))
    beamDiffCirc.append(beamMaps[b]-beamMapCirc[b])
    beamMapCirc2.append(beamProf2Interp(rArr))
    beamDiffCirc2.append(beamMaps[b]-beamMapCirc2[b])

    beamDiffRel.append(beamDiffCirc[b]/beamMaps[b])
    beamDiffRel2.append(beamDiffCirc2[b]/beamMaps[b])

    fileCirc.append('../Inputs/beamCirc_%s.fits'%(bandstr[b]))
    fileCirc2.append('../Inputs/beamCirc2_%s.fits'%(bandstr[b]))
    fileDiff.append('../Inputs/beamDiff_%s.fits'%(bandstr[b]))
    fileDiff2.append('../Inputs/beamDiff2_%s.fits'%(bandstr[b]))
    fileDiffRel.append('../Inputs/beamDiffRel_%s.fits'%(bandstr[b]))
    fileDiffRel2.append('../Inputs/beamDiffRel2_%s.fits'%(bandstr[b]))
    fileBeam.append('../Inputs/spire_beams_measured/%s_beam_1arcsec.fits'%(bandstr[b]))
    hdr=pyfits.getheader(fileBeam[b],0)
    hdr.update('CRPIX1',brad)
    hdr.update('CRPIX2',brad)
    pyfits.writeto(fileCirc[b],beamMapCirc[b],clobber=True,header=hdr)
    pyfits.writeto(fileCirc2[b],beamMapCirc2[b],clobber=True,header=hdr)
    pyfits.writeto(fileDiff[b],beamDiffCirc[b],clobber=True,header=hdr)
    pyfits.writeto(fileDiff2[b],beamDiffCirc2[b],clobber=True,header=hdr)
    pyfits.writeto(fileDiffRel[b],beamDiffRel[b],clobber=True,header=hdr)
    pyfits.writeto(fileDiffRel2[b],beamDiffRel2[b],clobber=True,header=hdr)

fig=aplpy.FITSFigure(fileBeam[0])
fig.show_grayscale(stretch='log',vmin=1.e-5,vmax=1.)
fig.recenter(0,0,2./60)
fig.ticks.set_yspacing(30./3600)
fig.tick_labels.set_yformat('DD:mm:ss')
fig.ticks.set_xspacing(30./3600)
fig.tick_labels.set_xformat('DD:mm:ss')

plot.show()

#fig.show_contour(fileDiffRel[0])


###define source map
#srcElip=2. #ratio of source FWHMs
#srcFWHM1 = 50. #FWHM of axis 1 in arcsec
#srcFWHM2 = srcFWHM1 * srcElip #FWHM of axis 2 in arcsec
#
#srcSigma1=srcFWHM1/sqrt(2.*log(2.))
#srcSigma2=srcFWHM2/sqrt(2.*log(2.))
#
#srcArr=exp(-xArr**2./(2.*srcSigma1**2.) - yArr**2./(2.*srcSigma2**2.))
#srcProf=beam_azsm(xArr,yArr,srcArr,radArr)
#srcArea=sum(srcArr)
#print 'srcArea: %.3f'%srcArea
#
#plot.figure(1)
#plot.clf()
#plot.plot(xArr[:,0],srcArr[:,brad],'r-')
#plot.plot(yArr[0,:],srcArr[brad,:],'b-')
#plot.plot(radArr,srcProf,'g-')
#plot.plot(-radArr,srcProf,'g-')
#plot.show()
