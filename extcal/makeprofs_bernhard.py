import pyfits
from numpy import arange,zeros_like,sqrt,floor,isfinite
from scipy import where
from beam import beam2prof,beam2prof2,beamarea_az,modbeam_area

###############################################################
### Read in Bernhard's Beam Maps
### Produce beam profiles
###############################################################

######################################
### Set filenames                  ###
######################################

fitsIn=['../Bernhard/0x5000241aL_PSW_pmcorr_1arcsec.fits', \
   '../Bernhard/0x5000241aL_PMW_pmcorr_1arcsec.fits', \
   '../Bernhard/0x5000241aL_PLW_pmcorr_1arcsec.fits']

fitsClean=['../Bernhard/0x5000241aL_PSW_pmcorr_1arcsec_cln_bgsub.fits', \
    '../Bernhard/0x5000241aL_PMW_pmcorr_1arcsec_cln_bgsub.fits', \
    '../Bernhard/0x5000241aL_PLW_pmcorr_1arcsec_cln_bgsub.fits']

fitsNorm=['../Bernhard/0x5000241aL_PSW_pmcorr_1arcsec_norm_beam.fits', \
    '../Bernhard/0x5000241aL_PMW_pmcorr_1arcsec_norm_beam.fits', \
    '../Bernhard/0x5000241aL_PLW_pmcorr_1arcsec_norm_beam.fits']

bands=['PSW','PMW','PLW']

##################################################
### Set background level and central positions ###
##################################################

bgLev=[-8.e-6,-2.e-5,-2.e-5]
peakNorm=[150.49599008,96.70399888,58.608041759]
cPxs=[[965,1119],[966,1139],[977,1138]]

########################################
### Read in maps & calculate radii   ###
########################################
    
mapsIn=[]
mapsClean=[]
mapsNorm=[]
masksIn=[]
masksClean=[]
xArrs=[]
yArrs=[]
radArrs=[]
radIntArrs=[]
for b in range(3):
    print 'Reading maps for band %d'%(b+1)
    mapsIn.append(pyfits.getdata(fitsIn[b]))
    mapsClean.append(pyfits.getdata(fitsClean[b]))
    mapsNorm.append(pyfits.getdata(fitsNorm[b]))

    masksIn.append(where(isfinite(mapsIn[b]),1.,0.))
    mapsIn[b]=where(isfinite(mapsIn[b]),mapsIn[b],0.)

    masksClean.append(where(isfinite(mapsClean[b]),1.,0.))
    mapsClean[b]=where(isfinite(mapsClean[b]),mapsClean[b],0.)
    
    beamx=zeros_like(mapsIn[b])
    beamy=zeros_like(mapsIn[b])

    hdr=pyfits.getheader(fitsIn[b],1)
    nx0=hdr.get('NAXIS1')
    ny0=hdr.get('NAXIS2')
    dx0=abs(hdr.get('CDELT1'))*3600. #in arcsec

    ##################################################
    ###   Calculate x,y position
    ##################################################
        
    for x in range(nx0):
        beamx[:,x]=(x - cPxs[b][0]) * dx0
        beamy[:,x]=(arange(ny0) - cPxs[b][1]) * dx0
    beamrad = sqrt(beamx**2 + beamy**2)
    beamradint=floor(beamrad)
    xArrs.append(beamx)
    yArrs.append(beamy)
    radArrs.append(beamrad)
    radIntArrs.append(beamradint)

#############################################
### Make beam profiles                    ###
#############################################

maxRad=1000.
radList=arange(maxRad)
nRad=len(radList)

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

outFiles=['../Inputs/beamProfs_Bernhard_cln_bgsub_PSW.dat', \
        '../Inputs/beamProfs_Bernhard_cln_bgsub_PMW.dat', \
        '../Inputs/beamProfs_Bernhard_cln_bgsub_PLW.dat']
for b in range(3):
    print 'Calculating beam profile for band %d (Input)'%(b+1)
    (prof,area,nPix,nVal)=beam2prof2(radArrs[b],mapsClean[b],radList,mapMask=masksClean[b],retAll=True,verbose=True)
    profsClean.append(prof)
    areasClean.append(area)
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
    line='#rad,In,Norm PSW,Area,nPix,nVal\n'
    oFile.write(line)
    for r in range(nRad):
        line='%.1f, %.9g, %.9g, %d, %d\n'%\
            (radList[r],profsClean[b][r],areasClean[b][r],nPixClean[b][r],nValClean[b][r])
        oFile.write(line)
    oFile.close()
