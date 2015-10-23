poolHsa = 'SPIRE_POOL'
poolPath = '/Users/Tom/HIPE/pool'

# Settings
outDir = "/Users/Tom/HIPE/plots/rotations/"

bands = ['PLW']
obsids = [1342249237]

beamSizes = {
    'PLW' : 80,
    'PMW' : 60,
    'PSW' : 40
}

angles = [0, 30, 60, 90, 120, 150, 180]

def loadBeams(obsIn):
    fullBeams = {}
    for band in bands:
        fullBeams[band] = obsIn.calibration.getPhot().getProduct('BeamProfList').getProduct(band, 'fine')
    return fullBeams

# Get the observation from the archive if not already stored in the pool locally
def loadObservation(obsid):
    try:
       print 'Trying to get observation from pool ...'
       obs = getObservation(obsid=obsid, poolName=poolHsa, instrument='SPIRE', poolLocation=poolPath)
       obs.refs.remove('level0')
       obs.refs.remove('level0_5')
    except:
        print 'Failed to get observation from pool, downloading from HSA ...'
        obs = getObservation(obsid=obsid, useHsa=True, instrument='SPIRE')
        obs.refs.remove('level0')
        obs.refs.remove('level0_5')
        spiaSaveObs(obs, Pool=poolHsa, nameTag=str(obsid), PoolPath=poolPath)
    return obs

def processHiRes(inLevel1, inArray, inBeam, inWcs, inMapMin, inMapMax, fluxOffsets):
    level1 = Level1Context()
    selectRefs = filter(lambda x:x.meta['bbid'].value>>16==0xa103, inLevel1.refs)
    for ref in selectRefs:
        level1.addProduct(ref.getProduct().copy())

    wcs = inWcs.copy()

    #Update wcs cdelt, naxis and crpix values
    keyMap = {'cdelt':[float, 0.5], 'naxis':[int, 2.], 'crpix':[float, 2.]}
    map(lambda x: wcs.meta[x].setValue(keyMap[x[:-1]][0](wcs.meta[x].value*keyMap[x[:-1]][1])) , \
        map(lambda x:keyMap.keys()[x/2]+str(x%2+1), range(2*len(keyMap.keys()))))
    maps = []
    for angle in angles:
        rotBeam = rotate(image=inBeam, angle=angle)
        beamSize = beamSizes[band]
        bcenter = rotBeam.image.dimensions
        bcenter[0] = bcenter[0] / 2
        bcenter[1] = bcenter[1] / 2
        beam = crop(rotBeam, int(bcenter[0] - beamSize) , int(bcenter[1] - beamSize), int(bcenter[0] + beamSize+1), int(bcenter[1] + beamSize+1))

        hiresImage, hiresBeam = hiresMapper(level1, array = inArray, beam=beam, wcs=wcs, fluxOffset=fluxOffsets)
        tempIndx = hiresImage.image.where((hiresImage.image>5*inMapMax).or(hiresImage.coverage<1e-10))
        hiresImage.image[tempIndx] = Double.NaN
        maps.append(hiresImage)
    return maps

for obsid in obsids:
    obs = loadObservation(obsid)
    fullBeams = loadBeams(obs)
    
    fullBeams = loadBeams(obs)

    for band in bands:
        level1 = obs.getLevel1()
        level2 = obs.getLevel2()

        chanRelGains = obs.calibration.getPhot().chanRelGain
        level1RelGains = Level1Context()
        level1RelGains.meta = level1.meta

        for i in range(level1.getCount()):
            psp = level1.getProduct(i)
            if psp.type=="PPT": psp.setType("PSP") #for old Level 1 contexts
            psp = applyRelativeGains(psp, chanRelGains)
            level1RelGains.addProduct(psp)

        diag = level2.getProduct('psrc%sdiag'%band)
        level1Corrected,mapZero,diagZero, p4,p5 = destriper(level1=level1RelGains, array=band, withMedianCorrected=True, startParameters=diag)

        wcs = level2.getProduct('extd%s'%band).wcs
        mapMax = MAX(level2.getProduct('extd%s'%band).image[\
            level2.getProduct('extd%s'%band).image.where(IS_FINITE)])
        mapMin = MIN(level2.getProduct('extd%s'%band).image[\
            level2.getProduct('extd%s'%band).image.where(IS_FINITE)])

        # Get flux offsets
        fluxOffsetsExtd = level2.getProduct('extd%s'%band).meta['zPointOffset'].value
        beamAreaPipSr = obs.calibration.getPhot().getProduct('ColorCorrBeam').meta['beamPipeline%sSr'%band.capitalize()].value
        fluxOffsetsPsrc = fluxOffsetsExtd * 1.e6 * beamAreaPipSr
        # Run hires with flux offsets
        tempMaps = processHiRes(level1Corrected, band, fullBeams[band], wcs, mapMin, mapMax, fluxOffsetsPsrc)
        
        rotMaps = []
        for map in tempMaps:
            # Remove flux offset
            mapHiresPsrcNew = imageSubtract( image1=map, scalar=fluxOffsetsPsrc)
            
            # Convert units to MJy/sr
            beamAreaPipArc = obs.calibration.getPhot().getProduct('ColorCorrBeam').meta['beamPipeline%sArc'%band.capitalize()].value
            mapHiresExtdNew = convertImageUnit(image=mapHiresPsrcNew,beamArea=beamAreaPipArc,newUnit='MJy/sr')
            k4P = obs.calibration.getPhot().getProduct('FluxConvList')[0].meta['k4P_%s'%band].value
            k4E = obs.calibration.getPhot().getProduct('FluxConvList')[0].meta['k4E_%s'%band].value
            mapHiresExtdNew = imageDivide(image1=mapHiresExtdNew,scalar=k4P)
            mapHiresExtdNew = imageMultiply(image1=mapHiresExtdNew,scalar=k4E)
            mapHiresExtdNew = imageAdd(image1=mapHiresExtdNew,scalar=fluxOffsetsExtd)
            rotMaps.append(mapHiresExtdNew)
        
        for i in range(len(angles)):
            diff = imageSubtract(image1=rotMaps[0], image2=rotMaps[i], ref=0)
            simpleFitsWriter(rotMaps[i], outDir + str(obsid) + '_HIRES_ROT_' + str(angles[i]) + '_' + band + '.fits')
            simpleFitsWriter(diff, outDir + str(obsid) + '_HIRES_DIFF0_ROT_' + str(angles[i]) + '_' + band + '.fits')
        
        