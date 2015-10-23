import datetime

# Settings
poolHsa = 'SPITZER_M74'
poolPath = '/Users/Tom/HIPE/pool'
outDirBase = "/Users/Tom/HIPE/plots/"
radialised=False

# M74
#spitzer24Src = '/Users/Tom/Projects/spire-hires-testing/spitzer-fits/m74.fits'
#obsid = 1342189427

# NGC 4151
obsid = 1342188588
spitzer24Src = '/Users/Tom/Projects/spire-hires-testing/spitzer-fits/ngc4151.fits'

# M81
#obsid = 1342185538
#spitzer24Src = '/Users/Tom/Projects/spire-hires-testing/spitzer-fits/m81.fits'


outDir = outDirBase + str(obsid) + '/'
bands = ['PLW', 'PMW', 'PSW']

# Get the observation from the archive if not already stored in the pool locally
def loadObservation(obsid):
    try:
        print 'Trying to get observation from pool ...'
        obsIn = getObservation(obsid=obsid, poolName=poolHsa, instrument='SPIRE', poolLocation=poolPath)
        obsIn.refs.remove('level0')
        obsIn.refs.remove('level0_5')
    except:
        print 'Failed to get observation from pool, downloading from HSA ...'
        obsIn = getObservation(obsid=obsid, useHsa=True, instrument='SPIRE')
        obsIn.refs.remove('level0')
        obsIn.refs.remove('level0_5')
        spiaSaveObs(obsIn, Pool=poolHsa, nameTag=str(obsid), PoolPath=poolPath)
    return obsIn

# Load in the spitzer image
def loadSpitzer():
    spitzer = fitsReader(spitzer24Src)
    # Convert spitzer image from MJ/sr to J/pixel
    sim = spitzer.getImage()
    sim *= 1e6
    sqdeg2sr = 3.0462e-4
    pix_size = abs(spitzer.wcs.cdelt1 * spitzer.wcs.cdelt2) * sqdeg2sr
    sim *= pix_size
    return spitzer

def loadBeams():
    fullBeams = {}
    for band in bands:
        
        fullBeams[band] = obsIn.calibration.getPhot().getProduct('BeamProfList').getProduct(band, 'fine')
        
        if radialised:
            radial_profile = obsIn.refs["calibration"].product.refs["Phot"].product.refs["RadialCorrBeam"].product["core"][band].data
            width = 2 * len(radial_profile) - 1
            beam = Double2d(width, width, Double.NaN)
            cx = cy = len(radial_profile) - 1
            
            for x in range(0, width):
                for y in range(0, width):
                    r = SQRT((x - cx)**2 + (y - cy)**2)
                    if r <= len(radial_profile)-1:
                        r1 = int(FLOOR(r))
                        r2 = int(CEIL(r))
                        val1 = radial_profile[r1]
                        val2 = radial_profile[r2]
                        # Linearly interpolate the radial profile - might want to bicubic this
                        if r2 == r1:
                            beam[x, y] = val1
                        else:        
                            beam[x, y] = val1 + (val2 - val1) * (r - r1) / (r2 - r1)
            fullBeams[band].image = beam
            fullBeams[band].wcs.naxis1 = width
            fullBeams[band].wcs.naxis2 = width
            fullBeams[band].wcs.crpix1 = cx
            fullBeams[band].wcs.crpix2 = cy
                    
        if radialised:
            simpleFitsWriter(fullBeams[band], outDirBase + 'RADIAL_BEAM_' + band + '.fits')
        else:
            simpleFitsWriter(fullBeams[band], outDirBase + 'BEAM_' + band + '.fits')
    return fullBeams

truthImages = {}
spitzerConvolvedImages = {}
def generateTruthImages():
    # Generate 'Truth' images and fake observation sources from SPITZER source
    for band in bands:
        beam = fullBeams[band]
         
        # Convolve spitzer image with beam
        print 'Convolving spitzer source with beam ...'
        
        # Generate Convolution beams
        nwcs = beam.wcs.copy()
        nwcs.cdelt1 = spitzer.wcs.cdelt1
        nwcs.cdelt2 = spitzer.wcs.cdelt2
        nwcs.naxis1 = spitzer.wcs.naxis1
        nwcs.naxis2 = spitzer.wcs.naxis2
        nwcs.crpix1 = (nwcs.naxis1 + 1) / 2
        nwcs.crpix2 = (nwcs.naxis2 + 1) / 2
        
        # Generate 'truth' beam
        twcs = beam.wcs.copy()
        twcs.cdelt1 = spitzer.wcs.cdelt1 * 2
        twcs.cdelt2 = spitzer.wcs.cdelt2 * 2
        twcs.naxis1 = spitzer.wcs.naxis1
        twcs.naxis2 = spitzer.wcs.naxis2
        twcs.crpix1 = (twcs.naxis1 + 1) / 2
        twcs.crpix2 = (twcs.naxis2 + 1) / 2
        
        # Regrid and fix NaNs
        nbeam = regrid(beam, wcs=nwcs)
        im = nbeam.getImage()
        im[im.where(IS_NAN(im))] = 0
        
        # Normalise back to 1.0 at peak 
        im *= 1/max(im)
        
        # Regrid truth beam
        tbeam = regrid(beam, wcs=twcs)
        tim = tbeam.getImage()
        tim[tim.where(IS_NAN(tim))] = 0
        tim *= 1/max(tim)
        
        convolvedSpitzer = imageConvolution(image=spitzer, kernel=nbeam)
        truthSpitzer = imageConvolution(image=spitzer, kernel=tbeam)
        
        # Remove background level
        im = convolvedSpitzer.getImage()
        im -= MEDIAN(NAN_FILTER(im))
        tim = truthSpitzer.getImage()
        tim -= MEDIAN(NAN_FILTER(tim))
        tim[tim.where(IS_NAN(tim))] = 0
        
        spitzerConvolvedImages[band] = convolvedSpitzer
        truthImages[band] = truthSpitzer
    return spitzerConvolvedImages, truthImages
   
def generateSpitzerObservations():
    newScans = {}
    
    # Create new level 1 observation with data from spitzer
    for band in bands:
        scanlines = obsIn.level1
        newScan = Level1Context()
        newScan.meta = scanlines.meta
        
        # Process each scan in the existing observation and generate new one with spitzer data and BG Noise
        for prodID in range(scanlines.count):
            print 'Processing scanline %d' % prodID
            line = scanlines.getProduct(prodID)
            # Fetch bolometer names for band we are processing
            bolometers = [name for name in line.getChannelNames() if name.startswith(band)]
            lines = scanlines.getProduct(prodID).copy()
               
            for b in bolometers:
                # Generate data for each bolometer
                ra = line.getRa(b)
                dec = line.getDec(b)
                sig = line.getSignal(b)
               
                # Loop over each value and get data from beam convolved spitzer image
                for i in range(len(sig)):
                    try:
                        # Load signal from spitzer
                        intensity = spitzerConvolvedImages[band].getIntensityWorldCoordinates(ra[i], dec[i])
                    except java.lang.IndexOutOfBoundsException:
                        instensity = Double.NaN
           
                    # Remove any NaNs from source data
                    if IS_NAN(intensity):
                        lines.getSignal(b)[i] = 0
                    else:
                        lines.getSignal(b)[i] = intensity
                       
            # Add scanline to new product
            newScan.addProduct(lines)
        newScans[band] = newScan
    return newScans
        
def processHiRes(inLevel1, inArray, inBeam, inWcs, inMapMin, inMapMax):
    level1 = Level1Context()
    selectRefs = filter(lambda x:x.meta['bbid'].value>>16==0xa103, inLevel1.refs)
    for ref in selectRefs: 
        level1.addProduct(ref.getProduct().copy())

    wcs = inWcs.copy()

    #Update wcs cdelt, naxis and crpix values
    keyMap = {'cdelt':[float, 0.5], 'naxis':[int, 2.], 'crpix':[float, 2.]}
    map(lambda x: wcs.meta[x].setValue(keyMap[x[:-1]][0](wcs.meta[x].value*keyMap[x[:-1]][1])) , \
        map(lambda x:keyMap.keys()[x/2]+str(x%2+1), range(2*len(keyMap.keys()))))

    hiresImage, hiresBeam = hiresMapper(level1, array = inArray, beam=inBeam, wcs=wcs)
    tempIndx = hiresImage.image.where((hiresImage.image>2*inMapMax).or(hiresImage.coverage<1e-10))
    hiresImage.image[tempIndx] = Double.NaN
    return hiresImage

# Load target observation and SPITZER image
obsIn = loadObservation(obsid)
spitzer = loadSpitzer()

# Load the beam images
fullBeams = loadBeams()

# Generate Truth Images
print 'Generating Truth Images ...'
spitzerConvolvedImages, truthImages  = generateTruthImages()

# Create fake observaitons
print 'Creating Fake Observation Data ...'
newScans = generateSpitzerObservations()

def simulationBeamSizes():
    # Generate HIRES maps using different beams
    beamSizes = {
        'PLW' : [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125, 150, 200],
        'PMW' : [15, 20, 25, 30, 35, 40, 45, 50, 60, 70, 80, 90, 100, 125],
        'PSW' : [5, 10, 15, 20, 25, 30, 35, 40, 50, 60, 75, 100],
    }
    
    runtime = []
    #beams = {}
    
    print 'Generating HiRes Images ..'
    hiresSample = {}
    for band in bands:
        for size in beamSizes[band]:
            t1 = datetime.datetime.now()
            bcenter = fullBeams[band].image.dimensions
            bcenter[0] = bcenter[0] / 2
            bcenter[1] = bcenter[1] / 2
            beam = crop(fullBeams[band], int(bcenter[0] - size) , int(bcenter[1] - size), int(bcenter[0] + size+1), int(bcenter[1] + size+1))
            #beams[band + str(size)] = beam
            
            level2 = obsIn.level2
            wcs = level2.getProduct('extd%s'%band).wcs
            mapMax = MAX(level2.getProduct('psrc%s'%band).image[\
                level2.getProduct('psrc%s'%band).image.where(IS_FINITE)])
            mapMin = MIN(level2.getProduct('psrc%s'%band).image[\
                level2.getProduct('psrc%s'%band).image.where(IS_FINITE)])
            tempMap = processHiRes(newScans[band], band, beam, wcs, mapMin, mapMax)
            simpleFitsWriter(tempMap, outDir + 'HIRES_' + band +'_BEAMHSIZE_' + str(size) + '.fits')
            hiresSample[band] = tempMap
            t2 = datetime.datetime.now()
            runtime.append( (t2-t1).seconds )
            
def number_iterations():
    pass

simulationBeamSize()

# Save out Truth images
for band in bands:
    truthImages[band] = regrid(source=truthImages[band], target=hiresSample[band])
    simpleFitsWriter(truthImages[band], outDir + 'TRUTH_' + band + '.fits')
    
for band in bands:
    wcs = obsIn.level2.getProduct("psrc"+band).wcs.copy()
    newmap = naiveScanMapper(newScans[band], array=band, wcs=wcs)
    newmap = regrid(source=newmap, target=hiresSample[band])
    simpleFitsWriter(newmap, outDir + 'ORIGINAL_' + band + '.fits')
    
  
