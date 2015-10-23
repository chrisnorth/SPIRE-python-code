import datetime

# Settings
obsid = 1342249237
band = 'PLW'
beamSize = 100

obs = getObservation(obsid=obsid, useHsa=True, instrument='SPIRE')
beam = obs.calibration.getPhot().getProduct('BeamProfList').getProduct(band, 'fine')
bcenter = beam.image.dimensions
bcenter[0] = bcenter[0] / 2
bcenter[1] = bcenter[1] / 2
beam = crop(beam, int(bcenter[0] - beamSize) , int(bcenter[1] - beamSize), int(bcenter[0] + beamSize+1), int(bcenter[1] + beamSize+1))

iters = range(30)
print(iters)
wcs = obs.level2.getProduct('extd%s'%band).wcs.copy()

#Update wcs cdelt, naxis and crpix values
keyMap = {'cdelt':[float, 0.5], 'naxis':[int, 2.], 'crpix':[float, 2.]}
map(lambda x: wcs.meta[x].setValue(keyMap[x[:-1]][0](wcs.meta[x].value*keyMap[x[:-1]][1])) , \
    map(lambda x:keyMap.keys()[x/2]+str(x%2+1), range(2*len(keyMap.keys()))))
  
times = []
l1 = obs.level1.copy()
for iter in iters:
    t1 = datetime.datetime.now()
    map, _ = hiresMapper(l1, maxIter=iter+1, beam=beam, wcs=wcs)
    t2 = datetime.datetime.now()
    times.append( (t2-t1).seconds )
    
print times
