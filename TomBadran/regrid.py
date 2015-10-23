# Settings
outDir = "/Users/Tom/HIPE/plots/SPIRE/"

obsids = [1342249237, 1342227726, 1342210936, 1342216940]
bands = ['PLW'] #, 'PMW', 'PSW']

beamSizes = {
    'PLW' : [20],
    'PMW' : [15],
    'PSW' : [10],
}

for obsid in obsids:
    for band in bands:
        beamSize = beamSizes[band][0]
        nominal = simpleFitsReader(outDir + '%d_NOMINAL_%s.fits' % (obsid, band))
        hires = simpleFitsReader(outDir + str(obsid) + '_HIRES_' + band +'_BEAMHSIZE_' + str(beamSize) + '.fits')
        rnom = regrid(source=nominal, target=hires)
        rnom = imageMultiply(rnom, 4)
        simpleFitsWriter(rnom, outDir + '%d_NOMINAL_REGRID_%s.fits' % (obsid, band))

