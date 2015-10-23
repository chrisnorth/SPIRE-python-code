bands = ['PLW', 'PMW', 'PSW']

def gen_thumb(obsid, outdir):
    obsin = getObservation(obsid, useHsa=True)
    try:
        for band in bands:
            nominal = obsin.level2.getProduct("extd"+band)
            simpleFitsWriter(nominal, outdir +  '%d_NOMINAL_%s.fits' % (obsid, band))
    except:
        pass
    
#obs_ids = asciiTableReader(file='/Users/Tom/projects/spire-hires-testing/obs-lists/ids-abovethresh-lowsnr.csv')[0].getData()
#for obsid in obs_ids:
#    gen_thumb(obsid, '/Users/Tom/HIPE/plots/obs_thumbs/above_lowsnr/')
#    
#obs_ids = asciiTableReader(file='/Users/Tom/projects/spire-hires-testing/obs-lists/ids-belowthresh-highsnr.csv')[0].getData()
#for obsid in obs_ids:
#    gen_thumb(obsid, '/Users/Tom/HIPE/plots/obs_thumbs/below_highsnr/')
#    
obs_ids = asciiTableReader(file='/Users/Tom/projects/spire-hires-testing/obs-lists/nearby.csv')[0].getData()
for obsid in obs_ids:
    gen_thumb(obsid, '/Users/Tom/HIPE/plots/obs_thumbs/nearby/')
    
