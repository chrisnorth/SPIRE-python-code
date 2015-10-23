from pylab import *
from astropy.table import Table
ion()
data = Table.read('stats/Results-floats2.csv',format='ascii.csv')
types = Table.read('obs-lists/obs-types_orig.csv', format='ascii.csv')
#data = loadtxt('stats/Results-floats2.csv', delimiter=',', skiprows=1, unpack=True)
#types = loadtxt('obs-lists/obs-types.csv', delimiter=',', skiprows=1, unpack=True)

#extract columns
obsids=data['Observation ID']
plw_99=data['PLW 99S']
plw_pix=data['PLW Pix > 20']
pmw_99=data['PMW 99S']
pmw_pix=data['PMW Pix > 30']
psw_99=data['PSW 99S']
psw_pix=data['PSW Pix > 100']

#thresholds
plw_pix_th=100
pmw_pix_th=278
psw_pix_th=544
plw_pix_val=20
pmw_pix_val=30
psw_pix_val=100
plw_99_th=10
pmw_99_th=15
psw_99_th=30

#psw_extra = loadtxt('stats/results-PSW.csv', delimiter=',', skiprows=7, unpack=True)
#pmw_extra = loadtxt('stats/results-PMW.csv', delimiter=',', skiprows=7, unpack=True)

# Extra Galactic
# GAL = 1
# COS = 2

# Glactic
# STA = 3
# ISM = 4
# SOL = 5

# Everything else/unknown
# TOO = 6
# DDT = 7

# Obs Modes
# Small = Type 1
# Large = Type 2
# Parallel = Type 3
msize = 2

obsMode=types['obsMode']
propType=types['PropType']
obsidType=types['Obs. Ids']
propName=types['Proposal']

props={
    'VNGS':{'code':'KPGT_cwilso01_1','color':'k'},
    'KINGFISH':{'code':'KPOT_rkennicu_1','color':'r'},
    'HEROES':{'code':'GT1_mbaes_1','color':'b'},
    'CygnusA':{'code':'OT1_aedge_4','color':'g'},
    'NHEMESES':{'code':'OT1_bholwerd_1','color':'c'},
    #'HRS':{'code':'KPGT_seales01_1','color':'m'},
    'ETdust':{'code':'GT2_mbaes_2','color':'orange'},
    'HeViCS':{'code':'KPOT_jdavie01_1','color':'gray'}
}

for p in props:
    props[p]['ids']=obsidType[propName==props[p]['code']]
    props[p]['n']=len(props[p]['ids'])
    props[p]['filter']=np.in1d(obsids, props[p]['ids'])
    print('%s: %d'%(p,props[p]['n']))

exgal_test = logical_or(propType == 'GAL', propType == 'COS')
gal_test = logical_or(propType == 'STA', propType == 'ISM', propType == 'SOL')
other_test = logical_not(logical_or(exgal_test, gal_test))

gal_ids = obsidType[gal_test]
exgal_ids = obsidType[exgal_test]
other_ids = obsidType[other_test]

gal_filter = np.in1d(obsids, gal_ids)
exgal_filter = np.in1d(obsids, exgal_ids)
other_filter = np.in1d(obsids, other_ids)

print('Galactic: %d ; Extragalactic: %d; Other: %d (Total: %d)'%(len(gal_ids),len(exgal_ids),len(other_ids),len(gal_ids)+len(exgal_ids)+len(other_ids)))

small_ids = obsidType[obsMode == 'SpirePhotoSmallScan']
large_ids = obsidType[obsMode == 'SpirePhotoLargeScan']
parallel_ids = obsidType[obsMode == 'SpirePacsParallel']

small_filter = np.in1d(obsids, small_ids)
large_filter = np.in1d(obsids, large_ids)
parallel_filter = np.in1d(obsids, parallel_ids)

print('Large: %d ; Small: %d; Parallel: %d (Total: %d)'%(len(large_ids),len(small_ids),len(parallel_ids),len(large_ids)+len(small_ids)+len(parallel_ids)))

#small_filter_pmw = np.in1d(obsids, small_ids)
#large_filter_pmw = np.in1d(obsids, large_ids)
#parallel_filter_pmw = np.in1d(obsids, parallel_ids)

#small_filter_psw = np.in1d(obsids, small_ids)
#large_filter_psw = np.in1d(obsids, large_ids)
#parallel_filter_psw = np.in1d(obsids, parallel_ids)


#gal_filter_psw = np.in1d(obsids, gal_ids)
#exgal_filter_psw = np.in1d(obsids, exgal_ids)
#other_filter_psw = np.in1d(obsids, other_ids)

#gal_filter_pmw = np.in1d(obsids, gal_ids)
#exgal_filter_pmw = np.in1d(obsids, exgal_ids)
#other_filter_pmw = np.in1d(obsids, other_ids)

def run_filters():
    #PLW
    filt_plw = logical_and(plw_pix.data > plw_pix_th, plw_99.data > plw_99_th)
    obs_ids_pass = obsids[filt_plw]
    obs_ids_fail = obsids[logical_not(filt_plw)]
    print('PLW: %d pass test' % len(obs_ids_pass))
    print('PLW: %d fail test' % len(obs_ids_fail))
    savetxt('obs-to-hires/plw-pass-2.csv', obs_ids_pass, fmt='%d')

    filt_pmw = logical_and(pmw_pix.data > pmw_pix_th, pmw_99.data > pmw_99_th)
    obs_ids_pass = obsids[filt_pmw]
    obs_ids_fail = obsids[logical_not(filt_pmw)]
    print('PMW: %d pass test' % len(obs_ids_pass))
    print('PMW: %d fail test' % len(obs_ids_fail))
    savetxt('obs-to-hires/pmw-pass-2.csv', obs_ids_pass, fmt='%d')

    filt_psw = logical_and(psw_pix.data > psw_pix_th, psw_99.data > psw_99_th)
    obs_ids_pass = obsids[filt_psw]
    obs_ids_fail = obsids[logical_not(filt_psw)]
    print('PSW: %d pass test' % len(obs_ids_pass))
    print('PSW: %d fail test' % len(obs_ids_fail))
    savetxt('obs-to-hires/psw-pass-2.csv', obs_ids_pass, fmt='%d')

    ############################################################################

    print('Halving the area')
    filt = logical_and(plw_pix > plw_pix_th/2, plw_99 > plw_99_th)
    obs_ids_pass = obsids[filt]
    obs_ids_fail = obsids[logical_not(filt)]
    print('PLW: %d pass test' % len(obs_ids_pass))
    print('PLW: %d fail test' % len(obs_ids_fail))

    filt = logical_and(pmw_pix > pmw_pix_th/2, pmw_99 > pmw_99_th)
    obs_ids_pass = obsids[filt]
    obs_ids_fail = obsids[logical_not(filt)]
    print('PMW: %d pass test' % len(obs_ids_pass))
    print('PMW: %d fail test' % len(obs_ids_fail))

    filt = logical_and(psw_pix > psw_pix_th, psw_99 > psw_99_th)
    obs_ids_pass = obsids[filt]
    obs_ids_fail = obsids[logical_not(filt)]
    print('PSW: %d pass test' % len(obs_ids_pass))
    print('PSW: %d fail test' % len(obs_ids_fail))

    ############################################################################

    # Pass all
    pass_plw = logical_and(plw_pix.data > plw_pix_th, plw_99.data > plw_99_th)
    obs_ids_pass_plw = obsids[pass_plw]
    pass_pmw = logical_and(pmw_pix.data > psw_pix_th, pmw_99.data > pmw_99_th)
    obs_ids_pass_pmw = obsids[pass_pmw]
    pass_psw = logical_and(psw_pix.data > psw_pix_th, psw_99.data > psw_99_th)
    obs_ids_pass_psw = obsids[pass_psw]

    obs_pass_pmw_psw = obs_ids_pass_pmw[in1d(obs_ids_pass_pmw, obs_ids_pass_psw)]
    obs_pass_plw_psw = obs_ids_pass_plw[in1d(obs_ids_pass_plw, obs_ids_pass_psw)]
    obs_pass_plw_pmw = obs_ids_pass_plw[in1d(obs_ids_pass_plw, obs_ids_pass_pmw)]

    obs_pass_all = obs_ids_pass_plw[in1d(obs_ids_pass_plw, obs_pass_pmw_psw)]
    print('Pass all bands:', len(obs_pass_all))

    # Pass 2
    obs_pass_2_1 = obs_pass_pmw_psw[in1d(obs_pass_pmw_psw, obs_pass_all, invert=True)]
    obs_pass_2_2 = obs_pass_plw_psw[in1d(obs_pass_plw_psw, obs_pass_all, invert=True)]
    obs_pass_2_3 = obs_pass_plw_pmw[in1d(obs_pass_plw_pmw, obs_pass_all, invert=True)]
    obs_pass_2 = concatenate((obs_pass_2_1, obs_pass_2_2, obs_pass_2_3))
    print('Pass 2 Bands:', len(obs_pass_2))

    # Pass 1
    obs_pass_1_1 = obs_ids_pass_psw[in1d(obs_ids_pass_psw, obs_pass_all, invert=True)]
    obs_pass_1_2 = obs_ids_pass_pmw[in1d(obs_ids_pass_pmw, obs_pass_all, invert=True)]
    obs_pass_1_3 = obs_ids_pass_plw[in1d(obs_ids_pass_plw, obs_pass_all, invert=True)]
    obs_pass_1 = concatenate((obs_pass_1_1, obs_pass_1_2, obs_pass_1_3))
    obs_pass_1 = obs_pass_1[in1d(obs_pass_1, obs_pass_2, invert=True)]
    print('Pass 1 Bands:', len(obs_pass_1))

    # > threshold buts less than 5 SNR:
    filt = logical_and(plw_99.data < 5, plw_pix.data > plw_pix_th)
    obs_ids = obsids[filt]
    savetxt('obs-lists/ids-abovethresh-lowsnr.csv', obs_ids, fmt='%d')

    # < threshold but above 17 SNR
    filt = logical_and(plw_99.data > 17, plw_pix.data < plw_pix_th)
    obs_ids = obsids[filt]
    savetxt('obs-lists/ids-belowthresh-highsnr.csv', obs_ids, fmt='%d')

    # Nearby Thresholds
    filt = logical_and(plw_99.data > 5, plw_99.data < 15)
    filt2 = logical_and(plw_pix.data > plw_pix_th/2, plw_pix.data < plw_pix_th*2)
    filt = logical_and(filt, filt2)
    obs_ids = obsids[filt]
    #passes_pixcount = zeroes_like(obs_ids)

    savetxt('obs-lists/nearby.csv', obs_ids, fmt='%d')

    return(filt_psw,filt_pmw,filt_plw)

def do_plots(filt_psw=None,filt_pmw=None,filt_plw=None):
    rcParams.update({
        'axes.titlesize': 'x-large',
        'axes.labelsize': 'large',
        'axes.linewidth':2,
        'legend.fontsize':'large',
        'xtick.labelsize':'large',
        'xtick.major.width':2,
        'xtick.minor.width':1,
        'ytick.labelsize':'large',
        'ytick.major.width':2,
        'ytick.minor.width':1,
        'legend.markerscale':2,
        'lines.markeredgewidth':1
    })

    msize=3
    mstyle ='o'
    ma=0.5 #marker alpha
    colors={
        'exgal':(0,0.7,0,ma),'gal':(0,0,1,ma),'other':(1,0,0,ma),
        'large':(0,0.7,0,ma),'small':(0,0,1,ma),'parallel':(1,0,0,ma),
        'passed_out':(0.5,0.5,0.5,0),'passed_in':'w',
        'mec':'None',
        'threshold':(0.5,0.5,0.5)}

    figure(1,figsize=(10,12))
    clf()

    axhline(psw_pix_th, color=colors['threshold'],lw=2)
    axvline(psw_99_th, color=colors['threshold'],lw=2)

    loglog(psw_99[exgal_filter], psw_pix[exgal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='ExGal', c=colors['exgal'])
    loglog(psw_99[gal_filter], psw_pix[gal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Galactic', c=colors['gal'])
    loglog(psw_99[other_filter], psw_pix[other_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Other', c=colors['other'])

    for p in props:
        loglog(psw_99[props[p]['filter']], psw_pix[props[p]['filter']], 'o', markersize=msize*3, markeredgewidth=2,label=p, mec=props[p]['color'],mfc='None')

    xlabel('99th Percentile Signal')
    ylabel('Pixel Count > 100 MJy/sr')
    title('PSW')
    xlim(1,1e5)
    ylim(1,1e7)

    legend(loc='lower right', frameon=False, numpoints=1)
    savefig('doc/stats-both-thresholds-2-psw.pdf')

    ############################################################################

    figure(2,figsize=(16,7))
    clf()
    subplot(131)

    axhline(plw_pix_th, color=colors['threshold'], lw=2, ls='--')
    axvline(plw_99_th, color=colors['threshold'], lw=2, ls='--')

    #loglog(plw_99[filt_plw], plw_pix[filt_plw], 'o', markersize=msize, label='Passed', mec=colors['passed_out'],mfc=colors['passed_in'])
    #loglog(plw_99[filt_plw], plw_pix[filt_plw], 'o', markersize=msize*2, mec=colors['passed_out'],mfc=colors['passed_in'])

    loglog(plw_99[exgal_filter], plw_pix[exgal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='ExGal', c=colors['exgal'])
    loglog(plw_99[gal_filter], plw_pix[gal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Galactic', c=colors['gal'])
    loglog(plw_99[other_filter], plw_pix[other_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Other', c=colors['other'])

    xlabel('99th Percentile Signal')
    ylabel('Pixel Count > %d MJy/sr'%(plw_pix_val))
    title('PLW')
    xlim(1,1e4)
    ylim(1,1e7)
    annotate('%d pixels'%(plw_pix_th),[xlim()[1]/1.2,plw_pix_th*1.2],ha='right',color=colors['threshold'],fontweight='bold')
    annotate('%d MJy/sr'%(plw_99_th),[plw_99_th*1.2,ylim()[1]/1.2],va='top',color=colors['threshold'],rotation=90,fontweight='bold')
    legend(loc='lower right', frameon=False, numpoints=1)

    subplot(132)

    axhline(pmw_pix_th, color=colors['threshold'], lw=2, ls='--')
    axvline(pmw_99_th, color=colors['threshold'], lw=2, ls='--')
    #loglog(pmw_99[filt_pmw], pmw_pix[filt_pmw], 'o', markersize=msize, label='Passed', mec=colors['passed_out'],mfc=colors['passed_in'])
    #loglog(pmw_99[filt_pmw], pmw_pix[filt_pmw], 'o', markersize=msize*2, mec=colors['passed_out'],mfc=colors['passed_in'])

    loglog(pmw_99[exgal_filter], pmw_pix[exgal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='ExGal', c=colors['exgal'])
    loglog(pmw_99[gal_filter], pmw_pix[gal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Galactic', c=colors['gal'])
    loglog(pmw_99[other_filter], pmw_pix[other_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Other', c=colors['other'])

    xlabel('99th Percentile Signal')
    ylabel('Pixel Count > %d MJy/sr'%(pmw_pix_val))
    title('PMW')
    xlim(1,1e5)
    ylim(1,1e7)
    annotate('%d pixels'%(pmw_pix_th),[xlim()[1]/1.2,pmw_pix_th*1.2],ha='right',color=colors['threshold'],fontweight='bold')
    annotate('%d MJy/sr'%(pmw_99_th),[pmw_99_th*1.2,ylim()[1]/1.2],va='top',color=colors['threshold'],rotation=90,fontweight='bold')
    legend(loc='lower right', frameon=False, numpoints=1)

    subplot(133)

    axhline(psw_pix_th, color=colors['threshold'], lw=2, ls='--')
    axvline(psw_99_th, color=colors['threshold'], lw=2, ls='--')

    #loglog(psw_99[filt_psw], psw_pix[filt_psw], 'o', markersize=msize, label='Passed', mec=colors['passed_out'],mfc=colors['passed_in'])
    #loglog(psw_99[filt_psw], psw_pix[filt_psw], 'o', markersize=msize*2, mec=colors['passed_out'],mfc=colors['passed_in'])

    loglog(psw_99[exgal_filter], psw_pix[exgal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='ExGal', c=colors['exgal'])
    loglog(psw_99[gal_filter], psw_pix[gal_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Galactic', c=colors['gal'])
    loglog(psw_99[other_filter], psw_pix[other_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Other', c=colors['other'])

    xlabel('99th Percentile Signal')
    ylabel('Pixel Count > %d MJy/sr'%(psw_pix_val))
    title('PSW')
    xlim(1,1e5)
    ylim(1,1e7)
    annotate('%d pixels'%(psw_pix_th),[xlim()[1]/1.2,psw_pix_th*1.2],ha='right',color=colors['threshold'],fontweight='bold')
    annotate('%d MJy/sr'%(psw_99_th),[psw_99_th*1.2,ylim()[1]/1.2],va='top',color=colors['threshold'],rotation=90,fontweight='bold')
    legend(loc='lower right', frameon=False, numpoints=1)

    tight_layout()
    savefig('doc/snr-thresholds-2.pdf')

    ############################################################################

    # Plot same but for different split
    figure(3,figsize=(16,7))
    clf()
    subplot(131)

    axhline(plw_pix_th, color=colors['threshold'], lw=2, ls='--')
    axvline(plw_99_th, color=colors['threshold'], lw=2, ls='--')

    #loglog(plw_99[filt_plw], plw_pix[filt_plw], 'o', markersize=msize, label='Passed', mec=colors['passed_out'],mfc=colors['passed_in'])
    #loglog(plw_99[filt_plw], plw_pix[filt_plw], 'o', markersize=msize*2, mec=colors['passed_out'],mfc=colors['passed_in'])

    loglog(plw_99[large_filter], plw_pix[large_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Large', c=colors['large'])
    loglog(plw_99[small_filter], plw_pix[small_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Small', c=colors['small'])
    loglog(plw_99[parallel_filter], plw_pix[parallel_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Parallel', c=colors['parallel'])

    xlabel('99th Percentile Signal')
    ylabel('Pixel Count > %d MJy/sr'%(plw_pix_val))
    title('PLW')

    xlim(1,1e4)
    ylim(1,1e7)
    annotate('%d pixels'%(plw_pix_th),[xlim()[1]/1.2,plw_pix_th*1.2],ha='right',color=colors['threshold'],fontweight='bold')
    annotate('%d MJy/sr'%(plw_99_th),[plw_99_th*1.2,ylim()[1]/1.2],va='top',color=colors['threshold'],rotation=90,fontweight='bold')
    legend(loc='lower right', frameon=False, numpoints=1)

    subplot(132)
    axhline(pmw_pix_th, color=colors['threshold'], lw=2, ls='--')
    axvline(pmw_99_th, color=colors['threshold'], lw=2, ls='--')

    #loglog(pmw_99[filt_pmw], pmw_pix[filt_pmw], 'o', markersize=msize, label='Passed', mec=colors['passed_out'],mfc=colors['passed_in'])
    #loglog(pmw_99[filt_pmw], pmw_pix[filt_pmw], 'o', markersize=msize*2, mec=colors['passed_out'],mfc=colors['passed_in'])

    loglog(pmw_99[large_filter], pmw_pix[large_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Large', c=colors['large'])
    loglog(pmw_99[small_filter], pmw_pix[small_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Small', c=colors['small'])
    loglog(pmw_99[parallel_filter], pmw_pix[parallel_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Parallel', c=colors['parallel'])

    xlabel('99th Percentile Signal [MJy/sr]')
    ylabel('Pixel Count > %d MJy/sr'%(pmw_pix_val))
    title('PMW')
    xlim(1,1e5)
    ylim(1,1e7)

    annotate('%d pixels'%(pmw_pix_th),[xlim()[1]/1.2,pmw_pix_th*1.2],ha='right',color=colors['threshold'],fontweight='bold')
    annotate('%d MJy/sr'%(pmw_99_th),[pmw_99_th*1.2,ylim()[1]/1.2],va='top',color=colors['threshold'],rotation=90,fontweight='bold')
    legend(loc='lower right', frameon=False, numpoints=1)

    subplot(133)
    axhline(psw_pix_th, color=colors['threshold'], lw=2, ls='--')
    axvline(psw_99_th, color=colors['threshold'], lw=2, ls='--')
    #loglog(psw_99[filt_psw], psw_pix[filt_psw], 'o', markersize=msize, label='Passed', mec=colors['passed_out'],mfc=colors['passed_in'])
    #loglog(psw_99[filt_psw], psw_pix[filt_psw], 'o', markersize=msize*2, mec=colors['passed_out'],mfc=colors['passed_in'])

    loglog(psw_99[large_filter], psw_pix[large_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Large', c=colors['large'])
    loglog(psw_99[small_filter], psw_pix[small_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Small', c=colors['small'])
    loglog(psw_99[parallel_filter], psw_pix[parallel_filter], mstyle, mec=colors['mec'] , markersize=msize, label='Parallel', c=colors['parallel'])

    xlabel('99th Percentile Signal')
    ylabel('Pixel Count > %d MJy/sr'%(psw_pix_val))
    title('PSW')
    xlim(1,1e5)
    ylim(1,1e7)
    annotate('%d pixels'%(psw_pix_th),[xlim()[1]/1.2,psw_pix_th*1.2],ha='right',color=colors['threshold'],fontweight='bold')
    annotate('%d MJy/sr'%(psw_99_th),[psw_99_th*1.2,ylim()[1]/1.2],va='top',color=colors['threshold'],rotation=90,fontweight='bold')
    legend(loc='lower right', frameon=False, numpoints=1)

    tight_layout()
    savefig('doc/snr-thresholds-byobstype-2.pdf')
    show()

if __name__ == '__main__':
    filt_psw,filt_pmw,filt_plw = run_filters()
    do_plots(filt_psw,filt_pmw,filt_plw)
