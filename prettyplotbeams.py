### Make pretty plots of the beams

def plot():
    import aplpy as ap

    fitsIn=['../CalProducts/0x5000241aL_PSW_pmcorr_1arcsec_cln_bgsub.fits', \
        '../CalProducts/0x5000241aL_PMW_pmcorr_1arcsec_cln_bgsub.fits', \
        '../CalProducts/0x5000241aL_PLW_pmcorr_1arcsec_cln_bgsub.fits', \
        '../CalProducts/0x5000241aL_PSW_pmcorr_6arcsec_cln_bgsub.fits', \
        '../CalProducts/0x5000241aL_PMW_pmcorr_10arcsec_cln_bgsub.fits', \
        '../CalProducts/0x5000241aL_PLW_pmcorr_14arcsec_cln_bgsub.fits']
    pngOut=['../CalProducts/psw_prettybeam_cln_bgsub_1arcsec.png', \
        '../CalProducts/pmw_prettybeam_cln_bgsub_1arcsec.png', \
        '../CalProducts/plw_prettybeam_cln_bgsub_1arcsec.png', \
        '../CalProducts/psw_prettybeam_cln_bgsub_6arcsec.png', \
        '../CalProducts/pmw_prettybeam_cln_bgsub_10arcsec.png', \
        '../CalProducts/plw_prettybeam_cln_bgsub_14arcsec.png']
    epsOut=['../CalProducts/psw_prettybeam_cln_bgsub_1arcsec.eps', \
        '../CalProducts/pmw_prettybeam_cln_bgsub_1arcsec.eps', \
        '../CalProducts/plw_prettybeam_cln_bgsub_1arcsec.eps', \
        '../CalProducts/psw_prettybeam_cln_bgsub_6arcsec.eps', \
        '../CalProducts/pmw_prettybeam_cln_bgsub_10arcsec.eps', \
        '../CalProducts/plw_prettybeam_cln_bgsub_14arcsec.eps']

    labels=[r'250$\,\mathbf{\mu}$m (1" pixels)', \
        r'350$\,\mathbf{\mu}$m (1" pixels)', \
        r'500$\,\mathbf{\mu}$m (1" pixels)', \
        r'250$\,\mathbf{\mu}$m (6" pixels)', \
        r'350$\,\mathbf{\mu}$m (10" pixels)', \
        r'500$\,\mathbf{\mu}$m (14" pixels)']
            

    RaCtr=326.009629253
    DecCtr=-14.0731951967

##old values   
#    fitsIn=['../Inputs/spire_beams_measured/psw_beam_1arcsec.fits', \
#        '../Inputs/spire_beams_measured/pmw_beam_1arcsec.fits', \
#        '../Inputs/spire_beams_measured/plw_beam_1arcsec.fits']
#
#    pngOut=['../Docs/Paper/Figs/psw_prettybeam.png', \
#        '../Docs/Paper/Figs/pmw_prettybeam.png', \
#        '../Docs/Paper/Figs/plw_prettybeam.png']
#    epsOut=['../Docs/Paper/Figs/psw_prettybeam.eps', \
#        '../Docs/Paper/Figs/pmw_prettybeam.eps', \
#        '../Docs/Paper/Figs/plw_prettybeam.eps']
#    RaCtr=0.
#    DecCtr=0.
    
    for m in range(len(fitsIn)):
        print 'Plotting map %d: %s'%(m+1,fitsIn[m])
        fig=ap.FITSFigure(fitsIn[m])
        fig.show_grayscale(stretch='log',vmin=1e-4,vmax=1,invert=True)
    
        fig.recenter(RaCtr,DecCtr,width=15/60.,height=15/60.)
        
        fig.hide_tick_labels()
        fig.hide_axis_labels()
    
        fig.add_label(0.5,0.95,labels[m],relative='True',weight='normal',size=24)
    
        fig.save(pngOut[m])
        fig.save(epsOut[m])
    
if __name__=="__main__":
    plot()