import pyfits
from scipy import where
from numpy import sum,isnan,isfinite,zeros,array,floor,min,max,mean,std,nan
#def maincode():
    
band=range(3)

filesin=['../Maps/tmc_all_psw_deg0_image_align500.fits',
         '../Maps/tmc_all_pmw_deg0_image_align500.fits',
         '../Maps/tmc_all_plw_deg0_image.fits']
filesmaskout=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_mask.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_mask.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_mask.fits']
filesmeanout=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_mean.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_mean.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_mean.fits']
filesstdout=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_stdev.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_stdev.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_stdev.fits']
filesrngout=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_range.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_range.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_range.fits']
filessigout=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_sigma.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_sigma.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_sigma.fits']
filessig2out=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_sigma2.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_sigma2.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_sigma2.fits']
filesrmnout=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_relmean.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_relmean.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_relmean.fits']
filesout=['../Maps/tmc_all_psw_deg0_image_align500_smoothed_all.fits',
              '../Maps/tmc_all_pmw_deg0_image_align500_smoothed_all.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed_all.fits']

for b in band:
    print 'Band: %d'%(b)
    fitsfile=filesin[b]
    mapin=pyfits.getdata(fitsfile,0)
    hdrin=pyfits.getheader(fitsfile,0)
    nxin=hdrin.get('NAXIS2')
    nyin=hdrin.get('NAXIS1')
    cpx=hdrin.get('CRPIX2')
    cpy=hdrin.get('CRPIX1')
    cpxval=hdrin.get('CRVAL2')
    cpyval=hdrin.get('CRVAL1')
    dx=hdrin.get('CDELT2')*3600. #deg->arcsec
    dy=hdrin.get('CDELT1')*3600. #deg->arcsec
    print '  Map size: %d x %d pixels'%(nxin,nyin)
    print '  Map resolution: %d x %d arcsec'%(abs(dx),abs(dy))
    
    #mask_full=where(isnan(mapin),0.,1.)
    #pixmask=where(mask_full == 1)
    #npm=pixmask.size
    print '  Number of pixels: %d'%(mapin.size)
    #print '  Number of non-blank pixels: %d'%(sum(mask_full))

    lores=280.
    loker=array((abs(lores/dx),abs(lores/dy)))
    lop=int(abs(lores/dx)/2)
    #lomap=zeros((floor(2*nxin/lop),floor(2*nyin/lop)))
    
    mapsm=zeros((nxin,nyin,7)) #mask,mean,sd,range
    px=lop
    py=lop
    while px < nxin:
        print 'row %d of %d [%.2f%%]'%(px,nxin,100.*px/nxin)
        while py < nyin:
            subarr=mapin[px-lop:px+lop,py-lop:py+lop]
            #print '%d,%d: %.2f'%(px,py,sum(subarr))
            if isfinite(sum(subarr)):
                mapsm[px-lop:px+lop,py-lop:py+lop,0]=1.
                mapsm[px-lop:px+lop,py-lop:py+lop,1]=mean(subarr)
                mapsm[px-lop:px+lop,py-lop:py+lop,2]=std(subarr)
                mapsm[px-lop:px+lop,py-lop:py+lop,3]=max(subarr)-min(subarr)
                mapsm[px-lop:px+lop,py-lop:py+lop,4]=abs(mean(subarr)/std(subarr))
                mapsm[px-lop:px+lop,py-lop:py+lop,5]=(mean(subarr)/max(subarr)-min(subarr))
                mapsm[px-lop:px+lop,py-lop:py+lop,6]=subarr/mean(subarr)
            else:
                mapsm[px-lop:px+lop,py-lop:py+lop,:]=nan
            py = py + 2*lop
        px = px + 2*lop
        py=lop
    #for px in range(lop,nxin-lop):
    #    print 'row %d of %d [%.2f%%]'%(px,nxin-lop*2,100.*px/(nxin-lop*2))
    #    for py in range(lop,nyin-lop):
    #        subarr=mapin[px-lop:px+lop,py-lop:py+lop]
    #        if subarr.size != (2*lop)**2:
    #            print subarr.shape
    #            stop
    #        if isfinite(sum(subarr)):
    #            mapsm[px,py,0]=1.
    #            mapsm[px,py,1]=mean(subarr)
    #            mapsm[px,py,2]=std(subarr)
    #            mapsm[px,py,3]=max(subarr)-min(subarr)
    

    hdrori=hdrin.copy()
    hdrori.update('EXTNAME','InputMap',comment='Input Map')
    hdrin.update('SMTH_WID',str(lores),comment='Width of smoothed regions')
    
    hdrmsk=hdrin.copy()
    hdrmsk.update('EXTNAME','Mask',comment='Mask for valid regions')
    hdrmean=hdrin.copy()
    hdrmean.update('EXTNAME','Mean',comment='Mean of smoothed region')
    hdrstd=hdrin.copy()
    hdrstd.update('EXTNAME','Std. Dev',comment='Standard Deviation of smoothed region')
    hdrrng=hdrin.copy()
    hdrrng.update('EXTNAME','Range',comment='Range of values in smoothed region')
    hdrsig=hdrin.copy()
    hdrsig.update('EXTNAME','SigmaS',comment='(Mean / Std. dev) of values in smoothed region')
    hdrsig2=hdrin.copy()
    hdrsig2.update('EXTNAME','SigmaR',comment='(Mean / Range) of values in smoothed region')
    hdrrmn=hdrin.copy()
    hdrrmn.update('EXTNAME','RelMean',comment='Values relative to mean in smoothed region')
    
    prihdu=pyfits.PrimaryHDU()
    hduori = pyfits.ImageHDU(mapin[:,:],header=hdrori)
    hdumsk = pyfits.ImageHDU(mapsm[:,:,0],header=hdrmsk)
    hdumean = pyfits.ImageHDU(mapsm[:,:,1],header=hdrmean)
    hdustd = pyfits.ImageHDU(mapsm[:,:,2],header=hdrstd)
    hdurng = pyfits.ImageHDU(mapsm[:,:,3],header=hdrrng)
    hdusig = pyfits.ImageHDU(mapsm[:,:,4],header=hdrsig)
    hdusig2 = pyfits.ImageHDU(mapsm[:,:,5],header=hdrsig2)
    hdurmn = pyfits.ImageHDU(mapsm[:,:,6],header=hdrrmn)
    
    #hdulist=pyfits.HDUList([prihdu,hdumsk])
    #hdulist.writeto(filesmaskout[b],clobber=True)
    #hdulist=pyfits.HDUList([prihdu,hdumean])
    #hdulist.writeto(filesmeanout[b],clobber=True)
    #hdulist=pyfits.HDUList([prihdu,hdustd])
    #hdulist.writeto(filesstdout[b],clobber=True)
    #hdulist=pyfits.HDUList([prihdu,hdurng])
    #hdulist.writeto(filesrngout[b],clobber=True)
    #hdulist=pyfits.HDUList([prihdu,hdusig])
    #hdulist.writeto(filessigout[b],clobber=True)
    #hdulist=pyfits.HDUList([prihdu,hdusig2])
    #hdulist.writeto(filessig2out[b],clobber=True)
    #hdulist=pyfits.HDUList([prihdu,hdurmn])
    #hdulist.writeto(filesrmnout[b],clobber=True)
    
    hdulist=pyfits.HDUList([prihdu,hduori,hdumsk,hdumean,hdustd,hdurng,hdusig,hdusig2,hdurmn])
    hdulist.writeto(filesout[b],clobber=True)
    
#if __name__ == "__main__":
#    maincode()