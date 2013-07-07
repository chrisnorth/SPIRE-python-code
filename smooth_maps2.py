import pyfits
from scipy import where
from numpy import sum,isnan,isfinite,zeros,array,floor,min,max,mean,std,nan,sqrt
#def maincode():
    
band=range(3)

lores=280.
file_base_psw='../Maps/tmc_all_psw_deg0_image_align500_smoothed2_%d'%(int(lores))
file_base_pmw='../Maps/tmc_all_pmw_deg0_image_align500_smoothed2_%d'%(int(lores))
file_base_plw='../Maps/tmc_all_plw_deg0_image_smoothed2_%d'%(int(lores))

filesin=['../Maps/tmc_all_psw_deg0_image_align500.fits',
         '../Maps/tmc_all_pmw_deg0_image_align500.fits',
         '../Maps/tmc_all_plw_deg0_image.fits']
filesmaskout=['../Maps/tmc_all_psw'+file_base+'_mask.fits',
              '../Maps/tmc_all_pmw'+file_base+'_mask.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_mask.fits']
filesmeanout=['../Maps/tmc_all_psw'+file_base+'_mean.fits',
              '../Maps/tmc_all_pmw'+file_base+'_mean.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_mean.fits']
filesstdout=['../Maps/tmc_all_psw'+file_base+'_stdev.fits',
              '../Maps/tmc_all_pmw'+file_base+'_stdev.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_stdev.fits']
filesrngout=['../Maps/tmc_all_psw'+file_base+'_range.fits',
              '../Maps/tmc_all_pmw'+file_base+'_range.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_range.fits']
filessigout=['../Maps/tmc_all_psw'+file_base+'_sigma.fits',
              '../Maps/tmc_all_pmw'+file_base+'_sigma.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_sigma.fits']
filessig2out=['../Maps/tmc_all_psw'+file_base+'_sigma2.fits',
              '../Maps/tmc_all_pmw'+file_base+'_sigma2.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_sigma2.fits']
filesrmnout=['../Maps/tmc_all_psw'+file_base+'_relmean.fits',
              '../Maps/tmc_all_pmw'+file_base+'_relmean.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_relmean.fits']
filesout=['../Maps/tmc_all_psw'+file_base+'_all.fits',
              '../Maps/tmc_all_pmw'+file_base+'_all.fits',
              '../Maps/tmc_all_plw_deg0_image_smoothed2_all.fits']

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

    loker=array((abs(lores/dx),abs(lores/dy)))
    lop=int(abs(lores/dx)/2)
    #lomap=zeros((floor(2*nxin/lop),floor(2*nyin/lop)))
    
    mapsm=zeros((nxin,nyin,7)) #mask,mean,sd,range
#    xarr=zeros((2*lop,2*lop))
#    yarr=zeros((2*lop,2*lop))
    rarr=zeros((2*lop,2*lop))
    for px in range(lop,nxin-lop):
        print 'row %d of %d [%.2f%%]'%(px,nxin-lop*2,100.*px/(nxin-lop*2))
        for py in range(lop,nyin-lop):
            subarr=mapin[px-lop:px+lop,py-lop:py+lop]
#            for x in range(2*lop):
#                xarr[x,:]=px-lop + x
#            for y in range(2*lop):
#                yarr[:,y]=py-lop + y
            #for x in range(2*lop):
            #    for y in range(2*lop):
            #        rarr[x,y]=sqrt((px-lop + x)**2 + (py-lop + y)**2)
            #inr=where(rarr < lop)
            if isfinite(sum(subarr)):
                #mapsm[px,py,0]=1.
                #mapsm[px,py,1]=mean(subarr[inr])
                #mapsm[px,py,2]=std(subarr[inr])
                #mapsm[px,py,3]=max(subarr[inr])-min(subarr[inr])
                #mapsm[px,py,4]=abs(mapsm[px,py,1]/mapsm[px,py,2])
                #mapsm[px,py,5]=mapsm[px,py,1]/mapsm[px,py,3]
                #mapsm[px,py,6]=mapin[px,py]/mapsm[px,py,1]
            
                mapsm[px,py,0]=1.
                mapsm[px,py,1]=mean(subarr)
                mapsm[px,py,2]=std(subarr)
                mapsm[px,py,3]=max(subarr)-min(subarr)
                mapsm[px,py,4]=abs(mapsm[px,py,1]/mapsm[px,py,2])
                mapsm[px,py,5]=mapsm[px,py,1]/mapsm[px,py,3]
                mapsm[px,py,6]=mapin[px,py]/mapsm[px,py,1]
            else:
                mapsm[px,py,:]=nan
    
    hdrori=hdrin.copy()
    hdrori.update('EXTNAME','InputMap',comment='Input Map')
    hdrin.update('SMTH_WID',str(lores),comment='Radius of smoothed regions')
    
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
    
    hdulist=pyfits.HDUList([prihdu,hdumsk])
    hdulist.writeto(filesmaskout[b],clobber=True)
    hdulist=pyfits.HDUList([prihdu,hdumean])
    hdulist.writeto(filesmeanout[b],clobber=True)
    hdulist=pyfits.HDUList([prihdu,hdustd])
    hdulist.writeto(filesstdout[b],clobber=True)
    hdulist=pyfits.HDUList([prihdu,hdurng])
    hdulist.writeto(filesrngout[b],clobber=True)
    hdulist=pyfits.HDUList([prihdu,hdusig])
    hdulist.writeto(filessigout[b],clobber=True)
    hdulist=pyfits.HDUList([prihdu,hdusig2])
    hdulist.writeto(filessig2out[b],clobber=True)
    hdulist=pyfits.HDUList([prihdu,hdurmn])
    hdulist.writeto(filesrmnout[b],clobber=True)
    
    hdulist=pyfits.HDUList([prihdu,hduori,hdumsk,hdumean,hdustd,hdurng,hdusig,hdusig2,hdurmn])
    hdulist.writeto(filesout[b],clobber=True)
    
#if __name__ == "__main__":
#    maincode()