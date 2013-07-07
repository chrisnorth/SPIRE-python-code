def rmdg(map,hdr,ra_z,dec_z):

    import pyfits
    import pywcs as wcs
    
    wcs_map=wcs.WCS(hdr)
    pix_z0=wcs.wcs_sky2pix(ra_z,dec_z,0)
    
    nxin=hdr.get('NAXIS2')
    nyin=hdr.get('NAXIS1')
    cpx=hdr.get('CRPIX2')
    cpy=hdr.get('CRPIX1')
    cpxval=hdr.get('CRVAL2')
    cpyval=hdr.get('CRVAL1')
    dx=hdr.get('CDELT2')*3600. #deg->arcsec
    dy=hdr.get('CDELT1')*3600. #deg->arcsec
    
    