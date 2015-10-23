##Fill in blanks in beam

import pyfits
from numpy import zeros,copy,mean,array,isnan,isfinite
from scipy import interpolate as interp
from scipy import where

infiles=[
    '../Inputs/spire_beams_measured/psw_beam_1arcsec_withblanks.fits',\
    '../Inputs/spire_beams_measured/pmw_beam_1arcsec_withblanks.fits',\
    '../Inputs/spire_beams_measured/plw_beam_1arcsec_withblanks.fits']
nf=len(infiles)
outfiles=[
    '../Inputs/spire_beams_measured/psw_beam_1arcsec_noblanks.fits',\
    '../Inputs/spire_beams_measured/pmw_beam_1arcsec_noblanks.fits',\
    '../Inputs/spire_beams_measured/plw_beam_1arcsec_noblanks.fits']
    
for f in range(nf):
    infile=infiles[f]    
    print 'Reading file: %s'%infile
    beamin=pyfits.getdata(infile,1)
    hdr=pyfits.getheader(infile,1)
    nxin=hdr.get('NAXIS1')
    nyin=hdr.get('NAXIS2')
    
    beamout=copy(beamin)

    blank_all=where(beamin!=beamin)
    nbl_tot=len(blank_all)
    print 'Total of %i blank pixels in full beam'%nbl_tot
    xctr=1000
    yctr=1000
    crad=500
    print 'Using central square %i pixels wide'%crad*2
    beam_mid=beamin[xctr-crad:xctr+crad,yctr-crad:yctr+crad]
    xmid=range(xctr-crad,xctr+crad)
    ymid=range(yctr-crad,yctr+crad)
    blank_mid=where(beam_mid != beam_mid)
    nbl_mid=len(blank_mid)
    print 'Total of %i blank pixels in centre of beam'%nbl_mid
    nxmid=len(beam_mid[:,0])
    nymid=len(beam_mid[0,:])
    irad=10
    for x in range(nxmid):
        for y in range(nymid):    
            xb=xmid[x]
            yb=ymid[y]
            if beamin[xb,yb] != beamin[xb,yb]:
                #is blank pixel
                print 'Blank pixel at [%i,%i]'%(xb,yb)                

                ##calculate average of surrounding pixels
                vals=array([beamin[xb-1,yb-1],beamin[xb-1,yb],beamin[xb-1,yb+1],\
                    beamin[xb,yb-1],beamin[xb,yb+1],\
                    beamin[xb+1,yb-1],beamin[xb+1,yb],beamin[xb+1,yb+1]])
                vfin=isfinite(vals)
                beamout[xb,yb]=mean(vals[where(vfin)])

                ##INTERPOLATE BEAM
                #beam_sm=beamin[xb-irad:xb+irad,yb-irad:yb+irad]
                #xsm=range(xb-irad,xb+irad)
                #ysm=range(yb-irad,yb+irad)
                #intbm=interp.RectBivariateSpline(xsm,ysm,beam_sm)
                #beamout[xb,yb]=intbm(xb,yb)
    
    print 'Writing to file: %s'%outfiles[f]
    pyfits.writeto(outfiles[f],beamout,header=hdr,clobber=True)
    