def getbeam(beamtype,band,bpix,bwid,bzlim=None,bzval=0.,verbose=None):

    ## Get beam maps from FITS files, regridded on given grid
    ## Uses the 2D beam maps
    
    import sys
    from numpy import zeros,array,arange,sqrt,exp,log,pi,min,max,floor
    import pyfits
    from scipy import where
    import scipy as sci
    import scipy.interpolate as interp
    #import matplotlib
    #import matplotlib.pyplot as plot

    npx=int(bwid/bpix)
    xl=range(npx)
    yl=range(npx)
    xr=(array(xl)) - floor(npx/2.)
    yr=(array(yl)) - floor(npx/2.)
    
    beamarr=zeros((npx,npx),dtype=float)
    xarr=zeros((npx,npx),dtype=float)
    yarr=zeros((npx,npx),dtype=float)
    for x in xl:
        xarr[x,:]=xr[x]
    for y in yl:
        yarr[:,y]=yr[y]
    negx=False
    negy=False

    if beamtype == 'G':
        #symmetrical Gaussian
        if verbose: print 'Generating Symmetrical Gaussian beam'
        fwhm_all=array((17.6,23.9,35.1)) #geometrical mean
        width=fwhm_all[band]/sqrt(8*log(2))
        print 'Beam FWHM: %.2f arcsec' % fwhm_all[band]
        for x in xl:
            beamarr[x,:]=exp(-(xr[x]**2)/(2.*width**2))
        for y in yl:
            beamarr[:,y]=beamarr[:,y]*exp(-(yr[y]**2)/(2.*width**2))
            
    elif beamtype == 'E':
        #elliptical Gaussian
        if verbose: print 'Generating Elliptical Gaussian beam'
        fwhm_x=array((18.3,24.7,36.9)) #major axis
        fwhm_y=array((17.0,23.2,33.4)) #minor axis
        width_x=fwhm_x[band]/sqrt(8*log(2))
        width_y=fwhm_y[band]/sqrt(8*log(2))
        print 'Beam FWHM: %.2f x %.2f arcsec' % (fwhm_x[band],fwhm_y[band])
        for x in xl:
            beamarr[x,:]=exp(-(xr[x]**2)/(2.*width_x**2))
        for y in yl:
            beamarr[:,y]=beamarr[:,y]*exp(-(yr[y]**2)/(2.*width_y**2))
        
        
    elif beamtype == 'M':
        #measured beam
        if band == 0:
            #file='../Inputs/spire_beams_measured/psw_beam_1arcsec.fits'
            file='../CalProducts/0x5000241aL_PSW_pmcorr_1arcsec_cln_bgsub.fits'
        elif band == 1:
            #file='../Inputs/spire_beams_measured/pmw_beam_1arcsec.fits'
            file='../CalProducts/0x5000241aL_PMW_pmcorr_1arcsec_cln_bgsub.fits'
        elif band == 2:
            #file='../Inputs/spire_beams_measured/plw_beam_1arcsec.fits'
            file='../CalProducts/0x5000241aL_PLW_pmcorr_1arcsec_cln_bgsub.fits'
        if verbose: print 'Generating measured beam from %s'%file
        #read in data
        beamin=pyfits.getdata(file,0)
        #Set NaNs to 0
        beamin=sci.where(beamin==beamin,beamin,0.)
        if bzlim != None:
            beamin=sci.where(beamin<=bzlim,bzval,beamin)
        #read in header data
        hdr=pyfits.getheader(file,0)
        nxin=hdr.get('NAXIS2')
        nyin=hdr.get('NAXIS1')
        cpx=hdr.get('CRPIX2')
        cpy=hdr.get('CRPIX1')
        dx=hdr.get('CDELT1')*3600. #deg->arcsec
        dy=hdr.get('CDELT1')*3600. #deg->arcsec
        #print 'FITS header centre: %d,%d'%(cpx,cpy)
        ##Hard-coded locations of maximums from FITS files
        if band == 0:
            #cpx=1000
            #cpy=999
            cpy=965
            cpx=1118
        elif band == 1:
            #cpx=1005
            #cpy=1004
            cpy=965
            cpx=1138
        elif band == 2:
            #cpx=1013
            #cpy=1013
            cpy=976
            cpx=1139
        #print 'hard-coded centre: %d,%d'%(cpx,cpy)
        #print 'empirical centre: ',where(beamin == 1.)
        #interpolation requires monotically increasing grids
        if dx < 0: negx=True
        if dy < 0: negy=True
        xin=abs(dx)*(arange(0,nxin) - cpx)
        yin=abs(dy)*(arange(0,nyin) - cpy)
        xarrin=zeros((nxin,nyin))
        for x in range(0,nxin):
            xarrin[x,:]=float(x) - cpx
        yarrin=zeros((nxin,nyin))
        for y in range(0,nyin):
            yarrin[:,y]=float(y) - cpy

        #interpolate beam
        print xin.shape
        print yin.shape
        print beamin.shape
        intbm=interp.RectBivariateSpline(xin,yin,beamin)

        beamarr=intbm(xr,yr)

#        plot.figure(5+band)
#        plot.contourf(xarrin,yarrin,log10(beamin))
#        plot.colorbar()

#        plot.figure(8+band)
#        plot.contour(xarr,yarr,log10(beamarr))
#        plot.colorbar()
       
 
    elif beamtype == 'T':
        #measured beam
        if verbose: print 'Generating theoretical beam'
        if band == 0:
            file='../Inputs/spire_beams_theoretical/spire_psw_flight_psf.fits'
        elif band == 1:
            file='../Inputs/spire_beams_theoretical/spire_pmw_flight_psf.fits'
        elif band == 2:
            file='../Inputs/spire_beams_theoretical/spire_plw_flight_psf.fits'
        
        #read in data
        beamin=pyfits.getdata(file,0)
        #set NaNs to 0
        beamin=sci.where(beamin==beamin,beamin,0.)

        if bzlim != None:
            beamin=sci.where(beamin<=bzlim,bzval,beamin)

        #read in header data
        hdr=pyfits.getheader(file,0)
        nxin=hdr.get('NAXIS1')
        nyin=hdr.get('NAXIS2')
        cpx=hdr.get('CRPIX1')
        cpy=hdr.get('CRPIX2')
        dx=hdr.get('CDELT1')*3600. #deg->arcsec
        dy=hdr.get('CDELT1')*3600. #deg->arcsec

        if dx < 0: negx=True
        if dy < 0: negy=True
        xin=abs(dx)*(array(range(0,nxin)) - cpx)
        yin=abs(dy)*(array(range(0,nyin)) - cpy)
        xarrin=zeros((nxin,nyin))
        for x in range(0,nxin):
            xarrin[x,:]=float(x) - cpx
        yarrin=zeros((nxin,nyin))
        for y in range(0,nyin):
            yarrin[:,y]=float(y) - cpy

        #interpolate beam
        intbm=interp.RectBivariateSpline(xin,yin,beamin)

        beamarr=intbm(xr,yr)

        beamarr=where(xarr < min(xin),0.,beamarr)
        beamarr=where(xarr > max(xin),0.,beamarr)
        beamarr=where(yarr < min(yin),0.,beamarr)
        beamarr=where(yarr > max(yin),0.,beamarr)
#        plot.figure(5+band)
#        plot.contourf(xarrin,yarrin,log10(beamin))
#        plot.colorbar()

#        plot.figure(8+band)
#        plot.contour(xarr,yarr,log10(beamarr))
#        plot.colorbar()
       
 
    else:
        error='ERROR: Unknown beam type "%s". For beam map, must be one of G|E|M|T.' % beamtype
        sys.exit(error)

####Don't truncate beam
#    #Set max radius
#    if brad != 0.:
#        beamarr=sci.where(sqrt(xarr**2 + yarr**2) > brad,0.,beamarr)

    if negx == True:
        for x in xl:
            xarr[x,:]=xr[x]
    if negy == True:
        for y in yl:
            xarr[y,:]=yr[y]
    

    return(xarr,yarr,beamarr)

    
def getbeamprofile(beamtype,band,radarr,bzlim=None,bzval=0.,verbose=None):

    ## Read in beam profile for all three bands
    ## Uses a single 1d beam profile
    ## Superseded by getnewbeamprofile
    
    from numpy import zeros,array,arange,max,sqrt,log,exp
    from scipy import where
    from scipy.interpolate import interp1d
    import sys

    nrad=radarr.size    
    beam_sm=zeros(nrad)

    if beamtype == 'G':
        if verbose: print 'Generating Symmetrical Gaussian beam profile'
        fwhm_all=array((17.6,23.9,35.1)) #geometrical mean
        width=fwhm_all[band]/sqrt(8*log(2))
        print 'Beam FWHM: %.2f arcsec' % fwhm_all[band]
        for r in nrad:
            beam_sm[r]=exp(-(radarr[r]**2)/(2.*width**2))
            
    else:
        if beamtype == 'S':
            files=['../Inputs/spire_beams_spliced/psw_beamprofile_1arcsec_spliced.txt',
                   '../Inputs/spire_beams_spliced/pmw_beamprofile_1arcsec_spliced.txt',
                   '../Inputs/spire_beams_spliced/plw_beamprofile_1arcsec_spliced.txt']
               
        elif beamtype == 'M':
            files=['../Inputs/spire_beams_measured/psw_beamprofile_1arcsec.txt',
                   '../Inputs/spire_beams_measured/pmw_beamprofile_1arcsec.txt',
                   '../Inputs/spire_beams_measured/plw_beamprofile_1arcsec.txt']

        elif beamtype == 'T':
            files=['../Inputs/spire_beams_theoretical/psw_beamprofile_1arcsec.txt',
                   '../Inputs/spire_beams_theoretical/pmw_beamprofile_1arcsec.txt',
                   '../Inputs/spire_beams_theoretical/plw_beamprofile_1arcsec.txt']
        else:
            error='ERROR: Unknown beam type "%s". For beam profile, must be one of G|M|T|S.' % beamtype
            sys.exit(error)

    
        f = open(files[band], 'r')
        lines = array(f.readlines())
        nlines=lines.size
        #print nlines
        radin=arange(nlines)
        beamin=zeros(nlines)
        maxrad=radin[nlines-1]
        maxradarr=radarr[nrad-1]
        #print maxrad,maxradarr
        for i in range(nlines):
            beamin[i]=(float(lines[i]))
            beamint=interp1d(radin,beamin)
        if maxradarr <= maxrad:
            beam_sm=beamint(radarr)
        else:
            radcrop=where(radarr > maxrad,maxrad,radarr)
            beam_sm=beamint(radcrop)
            beam_sm=where(radarr > maxrad,0.,beam_sm)
    
    if bzlim != None:
        beam_sm=where(beam_sm < bzlim,bzval,beam_sm)
    
        
    return(beam_sm)

def getnewbeamprofile(band,radarr,bzlim=None,bzval=0.,verbose=None,type='M'):

    ## Read in beam profile for given band
    ## Uses two beam profile files (one core, one outer)
    
    from numpy import zeros,array,arange,max,sqrt,log,exp
    from scipy import where
    from scipy.interpolate import interp1d
    import sys

    nrad=radarr.size    
    beam_scl=zeros(nrad)
    beam_fix=zeros(nrad)

    if type=='M':
        files_scl=['../Inputs/PSW_Core.csv',
                   '../Inputs/PMW_Core.csv',
                   '../Inputs/PLW_Core.csv']    
                   
        files_fix=['../Inputs/PSW_Constant.csv',
                   '../Inputs/PMW_Constant.csv',
                   '../Inputs/PLW_Constant.csv']      
    elif type=='T':
        files_scl=['../Inputs/spire_beams_theoretical/psw_beamprofile_1arcsec.txt',
                   '../Inputs/spire_beams_theoretical/pmw_beamprofile_1arcsec.txt',
                   '../Inputs/spire_beams_theoretical/plw_beamprofile_1arcsec.txt']
                   
        files_fix=['../Inputs/spire_beams_theoretical/beamprofile_1arcsec_blank.txt',
                   '../Inputs/spire_beams_theoretical/beamprofile_1arcsec_blank.txt',
                   '../Inputs/spire_beams_theoretical/beamprofile_1arcsec_blank.txt']
    
    fscl = open(files_scl[band], 'r')
    ffix = open(files_fix[band], 'r')
    lscl = array(fscl.readlines())
    lfix = array(ffix.readlines())
    nlines=lscl.size
    #print nlines
    radin=arange(nlines)
    beamin_scl=zeros(nlines)
    beamin_fix=zeros(nlines)
    maxrad=radin[nlines-1]
    maxradarr=radarr[nrad-1]
    #print maxrad,maxradarr
    for i in range(nlines):
        #beamin_scl[i]=(10**float(lscl[i]))
        #beamin_fix[i]=(10**float(lfix[i]))
        beamin_scl[i]=float(lscl[i])
        beamin_fix[i]=float(lfix[i])
    if verbose:
        print 'Input Scaled range: %g-%g'%(min(beamin_scl),max(beamin_scl))
        print 'Input Fixed range:  %g-%g'%(min(beamin_fix),max(beamin_fix))
    beamint_scl=interp1d(radin,beamin_scl)
    beamint_fix=interp1d(radin,beamin_fix)
    if maxradarr <= maxrad:
        beam_scl=beamint_scl(radarr)
        beam_fix=beamint_fix(radarr)
    else:
        radcrop=where(radarr > maxrad,maxrad,radarr)
        beam_scl=beamint_scl(radcrop)
        beam_fix=beamint_fix(radcrop)
        beam_scl=where(radarr > maxrad,0.,beam_scl)
        beam_fix=where(radarr > maxrad,0.,beam_fix)
        if verbose:
            print 'Cropping at %f arcsec'%(maxrad)
            print radcrop

    if bzlim != None:
        beam_scl=where(beam_scl < bzlim,bzval,beam_scl)
        beam_fix=where(beam_fix < bzlim,bzval,beam_fix)
        if verbose:
            print 'Using zero-limit of %g -> %f'%(bzlim,bzval)
            print 'Output Scaled range: %g-%g'%(min(beam_scl),max(beam_scl))
            print 'Output Fixed range:  %g-%g'%(min(beam_fix),max(beam_fix))
    return(beam_scl,beam_fix)

def getbeammaps(src='Nep',regrid=None):

    from numpy import array,zeros,zeros_like,mean,max
    from scipy import where
    import pyfits
    
    crpix=[[0,0],[0,0],[0,0]]
    if src=='Nep':
        fitsfiles=['../CalProducts/0x5000241aL_PSW_pmcorr_1arcsec_cln_bgsub.fits', \
            '../CalProducts/0x5000241aL_PMW_pmcorr_1arcsec_cln_bgsub.fits', \
            '../CalProducts/0x5000241aL_PLW_pmcorr_1arcsec_cln_bgsub.fits']
        crpix=[[965,1118],[965,1138],[976,1139]]
    elif src=='Theory':
        fitsfiles=['../Inputs/spire_beams_theoretical/spire_psw_flight_psf.fits', \
            '../Inputs/spire_beams_theoretical/spire_pmw_flight_psf.fits', \
            '../Inputs/spire_beams_theoretical/spire_plw_flight_psf.fits']
    elif src=="Mars1":
        fitsfiles=['../Inputs/Mars_5000532c/1arcsec/Mars_mapPSW_5000532c_1arcsec-cutTrim.fits', \
            '../Inputs/Mars_5000532c/1arcsec/Mars_mapPMW_5000532c_1arcsec-cutTrim.fits', \
            '../Inputs/Mars_5000532c/1arcsec/Mars_mapPLW_5000532c_1arcsec-cutTrim.fits']
    elif src=="Mars2":
        fitsfiles=['../Inputs/Mars_50004e4b/1arcsec/Mars_mapPSW_50004e4b_1arcsec-cutTrim.fits', \
            '../Inputs/Mars_50004e4b/1arcsec/Mars_mapPMW_50004e4b_1arcsec-cutTrim.fits', \
            '../Inputs/Mars_50004e4b/1arcsec/Mars_mapPLW_50004e4b_1arcsec-cutTrim.fits']
    xarr_map=[]
    yarr_map=[]
    beam_map=[]

    if (regrid!=None):
        if (len(regrid)!=3):
            regrid=[6,10,14]
        #print 'regrid:',regrid
    for b in range(3):
        #print '---\nReading in beam map for band %d'%(b+1)
        file=fitsfiles[b]
        beamin=pyfits.getdata(file,0)
        beamin=where(beamin==beamin,beamin,0.)
        hdr=pyfits.getheader(file,0)
        nxin=hdr.get('NAXIS1')
        nyin=hdr.get('NAXIS2')
        dx=hdr.get('CDELT1')*3600. #deg->arcsec
        dy=hdr.get('CDELT1')*3600. #deg->arcsec
        if src=="Nep":
            cpx=crpix[b][0]
            cpy=crpix[b][1]
        else:
            cpx=hdr.get('CRVAL1')
            cpy=hdr.get('CRVAL2')
        xin=abs(dx)*(array(range(0,nxin)) - cpx)
        yin=abs(dy)*(array(range(0,nyin)) - cpy)
        xarrin=zeros((nyin,nxin))
        for x in range(0,nxin):
            xarrin[:,x]=float(x) - cpx
        yarrin=zeros((nyin,nxin))
        for y in range(0,nyin):
            yarrin[y,:]=float(y) - cpy
        #print '(%d,%d):%.4f'%(xarrin[where(beamin==max(beamin))],yarrin[where(beamin==max(beamin))],beamin[where(beamin==max(beamin))])
        xz=where(xin==0)[0]
        yz=where(yin==0)[0]
        #print 'beam[xmax,ymax]:',beamin[yz,xz]
        #print 'max(beam)',max(beamin)
    
        if (regrid==None):
            xarr_map.append(xarrin)
            yarr_map.append(yarrin)
            beam_map.append(beamin)
        else:            
            xr=720
            yr=xr
            #beamin=beam_map[b]
            #xarrin=xarr_map[b]
            #yarrin=yarr_map[b]
            cpx=crpix[b][0]
            cpy=crpix[b][1]
            ixlim=[cpx-xr-regrid[b]/2,cpx+xr+regrid[b]/2]
            iylim=[cpy-yr-regrid[b]/2,cpy+yr+regrid[b]/2]
            ixout=range(ixlim[0],ixlim[1],regrid[b])
            iyout=range(iylim[0],iylim[1],regrid[b])    
            nxout=len(ixout)
            nyout=len(iyout)
            beamout=zeros((nxout,nyout))
            xarrout=zeros_like(beamout)
            yarrout=zeros_like(beamout)
            for x in range(len(ixout)):
                for y in range(len(iyout)):
                    x0=ixout[x]
                    x1=x0+regrid[b]
                    y0=iyout[y]
                    y1=y0+regrid[b]
            
                    beampatch=beamin[x0:x1,y0:y1]
                    beamout[x,y]=mean(beamin[x0:x1,y0:y1])
                    xarrout[x,y]=mean(xarrin[x0:x1,y0:y1])
                    yarrout[x,y]=mean(yarrin[x0:x1,y0:y1])
                    #if max(beampatch) > 0.9:
                        #print 'max(%d:%d,%d:%d)=%.4g ; (%.1f,%.1f)=%.4g)'%(x0,x1,y0,y1,max(beampatch),xarrout[x,y],yarrout[x,y],beamout[x,y])
            #print 'max(beamin): %.5g'%max(beamin)
            #print 'max(beamout): %.5g'%max(beamout)
            beamout=beamout/max(beamout)
            xarr_map.append(xarrout)
            yarr_map.append(yarrout)
            beam_map.append(beamout)
    
    xarr_map=array(xarr_map)
    yarr_map=array(yarr_map)
    beam_map=array(beam_map)
    return(xarr_map,yarr_map,beam_map)
        
def getbeammaps_theory(regrid=None,maxxy=999):

    from numpy import array,zeros,zeros_like,mean,max
    from scipy import where
    import pyfits
    
    fitsfiles=['../Inputs/spire_beams_theoretical/spire_psw_flight_psf.fits', \
            '../Inputs/spire_beams_theoretical/spire_pmw_flight_psf.fits', \
            '../Inputs/spire_beams_theoretical/spire_plw_flight_psf.fits']
    
    xarr_map=[]
    yarr_map=[]
    beam_map=[]

    if (regrid!=None):
        if (len(regrid)!=3):
            regrid=[6,10,14]
        #print 'regrid:',regrid
    for b in range(3):
        #print '---\nReading in beam map for band %d'%(b+1)
        file=fitsfiles[b]
        beamin=pyfits.getdata(file,0)
        beamin=where(beamin==beamin,beamin,0.)
        hdr=pyfits.getheader(file,0)
        nxin=hdr.get('NAXIS1')
        nyin=hdr.get('NAXIS2')
        dx=hdr.get('CDELT1')*3600. #deg->arcsec
        dy=hdr.get('CDELT1')*3600. #deg->arcsec
        cpx=hdr.get('CRPIX1')
        cpy=hdr.get('CRPIX2')
        xin=abs(dx)*(array(range(0,nxin)) - cpx)
        yin=abs(dy)*(array(range(0,nyin)) - cpy)
        xarrin=zeros((nyin,nxin))
        for x in range(0,nxin):
            xarrin[:,x]=float(x) - cpx
        yarrin=zeros((nyin,nxin))
        for y in range(0,nyin):
            yarrin[y,:]=float(y) - cpy
        #print '(%d,%d):%.4f'%(xarrin[where(beamin==max(beamin))],yarrin[where(beamin==max(beamin))],beamin[where(beamin==max(beamin))])
        xz=where(xin==0)[0]
        yz=where(yin==0)[0]
        #print 'beam[xmax,ymax]:',beamin[yz,xz]
        #print 'max(beam)',max(beamin)
    
        if (regrid==None):
            xarr_map.append(xarrin)
            yarr_map.append(yarrin)
            beam_map.append(beamin)
        else:            
            xr=720
            yr=xr
            #beamin=beam_map[b]
            #xarrin=xarr_map[b]
            #yarrin=yarr_map[b]
            #cpx=crpix[b][0]
            #cpy=crpix[b][1]
            ixlim=[cpx-xr-regrid[b]/2,cpx+xr+regrid[b]/2]
            iylim=[cpy-yr-regrid[b]/2,cpy+yr+regrid[b]/2]
            ixout=range(ixlim[0],ixlim[1],regrid[b])
            iyout=range(iylim[0],iylim[1],regrid[b])    
            nxout=len(ixout)
            nyout=len(iyout)
            beamout=zeros((nxout,nyout))
            xarrout=zeros_like(beamout)
            yarrout=zeros_like(beamout)
            for x in range(len(ixout)):
                for y in range(len(iyout)):
                    x0=ixout[x]
                    x1=x0+regrid[b]
                    y0=iyout[y]
                    y1=y0+regrid[b]
            
                    beampatch=beamin[x0:x1,y0:y1]
                    beamout[x,y]=mean(beamin[x0:x1,y0:y1])
                    xarrout[x,y]=mean(xarrin[x0:x1,y0:y1])
                    yarrout[x,y]=mean(yarrin[x0:x1,y0:y1])
                    #if max(beampatch) > 0.9:
                        #print 'max(%d:%d,%d:%d)=%.4g ; (%.1f,%.1f)=%.4g)'%(x0,x1,y0,y1,max(beampatch),xarrout[x,y],yarrout[x,y],beamout[x,y])
            #print 'max(beamin): %.5g'%max(beamin)
            #print 'max(beamout): %.5g'%max(beamout)
            beamout=beamout/max(beamout)
            xarr_map.append(xarrout)
            yarr_map.append(yarrout)
            beam_map.append(beamout)
            
    xarr_map=array(xarr_map)
    yarr_map=array(yarr_map)
    beam_map=array(beam_map)
    return(xarr_map,yarr_map,beam_map)


def comb_beam(beam_scl,beam_fix):

    ##Produce a combined beam file
    ##Takes core and constant beam maps

    from scipy import where
    beam_cmb=where(beam_fix > beam_scl,beam_fix,beam_scl)

    return(beam_cmb)
    
def beam_regrid(beamx,beamy,beamin,nu0,nu,ind=1.0):
    
    ## Regrid the beam based on a new frequency
    ## Uses 2d beam maps
    
    import scipy.interpolate as interp
    
    xin=beamx[:,0]
    yin=beamy[0,:]
    
    mult=(nu/nu0)**ind
    
    #generate new grid
    
    xout = xin * mult
    yout = yin * mult
    
    #make interpolation object for beam map
    intbeam=interp.RectBivariateSpline(xin,yin,beamin)
    
    beamout=intbeam(xout,yout)

    return(beamout)

def beam_regrid_xy(intbeam,xin,yin,nu0,nu,ind=1.0):
    
    #uses interpolation object to regrid beam map
    
    #import scipy.interpolate as interp
    
    mult=(nu/nu0)**ind
    
    xout = xin * mult
    yout = yin * mult
    
    beamout=intbeam(xout,yout)

    return(beamout)
    
def beamarea(beamx,beamy,beamin,brad=None,pix=1.):
    
    ## Calculate the effective area of a 2-d gridded beam over a fully extended source
    ## Requires pixels to be 1arcsec on a side, unless pix keyword is set
    
    from numpy import sqrt,sum,shape
    from scipy import where
   
    beamrad=sqrt(beamx**2 + beamy**2)
    
    if brad != None:
        beam_crop=where(beamrad <= brad,beamin,0.)
        area=sum(beam_crop)
    else:
        area=sum(beamin)
    
    area=area * pix**2
    return(area)
    
def beamarea_th(beamx,beamy,beamin,theta):

    ## Calculate the effective area of a 2-d gridded beam over a source
    ## source profile is Gaussian of FWHM theta
    ##Requires pixels to be 1arcmin on a side
    
    from numpy import sqrt,sum,log,exp
    
    rad=sqrt(beamx**2 + beamy**2)
    swid=theta/(2.*sqrt(log(2.)))
    
    src=exp(-(rad/swid)**2)
    area=sum(beamin*src)
    
    return(area)

def beamarea_az(radarr,beamaz,brad=None):
    
    ## Calculate the effective beam area of an azimuthally-smoothed beam
    ## Requires pixels to be 1arcmin on a side
    
    from numpy import pi
    from scipy import where
    from scipy.integrate import trapz

    if brad != None:
        beamcrop=where(radarr > brad,0.,beamaz)
        area=trapz(beamcrop*2.*pi*radarr,radarr)
    else:
        area=trapz(beamaz*2.*pi*radarr,radarr)
    
    return(area)

def beamarea_az_th(radarr,beamaz,theta):
    
    ## Calculate the effective beam area of an azimuthally-smoothed beam over a source
    ## source is Gaussian of FWHM theta
    ## Requires pixels to be 1arcmin on a side
    
    from numpy import sqrt,log,exp,pi
    from scipy.integrate import trapz
    
    swid=theta/(sqrt(8.*log(2.)))
    ##swid is in Gaussian widths
    
    src=exp(-radarr**2/(2.*swid**2))
    
    area=trapz(beamaz*2.*pi*radarr*src,radarr)
    
    return(area)

def beam_azsm(beamx,beamy,beamin,radarr,retall=None,pess=None,nsamp=360):

    ## Azimuthally-smooth a 3d beam map to produce a radial beam profile
    
    from numpy import zeros,arange,pi,cos,sin,average,min,max,std
    import scipy.interpolate as interpolate
    from scipy import where
    
    nrad=radarr.size
    interp=interpolate.RectBivariateSpline(beamx[:,0],beamy[0,:],beamin)
    
    beam_sm=zeros(nrad)
    beam_max=zeros(nrad)
    beam_min=zeros(nrad)
    beam_sd=zeros(nrad)
    for r in arange(nrad):
        thetalist=arange(nsamp)*(pi/180.) * (360./nsamp)# full circle in radians
        xlist = radarr[r] * cos(thetalist)
        ylist = -radarr[r] * sin(thetalist)
        beam_r=interp.ev(xlist,ylist)
        beam_sm[r]=average(beam_r)
        beam_min[r]=min(beam_r)
        beam_max[r]=max(beam_r)
        beam_sd[r]=std(beam_r)

    if pess:
        inmap=where(beamin > 0)
        maxx=max(abs(beamx[inmap]))
        maxy=max(abs(beamy[inmap]))
        maxrad=min((maxx,maxy))
        print 'Limiting beam profile to %d arcsec'%(maxrad)
        beam_sm=where(radarr > maxrad,0.,beam_sm)
    
    if retall == True:
        return(beam_sm,beam_min,beam_max,beam_sd)
    else:
        return(beam_sm)

def beam_azsm2(beamx,beamy,beamin,radarr,retall=None,pess=None,nsamp=360):

    ## Azimuthally-smooth the beam (alternative method)
    
    from numpy import zeros,arange,pi,cos,sin,average,min,max,std,reshape,array
    import scipy.interpolate as interpolate
    
    nrad=radarr.size
    npt=beamx.size
    x1=reshape(beamx,npt)
    y1=reshape(beamy,npt)
    beamin1=reshape(beamin,npt)
    
    pos1=zeros((npt,2))
    pos1[:,0]=x1
    pos1[:,1]=y1
    
    #interp=interpolate.interp2d(x1,y1,beamin1)
    interp=interpolate.LinearNDInterpolator(pos1,beamin1)
    #interp=interpolate.RectBivariateSpline(beamx[:,0],beamy[0,:],beamin)
    

    beam_sm=zeros(nrad)
    beam_max=zeros(nrad)
    beam_min=zeros(nrad)
    beam_sd=zeros(nrad)
    for r in arange(nrad):
        thetalist=arange(nsamp)*pi/180. * (360/nsamp) # full circle in radians
        xlist = radarr[r] * cos(thetalist)
        ylist = -radarr[r] * sin(thetalist)
        beam_r=interp.ev(xlist,ylist)
        beam_sm[r]=average(beam_r)
        beam_min[r]=min(beam_r)
        beam_max[r]=max(beam_r)
        beam_sd[r]=std(beam_r)
    
    if retall == True:
        return(beam_sm,beam_min,beam_max,beam_sd)
    else:
        return(beam_sm)

def measbeam(radarr,beam_sm,nuarr,rsrfarr,nu0,alpha_nep,ind=1.0):

    ## Calculate the measured beam area over Neptune
    ## Based on given effective frequency and Neptune spectral index
    ## Uses single radial beam profile

    from numpy import zeros,arange,pi
    import scipy
    import scipy.interpolate as interp
    import scipy.integrate as integrate
    
    nrad=radarr.size
    nnu=nuarr.size
    beammeas=zeros(nrad)
    beam_sm_nu=zeros(nnu)
    
    beam_sm_int=interp.interp1d(radarr,beam_sm)
    denom_arg = rsrfarr * nuarr**alpha_nep
    for r in arange(nrad):
        r_nu=radarr[r]*(nuarr/nu0)**ind
        inr=scipy.where(r_nu < max(radarr))
        outr=scipy.where(r_nu >= max(radarr))
        beam_sm_nu[inr]=beam_sm_int(r_nu[inr])
        beam_sm_nu[outr]=0.
        num_arg = rsrfarr * nuarr**alpha_nep * beam_sm_nu
        beammeas[r]=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)
    
    #calculate beam area    
    areameas=integrate.trapz(beammeas*2.*pi*radarr,radarr)

    return(areameas)            
    
def measbeam_new(radarr,beam_scl,beam_fix,nuarr,rsrfarr,nu0,alpha,ind=1.0,brad=None):

    ## Calculate the measured beam area over source of spectral index alpha
    ## based on given effective frequency and Neptune spectral index
    ## sing the scaled and fixed beams

    from numpy import zeros,arange,pi,max
    import scipy
    import scipy.interpolate as interp
    import scipy.integrate as integrate
    
    nrad=radarr.size
    nnu=nuarr.size
    beammeas=zeros(nrad)
    beam_scl_nu=zeros(nnu)
    beam_fix_nu=zeros(nnu)
    
    beam_scl_int=interp.interp1d(radarr,beam_scl)
    beam_fix_int=interp.interp1d(radarr,beam_fix)
    denom_arg = rsrfarr * nuarr**alpha
    for r in arange(nrad):
        ##integrate over the band
        r_nu=radarr[r]*(nuarr/nu0)**ind ##effective radius over the frequency range
        inr=scipy.where(r_nu < max(radarr))
        outr=scipy.where(r_nu >= max(radarr))
        beam_scl_nu[inr]=beam_scl_int(r_nu[inr])
        beam_scl_nu[outr]=0.
        beam_fix_nu[inr]=beam_fix_int(r_nu[inr])
        beam_fix_nu[outr]=0.
        
        ##add on fixed beam
        beam_cmb_nu=scipy.where(beam_scl_nu > beam_fix[r],beam_scl_nu,beam_fix[r])
        num_arg = rsrfarr * nuarr**alpha * beam_cmb_nu
        beammeas[r]=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)
    
    if brad==None:
        brad=max(radarr)
    beammeas_int=scipy.where(radarr <= brad,beammeas,0.)
    #calculate beam area
    areameas=integrate.trapz(beammeas_int*2.*pi*radarr,radarr)

    return(areameas)            

def measbeam_new_wrong(radarr,beam_scl,beam_fix,nuarr,rsrfarr,nu0,alpha,ind=1.0,brad=None):

    ## Calculate the measured beam area over source of spectral index alpha
    ## based on given effective frequency and Neptune spectral index
    ## sing the scaled and fixed beams

    from numpy import zeros,arange,pi,max
    import scipy
    import scipy.interpolate as interp
    import scipy.integrate as integrate
    
    nrad=radarr.size
    nnu=nuarr.size
    beammeas=zeros(nrad)
    beam_scl_nu=zeros(nnu)
    beam_fix_nu=zeros(nnu)
    
    beam_scl_int=interp.interp1d(radarr,beam_scl)
    beam_fix_int=interp.interp1d(radarr,beam_fix)
    denom_arg = rsrfarr * nuarr**alpha
    for r in arange(nrad):
        ##integrate over the band
        #################### THIS STEP IS DELIBERATELY WRONG ###############
        r_nu=radarr[r]/(nuarr/nu0)**ind ##effective radius over the frequency range
        ####################################################################
        inr=scipy.where(r_nu < max(radarr))
        outr=scipy.where(r_nu >= max(radarr))
        beam_scl_nu[inr]=beam_scl_int(r_nu[inr])
        beam_scl_nu[outr]=0.
        beam_fix_nu[inr]=beam_fix_int(r_nu[inr])
        beam_fix_nu[outr]=0.
        
        ##add on fixed beam
        beam_cmb_nu=scipy.where(beam_scl_nu > beam_fix[r],beam_scl_nu,beam_fix[r])
        num_arg = rsrfarr * nuarr**alpha * beam_cmb_nu
        beammeas[r]=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)
    
    if brad==None:
        brad=max(radarr)
    beammeas_int=scipy.where(radarr <= brad,beammeas,0.)
    #calculate beam area
    areameas=integrate.trapz(beammeas_int*2.*pi*radarr,radarr)

    return(areameas)            

def measbeam_simple(area_nu0,nuarr,rsrfarr,nuc,alpha,ind=0.85):
    
    import scipy.integrate as integrate
    
    num_arg = rsrfarr * nuarr**alpha * area_nu0 * (nuarr/nuc)**(-2*ind)
    denom_arg = rsrfarr * nuarr**alpha
    
    area_eff = integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)
    
    return(area_eff)

    
def measbeamtemp_new(radarr,beam_scl,beam_fix,nuarr,rsrfarr,nu0,temp,beta,ind=1.0):

    ## Calculate the measured beam area over source of given temperature
    ## based on given effective frequency and Neptune spectral index
    ## sing the scaled and fixed beams

    from numpy import zeros,arange,pi
    import scipy
    import scipy.interpolate as interp
    import scipy.integrate as integrate
    from calc_conv import greybody
    
    nrad=radarr.size
    nnu=nuarr.size
    beammeas=zeros(nrad)
    beam_scl_nu=zeros(nnu)
    beam_fix_nu=zeros(nnu)
    
    beam_scl_int=interp.interp1d(radarr,beam_scl)
    beam_fix_int=interp.interp1d(radarr,beam_fix)
    denom_arg = rsrfarr * greybody(nuarr,temp,beta,nu0)
    for r in arange(nrad):
        ##integrate over the band
        r_nu=radarr[r]*(nuarr/nu0)**ind ##effective radius over the frequency range
        inr=scipy.where(r_nu < max(radarr))
        outr=scipy.where(r_nu >= max(radarr))
        beam_scl_nu[inr]=beam_scl_int(r_nu[inr])
        beam_scl_nu[outr]=0.
        beam_fix_nu[inr]=beam_fix_int(r_nu[inr])
        beam_fix_nu[outr]=0.
        
        ##add on fixed beam
        beam_cmb_nu=scipy.where(beam_scl_nu > beam_fix[r],beam_scl_nu,beam_fix[r])
        num_arg = rsrfarr * greybody(nuarr,temp,beta,nu0) * beam_cmb_nu
        beammeas[r]=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)
    
    #calculate beam area    
    areameast=integrate.trapz(beammeas*2.*pi*radarr,radarr)

    return(areameast)            

def get_effnu(beamx,beamy,beams,rsrfarr,nuarr,nuc,brad,alpha_nep,aprec=0.0001,verbose=None,ind=1.0):
    
    ## Find the effective frequency of the beam, based on measurements of Neptune
    ## Takes 2d beam maps
    
    from numpy import array,zeros,arange    
    from beam import measbeam,beamarea

    band=arange(3)

    area=zeros(3)
    for b in band:
        area[b]=beamarea(beamx,beamy,beams[:,:,b],brad=brad)

    if verbose: print 'Beam Areas (brad): [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])

    #Spectral index of neptune
    #alpha_nep=array((1.26,1.39,1.45))

    a_init=array([[1.,1.01],[1.,1.02],[1.,1.03]])

    maxit=10.
    radarr=arange(brad)
    a_fin=zeros(3)
    ameas_fin=zeros(3)
    arel_fin=zeros(3)

    for b in band:
        if verbose: print 'Band %d'%b
        a_arr=a_init[b,:]
        ameas=zeros(2)
        beam_sm=beam_azsm(beamx,beamy,beams[:,:,b],radarr)
        ameas[0]=measbeam(radarr,beam_sm,nuarr,rsrfarr[:,b],nuc[b]*a_arr[0],alpha_nep[b])
        ameas[1]=measbeam(radarr,beam_sm,nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b])
        adiff=ameas-area[b]
        arel=(ameas-area[b])/area[b]
        if verbose:
            print 'a=%.4f, A=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[0],ameas[0],adiff[0],arel[0])
            print 'a=%.4f, A=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
        it=1
        done=False
        while done==False:
            grad=(adiff[1]-adiff[0])/(a_arr[1]-a_arr[0])
            a_arrnew=a_arr[1] - adiff[1]/grad
            a_arr[0]=a_arr[1]
            ameas[0]=ameas[1]
            a_arr[1]=a_arrnew
            ameas[1]=measbeam(radarr,beam_sm,nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b],ind=ind)
            adiff=ameas-area[b]        
            arel=(ameas-area[b])/area[b]
            if verbose:            
                print 'a=%.4f, Area=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
            it=it+1
            if abs(arel[1]) < aprec:
                done=True
            elif it >= maxit:
                done=True
        a_fin[b]=a_arr[1]
        ameas_fin[b]=ameas[1]
        arel_fin[b]=arel[1]
        if verbose:        
            print 'Final a=%.4f, Area=%.2f (Rel. error=%.6f)' % (a_fin[b],ameas_fin[b],arel_fin[b])
            print 'Nu(eff)=%.2f GHz'%(nuc[b]*a_fin[b]/1.e9)

    if verbose:
        print 'Meas areas: [%.2f, %.2f, %.2f]' % (ameas_fin[0],ameas_fin[1],ameas_fin[2])
        print 'True areas: [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])

    nueff=nuc*a_fin
    print 'Effective frequencies: [%.2f, %.2f, %.2f] GHz' % (nueff[0]/1.e9,nueff[1]/1.e9,nueff[2]/1.e9)
    return(nueff)
    
def get_effnu_az(radarr,beam_sm,rsrfarr,nuarr,nuc,brad,alpha_nep,aprec=0.0001,verbose=None,ind=1.0):
    
    ## Find the effective frequency of the beam, based on measurements of Neptune
    ## Takes 1d beam profile
    
    from numpy import array,zeros,arange    
    from beam import measbeam,beamarea_az

    band=arange(3)

    area=zeros(3)
    for b in band:
        area[b]=beamarea_az(radarr,beam_sm[:,b],brad=brad)

    if verbose: print 'Beam Areas (brad): [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])

    #Spectral index of neptune
    #alpha_nep=array((1.26,1.39,1.45))

    a_init=array([[1.,1.01],[1.,1.02],[1.,1.03]])

    maxit=10.
    radarr=arange(brad)
    a_fin=zeros(3)
    ameas_fin=zeros(3)
    arel_fin=zeros(3)

    for b in band:
        if verbose: print 'Band %d'%b
        a_arr=a_init[b,:]
        ameas=zeros(2)
        ameas[0]=measbeam(radarr,beam_sm[:,b],nuarr,rsrfarr[:,b],nuc[b]*a_arr[0],alpha_nep[b],ind=ind)
        ameas[1]=measbeam(radarr,beam_sm[:,b],nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b],ind=ind)
        adiff=ameas-area[b]
        arel=(ameas-area[b])/area[b]
        if verbose:
            print 'a=%.4f, A=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[0],ameas[0],adiff[0],arel[0])
            print 'a=%.4f, A=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
        it=1
        done=False
        while done==False:
            grad=(adiff[1]-adiff[0])/(a_arr[1]-a_arr[0])
            a_arrnew=a_arr[1] - adiff[1]/grad
            a_arr[0]=a_arr[1]
            ameas[0]=ameas[1]
            a_arr[1]=a_arrnew
            ameas[1]=measbeam(radarr,beam_sm[:,b],nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b],ind=ind)
            adiff=ameas-area[b]        
            arel=(ameas-area[b])/area[b]
            if verbose:            
                print 'a=%.4f, Area=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
            it=it+1
            if abs(arel[1]) < aprec:
                done=True
            elif it >= maxit:
                done=True
        a_fin[b]=a_arr[1]
        ameas_fin[b]=ameas[1]
        arel_fin[b]=arel[1]
        if verbose:        
            print 'Final a=%.4f, Area=%.2f (Rel. error=%.6f)' % (a_fin[b],ameas_fin[b],arel_fin[b])
            print 'Nu(eff)=%.2f GHz'%(nuc[b]*a_fin[b]/1.e9)

    if verbose:
        print 'Meas areas: [%.2f, %.2f, %.2f]' % (ameas_fin[0],ameas_fin[1],ameas_fin[2])
        print 'True areas: [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])

    nueff=nuc*a_fin
    print 'Effective frequencies: [%.2f, %.2f, %.2f] GHz' % (nueff[0]/1.e9,nueff[1]/1.e9,nueff[2]/1.e9)
    return(nueff)
    
def get_effnu_az_newbeam(radarr,beam_scl,beam_fix,rsrfarr,nuarr,nuc,brad,alpha_nep,aprec=0.0001,verbose=None,ind=1.0):
    
    ## Find the effective frequency of the beam, based on measurements of Neptune
    ## Takes constant and core sections of beam
    
    from numpy import array,zeros,arange    
    from beam import beamarea_az
    #from beam import measbeam_new_wrong as measbeam_new
    #print "**USING INCORRECT BEAM SCALING**"
    from beam import measbeam_new
    band=arange(3)
    beam_cmb=zeros((brad,3))
    area=zeros(3)
    for b in band:
        beam_cmb=comb_beam(beam_scl,beam_fix)
        area[b]=beamarea_az(radarr,beam_cmb[:,b],brad=brad)

    if verbose: print 'Beam Areas (brad): [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])

    #Spectral index of neptune
    #alpha_nep=array((1.26,1.39,1.45))

    a_init=array([[1.,1.01],[1.,1.02],[1.,1.03]])

    maxit=10.
    radarr=arange(brad)
    a_fin=zeros(3)
    ameas_fin=zeros(3)
    arel_fin=zeros(3)

    for b in band:
        if verbose: print 'Band %d'%b
        a_arr=a_init[b,:]
        ameas=zeros(2)
        ameas[0]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nuc[b]*a_arr[0],alpha_nep[b],ind=ind)
        ameas[1]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b],ind=ind)
        adiff=ameas-area[b]
        arel=(ameas-area[b])/area[b]
        if verbose:
            print 'a=%.4f, Area=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[0],ameas[0],adiff[0],arel[0])
            print 'a=%.4f, Area=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
        it=1
        done=False
        while done==False:
            grad=(adiff[1]-adiff[0])/(a_arr[1]-a_arr[0])
            a_arrnew=a_arr[1] - adiff[1]/grad
            a_arr[0]=a_arr[1]
            ameas[0]=ameas[1]
            a_arr[1]=a_arrnew
            ameas[1]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nuc[b]*a_arr[1],alpha_nep[b],ind=ind)
            adiff=ameas-area[b]        
            arel=(ameas-area[b])/area[b]
            if verbose:            
                print 'a=%.4f, Area=%.2f (Diff=%.3f ; Rel=%.6f)' % (a_arr[1],ameas[1],adiff[1],arel[1])
            it=it+1
            if abs(arel[1]) < aprec:
                done=True
            elif it >= maxit:
                done=True
        a_fin[b]=a_arr[1]
        ameas_fin[b]=ameas[1]
        arel_fin[b]=arel[1]
        if verbose:        
            print 'Final a=%.4f, Area=%.2f (Rel. error=%.6f)' % (a_fin[b],ameas_fin[b],arel_fin[b])
            print 'Nu(eff)=%.2f GHz'%(nuc[b]*a_fin[b]/1.e9)

    if verbose:
        print 'Meas areas: [%.2f, %.2f, %.2f]' % (ameas_fin[0],ameas_fin[1],ameas_fin[2])
        print 'True areas: [%.2f, %.2f, %.2f]' % (area[0],area[1],area[2])

    nueff=nuc*a_fin
    print 'Effective frequencies: [%.2f, %.2f, %.2f] GHz' % (nueff[0]/1.e9,nueff[1]/1.e9,nueff[2]/1.e9)
    return(nueff)

def beamarea_nu(nuarr,beamx,beamy,beams,nueff,ilim,dnu_a=30,brad=None,verbose=None,ind=1.0):
    
    ## Calculate the beam area for a range of frequencies
    ## Uses a 2-d gridded beam
    ## Assumes fully extended source
    from numpy import arange,zeros
    from scipy.interpolate import RectBivariateSpline,interp1d
    
    band=arange(3)
    nnu=nuarr.size
    area_a=zeros((nnu,3))
    xlist=beamx[:,0]
    ylist=beamy[0,:]

    for b in band:
        nu_a=nuarr[ilim[b,0]:ilim[b,1]+dnu_a:dnu_a]
        nnu_a=nu_a.size
        area_bm=zeros(nnu_a)
        if verbose: print 'Band %d...' % b
        bm_sp=RectBivariateSpline(xlist,ylist,beams[:,:,b])
        for n in arange(nnu_a):
            regrid_x=xlist*(nu_a[n]/nueff[b])**ind
            regrid_y=ylist*(nu_a[n]/nueff[b])**ind
            regrid_bm=bm_sp(regrid_x,regrid_y)
            area_bm[n]=beamarea(beamx,beamy,regrid_bm,brad=brad)
            
        area_sp=interp1d(nu_a,area_bm)
        area_a[ilim[b,0]:ilim[b,1],b]=area_sp(nuarr[ilim[b,0]:ilim[b,1]])
    
    return(area_a)    

def beamarea_th_nu(nuarr,thetarr,beamx,beamy,beamin,nueff,ilim,dnu_a=30,brad=None,verbose=None,ind=1.0):
    
    ## Calculate the beam area for a range of frequencies and source widths
    ## Uses the 2-d gridded beam
    
    from numpy import arange,zeros
    from scipy.interpolate import RectBivariateSpline,interp1d
    
    band=arange(3)
    nth=thetarr.size
    nnu=nuarr.size
    area_a=zeros((nnu,nth,3))
    xlist=beamx[:,0]
    ylist=beamy[0,:]
    for b in band:
        nu_a=nuarr[ilim[b,0]:ilim[b,1]+dnu_a:dnu_a]
        nnu_a=nu_a.size
        area_bm=zeros(nnu_a)
        if verbose: print 'Band %d...' % b
        bm_sp=RectBivariateSpline(xlist,ylist,beamin[:,:,b])
        for t in arange(nth):
            theta=thetarr[t]
            if verbose: print 'theta %f...'%theta
            if theta < 1e-2:
                beamx1 = beamx * 1.e-3
                beamy1 = beamy * 1.e-3
                xlist1 = xlist * 1.e-3
                ylist1 = ylist * 1.e-3
            elif theta < 0.1:
                beamx1 = beamx * 1.e-2
                beamy1 = beamy * 1.e-2
                xlist1 = xlist * 1.e-2
                ylist1 = ylist * 1.e-2
            elif theta < 1.:
                beamx1 = beamx * 1.e-1
                beamy1 = beamy * 1.e-1
                xlist1 = xlist * 1.e-1
                ylist1 = ylist * 1.e-1
            else:
                beamx1 = beamx
                beamy1 = beamy
                xlist1 = xlist
                ylist1 = ylist
            
            for n in arange(nnu_a):
                regrid_x=xlist1*(nu_a[n]/nueff[b])**ind
                regrid_y=ylist1*(nu_a[n]/nueff[b])**ind
                regrid_bm=bm_sp(regrid_x,regrid_y)
                            
                area_bm[n]=beamarea_th(beamx1,beamy1,regrid_bm,theta)
            
            area_sp=interp1d(nu_a,area_bm)
            area_a[ilim[b,0]:ilim[b,1],t,b]=area_sp(nuarr[ilim[b,0]:ilim[b,1]])
    
    return(area_a)
    
def beamarea_az_th_nu(nuarr,thetarr,radarr,beam_sm,nueff,ilim,dnu_a=30,brad=None,verbose=None,ind=1.0):

    ## Calculate the beam area for a range of frequencies and source widths
    ## Uses azimuthally-smoothed beam (single file)
    
    from numpy import arange,zeros
    from scipy.interpolate import interp1d
    from scipy import where
    
    band=arange(3)
    nth=thetarr.size
    nnu=nuarr.size
    area_a=zeros((nnu,nth,3))

    for b in band:
        nu_a=nuarr[ilim[b,0]:ilim[b,1]+dnu_a:dnu_a]
        nnu_a=nu_a.size
        area_bm=zeros(nnu_a)
        if verbose: print 'Band %d...' % b
        bm_sp=interp1d(radarr,beam_sm[:,b])
        for t in arange(nth):
            theta=thetarr[t]
#            #Adjust max radius
#            if theta < 1e-2:
#                radarr1 = radarr*1.e-3
#            elif theta < 0.1:
#                radarr1 = radarr*1.e-2
#            elif theta < 1.:
#                radarr1 = radarr*1.e-1
#            else:
#                radarr1 = radarr

            #Keep max radius constant
            maxrad = theta*5 #set max radius
            if maxrad > max(radarr): maxrad=max(radarr) #limit to max of input
            radarr1=arange(0,maxrad,maxrad/radarr.size) # max new radius array
            
            for n in arange(nnu_a):
                regrid_rad=radarr1*(nu_a[n]/nueff[b])**ind #regrid to new frequency values
                #print 'nu: %.2f GHz'%(nu_a[n]*1.e-9)
                #print 'rad (regrid): [%.2f : %.2f]'%(min(regrid_rad),max(regrid_rad))
                ingrid=where(regrid_rad < max(radarr)) #check which ones are within beammap
                regrid_bm=zeros(radarr1.size) #make new beam array
                regrid_bm[ingrid]=bm_sp(regrid_rad[ingrid]) #fill array with new beam area
                            
                area_bm[n]=beamarea_az_th(radarr1,regrid_bm,theta) #calculate area
            
            area_sp=interp1d(nu_a,area_bm)
            area_a[ilim[b,0]:ilim[b,1],t,b]=area_sp(nuarr[ilim[b,0]:ilim[b,1]])
    
    return(area_a)    

def beamarea_az_th_nu_new(nuarr,thetarr,radarr,beam_scl,beam_fix,nueff,ilim,dnu_a=30,brad=None,verbose=None,ind=1.0):

    ## Calculate the beam area for a range of frequencies and source widths
    ## Uses the azimuthally-smoothed beam
    ## Uses core and constant beam profiles
    
    from numpy import arange,zeros
    from scipy.interpolate import interp1d
    from scipy import where
    
    band=arange(3)
    nth=thetarr.size
    nnu=nuarr.size
    area_a=zeros((nnu,nth,3))

    for b in band:
        nu_a=nuarr[ilim[b,0]:ilim[b,1]+dnu_a:dnu_a]
        nnu_a=nu_a.size
        area_bm=zeros(nnu_a)
        if verbose: print 'Band %d...' % b
        bm_scl_sp=interp1d(radarr,beam_scl[:,b])
        bm_fix_sp=interp1d(radarr,beam_fix[:,b])
        for t in arange(nth):
            theta=thetarr[t]
            ###theta is in FWHM
            
#            #Adjust max radius
#            if theta < 1e-2:
#                radarr1 = radarr*1.e-3
#            elif theta < 0.1:
#                radarr1 = radarr*1.e-2
#            elif theta < 1.:
#                radarr1 = radarr*1.e-1
#            else:
#                radarr1 = radarr

            #set max radius
            maxrad = theta*5
            if maxrad > max(radarr): maxrad=max(radarr) #limit to max of input
            radarr1=arange(0,maxrad,maxrad/radarr.size) # max new radius array
            
            for n in arange(nnu_a):
                regrid_rad=radarr1*(nu_a[n]/nueff[b])**ind #regrid to new frequency values
                #print 'nu: %.2f GHz'%(nu_a[n]*1.e-9)
                #print 'rad (regrid): [%.2f : %.2f]'%(min(regrid_rad),max(regrid_rad))
                ingrid=where(regrid_rad < max(radarr)) #check which ones are within beammap
                #make new beam arrays
                regrid_bm_scl=zeros(radarr1.size) 
                regrid_bm_fix=zeros(radarr1.size)
                regrid_bm_scl[ingrid]=bm_scl_sp(regrid_rad[ingrid]) #fill array with new beam profile
                regrid_bm_fix[ingrid]=bm_fix_sp(regrid_rad[ingrid])
                regrid_bm_cmb=comb_beam(regrid_bm_scl,regrid_bm_fix)
                
                area_bm[n]=beamarea_az_th(radarr1,regrid_bm_cmb,theta) #calculate area
            
            area_sp=interp1d(nu_a,area_bm)
            area_a[ilim[b,0]:ilim[b,1],t,b]=area_sp(nuarr[ilim[b,0]:ilim[b,1]])
    
    return(area_a)

def beamarea_az_nu(nuarr,radarr,beam_sm,nueff,ilim,dnu_a=30,brad=None,verbose=None,ind=1.0):

    ## Calculate the beam area for a range of frequencies over fully extended source
    ## Uses the azimuthally-smoothed beam (single file)
    
    from numpy import zeros
    from scipy.interpolate import interp1d
    from scipy import where
    
    band=range(3)
    nnu=nuarr.size
    area_a=zeros((nnu,3))

    for b in band:
        nu_a=nuarr[ilim[b,0]:ilim[b,1]+dnu_a:dnu_a]
        nnu_a=nu_a.size
        area_bm=zeros(nnu_a)
        if verbose: print 'Band %d...' % b
        bm_sp=interp1d(radarr,beam_sm[:,b])
        for n in range(nnu_a):
            regrid_rad=radarr*(nu_a[n]/nueff[b])**ind #regrid to new frequency values
                #print 'nu: %.2f GHz'%(nu_a[n]*1.e-9)
                #print 'rad (regrid): [%.2f : %.2f]'%(min(regrid_rad),max(regrid_rad))
            ingrid=where(regrid_rad < max(radarr)) #check which ones are within beammap
            regrid_bm=zeros(radarr.size) #make new beam array
            regrid_bm[ingrid]=bm_sp(regrid_rad[ingrid]) #fill array with new beam area
                            
            area_bm[n]=beamarea_az(radarr,regrid_bm) #calculate area
            
        area_sp=interp1d(nu_a,area_bm)
        area_a[ilim[b,0]:ilim[b,1],b]=area_sp(nuarr[ilim[b,0]:ilim[b,1]])
    
    return(area_a)    
        
def beamarea_az_nu_new(nuarr,radarr,beam_scl,beam_fix,nueff,ilim,dnu_a=30,brad=None,verbose=None,ind=1.0):

    ## Calculate the beam area for a range of frequencies over fully extended source
    ## Uses the azimuthally-smoothed beam
    ## Uses core and constant profiles
    
    from numpy import zeros
    from scipy.interpolate import interp1d
    from scipy import where
    
    band=range(3)
    nnu=nuarr.size
    area_a=zeros((nnu,3))

    for b in band:
        nu_a=nuarr[ilim[b,0]:ilim[b,1]+dnu_a:dnu_a]
        nnu_a=nu_a.size
        area_bm=zeros(nnu_a)
        if verbose: print 'Band %d...' % b
        bm_scl_sp=interp1d(radarr,beam_scl[:,b])
        bm_fix_sp=interp1d(radarr,beam_fix[:,b])
        for n in range(nnu_a):
            regrid_rad=radarr*(nu_a[n]/nueff[b])**ind #regrid to new frequency values
                #print 'nu: %.2f GHz'%(nu_a[n]*1.e-9)
                #print 'rad (regrid): [%.2f : %.2f]'%(min(regrid_rad),max(regrid_rad))
            ingrid=where(regrid_rad < max(radarr)) #check which ones are within beammap
            #make new beam array
            regrid_bm_scl=zeros(radarr.size)
            regrid_bm_fix=zeros(radarr.size)
            #fill arrays with new beam profile
            regrid_bm_scl[ingrid]=bm_scl_sp(regrid_rad[ingrid])
            regrid_bm_fix[ingrid]=bm_fix_sp(regrid_rad[ingrid])
            regrid_bm_cmb=comb_beam(regrid_bm_scl,regrid_bm_fix)
            area_bm[n]=beamarea_az(radarr,regrid_bm_cmb) #calculate area
            
        area_sp=interp1d(nu_a,area_bm)
        area_a[ilim[b,0]:ilim[b,1],b]=area_sp(nuarr[ilim[b,0]:ilim[b,1]])
    
    return(area_a)

def arcsec2sr(arcsec):
    
    ## Convert from arcsec^2 to steradians
    
    from numpy import pi
    
    conv=(60.*60.*180/pi)**2
    area_sr = arcsec / conv
    
    return(area_sr)

def sr2arcsec(area_sr):
    
    ## Convert from arcsec^2 to steradians
    
    from numpy import pi
    
    conv=1./(60.*60.*180/pi)**2
    area_arcsec = area_sr / conv
    
    return(area_arcsec)
    
def avgbeam(nuarr,area_bm):

    ## Calculate the mean beam area and frequency

    from scipy.integrate import trapz
    from scipy.interpolate import interp1d
    
    avgnu=trapz(nuarr*area_bm,nuarr)/trapz(area_bm,nuarr)
    intbm=interp1d(nuarr,area_bm)
    avgbeam=intbm(avgnu)
    
    return(avgnu,avgbeam)

def beamarea_eff(alpharr,nu0,nuarr,area_bm,rsrfarr):

    ## Calculate the effective beam area over the band for a range of spectral indices
    ## Requires beam area to already be calculated
    
    from numpy import arange,zeros    
    from scipy.integrate import trapz

    nalph=alpharr.size
    effarea=zeros(nalph)
    for a in arange(nalph):
        numint=trapz(area_bm[:] * nuarr**alpharr[a] * rsrfarr[:],nuarr)
        denomint=trapz(nu0**alpharr[a] * rsrfarr[:],nuarr)
        effarea[a]=numint/denomint
    
    return(effarea)

def modbeam(radarr,beam_sm,nuarr,rsrfarr,nu0,alpha_nep,ind=1.0):

    ## Calculate modelled beam profile, integrating over band with Neptune spectral index
    ## Uses single beam profile
    
    from numpy import zeros,arange
    from scipy import where
    from scipy.interpolate import interp1d
    from scipy.integrate import trapz
    
    nrad=radarr.size
    nnu=nuarr.size
    beammod=zeros(nrad)
    beam_sm_nu=zeros(nnu)
    
    beam_sm_int=interp1d(radarr,beam_sm)
    denom_arg = rsrfarr * nuarr**alpha_nep
    for r in arange(nrad):
        r_nu=radarr[r]*(nuarr/nu0)**ind
        inr=where(r_nu < max(radarr))
        outr=where(r_nu >= max(radarr))
        beam_sm_nu[inr]=beam_sm_int(r_nu[inr])
        beam_sm_nu[outr]=0.
        num_arg = rsrfarr * nuarr**alpha_nep * beam_sm_nu
        beammod[r]=trapz(num_arg,nuarr)/trapz(denom_arg,nuarr)
    
    return(beammod)

def modbeam_new(radarr,beam_scl,beam_fix,nuarr,rsrfarr,nu0,alpha_nep,ind=1.0):

    ## Calculate modelled beam profile over Neptune
    ## integrating over band with Neptune spectral index
    ## Uses core and constant profiles
    
    from numpy import zeros,arange
    from scipy import where
    from scipy.interpolate import interp1d
    from scipy.integrate import trapz
    
    #print radarr[0:10]
    #print beam_scl[0:10]
    #print beam_fix[0:10]
    #print nuarr[0:10]
    #print rsrfarr[0:10]
    #print nu0
    #print alpha_nep
    #print ind
    
    nrad=radarr.size
    nnu=nuarr.size
    beammod=zeros(nrad)
    beam_scl_nu=zeros(nnu)
    beam_fix_nu=zeros(nnu)
    
    beam_scl_int=interp1d(radarr,beam_scl)
    beam_fix_int=interp1d(radarr,beam_fix)
    denom_arg = rsrfarr * nuarr**alpha_nep
    for r in arange(nrad):
        r_nu=radarr[r]*(nuarr/nu0)**ind
        inr=where(r_nu < max(radarr))
        outr=where(r_nu >= max(radarr))
        beam_scl_nu[inr]=beam_scl_int(r_nu[inr])
        beam_scl_nu[outr]=0.
        beam_fix_nu[inr]=beam_fix_int(r_nu[inr])
        beam_fix_nu[outr]=0.
        beam_cmb_nu=where(beam_scl_nu < beam_fix_nu,beam_fix_nu,beam_scl_nu)
        num_arg = rsrfarr * nuarr**alpha_nep * beam_cmb_nu
        beammod[r]=trapz(num_arg,nuarr)/trapz(denom_arg,nuarr)
    
    #print beammod[0:10]
    return(beammod)

def beamprof_nu(radarr,beam_scl,beam_fix,nuarr,nu0,ind=1.0):

    ## Calculate modelled beam profile at a range of frequencies
    ## integrating over band
    
    from numpy import zeros,arange
    from scipy import where
    from scipy.interpolate import interp1d
    
    nrad=radarr.size
    nnu=nuarr.size
    beam_nu=zeros((nrad,nnu))
    beam_scl_nu=zeros(nnu)
    beam_fix_nu=zeros(nnu)
    
    beam_scl_int=interp1d(radarr,beam_scl)
    for r in arange(nrad):
        r_nu=radarr[r]*(nuarr/nu0)**ind
        inr=where(r_nu < max(radarr))
        outr=where(r_nu >= max(radarr))
        beam_scl_nu[inr]=beam_scl_int(r_nu[inr])
        beam_scl_nu[outr]=1.e-8
        beam_fix_nu[:]=beam_fix[r]
        beam_nu[r,:]=where(beam_scl_nu < beam_fix_nu,beam_fix_nu,beam_scl_nu)
        
    return(beam_nu)


def modbeam_area(radarr,beam_sm,area_rarr):

    ## Calculate area of az-smoothed beam over range of radii
    ## Uses single beam profile
    
    from scipy import where
    from scipy.integrate import trapz
    from numpy import zeros,arange,pi,array
    
    nrad=area_rarr.size
    area_mod=zeros(nrad)
    for r in arange(nrad):
        inr=where(radarr <= area_rarr[r])
        if array(inr).size > 0:
            area_mod[r]=trapz(2.*pi*radarr[inr]*beam_sm[inr],radarr[inr])
        #print r, area_rarr[r],array(inr).size,area_mod[r]
    return(area_mod)

#def rotbeam(beamx,beamy,ang_deg):
#
#    ## Rotate beam by angle
#    
#    from numpy import sin,cos,pi
#    #from scipy.interpolate import RectBivariateSpline
#    
#    ang=ang_deg*pi/180.
#    xnew=beamx*cos(ang)-beamy*sin(ang)
#    ynew=beamx*sin(ang)+beamy*cos(ang)
#    
#    return(xnew,ynew)

def srcarea(fwhm):

    ##calculate source area for given Gaussian FWHM

    from scipy.integrate import trapz
    from numpy import arange,exp,pi,sqrt,log
    
    ntheta=100
    theta0=fwhm/(2.*sqrt(log(2.)))
    thetamax=5.*theta0
    thetalist=arange(0,thetamax,thetamax/ntheta)
    area=trapz(exp(-(thetalist/theta0)**2)*2.*pi*thetalist,thetalist)
    
    return(area)

def getbeamareauni(fileBeamArea):
    
    from csv import reader
    from numpy import array
    import string
    from scipy.interpolate import interp1d
    
    fileBeamArea='../Inputs/freq_beam-area.csv'
    fBeam=reader(open(fileBeamArea))
    nuArr=[]
    beamArea=[]
    for row in fBeam:
        if string.find(row[0],'#') < 0:
            nuArr.append(float(row[0]))
            beamArea.append(float(row[1]))
    nuArr=array(nuArr)
    beamArea=array(beamArea)
    
    beamAreaInt=interp1d(nuArr,beamArea)
    
    return(beamAreaInt)

def beam2prof(radIntArr,beamMap,radList,mapMask=None):
    
    ## Turn beam map into beam profile
    
    from numpy import zeros,zeros_like,ones_like,average,size
    from scipy import where
    nRad=len(radList)
    mapOnes=ones_like(beamMap)
    beamProf=zeros_like(radList)
    weights=ones_like(radList)
    
    if mapMask == None:
        doWeight=False
    else:
        doWeight=True
    
    for r in range(nRad):
        beamProf[r]=average(beamMap[where(radIntArr == radList[r])])
        if doWeight:
            weights[r]=sum(mapOnes[where(radIntArr == radList[r])])/ \
                sum(mapMask[where(radIntArr == radList[r])])
    
    beamProf = beamProf / weights
    
    return(beamProf)

def beam2prof2(radArr,beamMap,radList,mapMask=None,retAll=False,verbose=False):
    
    ## Turn beam map into beam profile
    
    from numpy import zeros,zeros_like,ones_like,average,size
    from scipy import where

    if mapMask == None:
        doWeight=False
    else:
        doWeight=True
    
    nRad=len(radList)
    dR=radList[2]-radList[1]
    mapOnes=ones_like(beamMap)
    beamProf=zeros_like(radList)
    beamArea=zeros_like(radList)
    nPixTot=zeros_like(radList)
    nValTot=zeros_like(radList)
    nPix=zeros_like(radList)
    nVal=zeros_like(radList)
    
    for r in range(nRad):
        ##Calculate beam area inside radius
        beamArea[r]=sum(beamMap[where(radArr <= (radList[r]+dR))])
        ##Calcualte total number of pixels inside radius
        nPixTot[r]=sum(mapOnes[where(radArr <= (radList[r]+dR))])
        ##Calculate number of valid pixels inside radius
        if doWeight:
            nValTot[r]=sum(mapMask[where(radArr <= (radList[r]+dR))])
        else:
            nValTot[r]=nPixTot[r]
        
        if r==0:
            nPix[r]=nPixTot[r]
            nVal[r]=nValTot[r]
            beamProf[r]=beamArea[r]/nVal[r]
        else:
            nPix[r] = nPixTot[r] - nPixTot[r-1]
            nVal[r] = nValTot[r] - nValTot[r-1]
            if nVal[r] == 0:
                beamProf[r]=beamProf[r-1]
            else:
                beamProf[r]=(beamArea[r]-beamArea[r-1])/nVal[r]
        if verbose:
            print radList[r],beamArea[r],beamProf[r],nVal[r],nPix[r]
    if retAll:
        #return all parameters
        return(beamProf,beamArea,nPix,nVal)
    else:
        return(beamProf)
        
def beam_fwhm(fwhm0,nuarr,rsrfarr,nueff,alpha,ind):
    
    from scipy import integrate

    denom_arg = rsrfarr * nuarr**alpha
    num_arg = fwhm0*(nuarr/nueff)**ind * rsrfarr * nuarr**alpha
    beamfwhm=integrate.trapz(num_arg,nuarr)/integrate.trapz(denom_arg,nuarr)

    return(beamfwhm)
