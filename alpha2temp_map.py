from numpy import zeros,ones,isnan,nan,array,arange,log,min,max,sqrt
from scipy import where
import pyfits
from csv import reader
from sys import exit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plot

def maincode():
    
    rel_err=0.03
    beta=2.

    path0='/data/Herschel/Calibration/'
    
    #Fit temperatures for proposed scheme
    pathin=path0+'Maps/Intensity maps/Proposed Photometry Scheme'
    pathout=path0+'Maps/Dust maps/Proposed Photometry Scheme'
    Tempr_prop,Tempr_prop_err = gettemp(path0,pathin,pathout,beta,rel_err)

    #Fit temperatures for current scheme
    pathin=path0+'Maps/Intensity maps/Current Photometry Scheme'
    pathout=path0+'Maps/Dust maps/Current Photometry Scheme'
    Tempr_curr,Tempr_curr_err = gettemp(path0,pathin,pathout,beta,rel_err)

    #get basic FITS header
    filehdr=pathin+'/PSW.fits'
    pathout=path0+'Maps/Dust maps/Current vs Proposed'
    filesout=[pathout+'/Temp_ar_map.fits',pathout+'/Temp-err_ar_map.fits']
    
    hdr0_in=pyfits.getheader(filehdr)
    hdr1_in=pyfits.getheader(filehdr,ext=1)
    
    dT = Tempr_curr - Tempr_prop
    dT_err = sqrt(Tempr_curr_err**2 + Tempr_prop_err**2)
    
    prihdu=pyfits.PrimaryHDU()
    hdrdt=hdr1_in.copy()
    hdrdt.update('UNITS','K',comment='Units of map')
    hdrdt.update('DATATYPE','T_grey_fit Difference',comment='Best-fit greybody-temperature difference')
    hdrdt.update('FITPARAM','Alpha_250-300/Alpha_350-500','Parameter used for fitting')
    hdrdt.update('CALIB','Current-Proposed','Current - Proposed photomoetry calibration scheme')
    hdrdt.update('BETA',beta,comment='Emmisivity spectral index assumed')
    hdrdt.update('REL_ERR',rel_err,comment='Relative error on input data')
    hdudt=pyfits.ImageHDU(dT,header=hdrdt)
    hdulist=pyfits.HDUList([prihdu,hdudt])
    hdulist.writeto(filesout[0],clobber=True)
    
    prihdu=pyfits.PrimaryHDU()
    hdrdte=hdrdt.copy()
    hdrdte.update('DATATYPE','T_grey_fit Difference Error',comment='Error on best-fit greybody-temperature difference')
    hdudte=pyfits.ImageHDU(dT,header=hdrdte)
    hdulist=pyfits.HDUList([prihdu,hdudte])
    hdulist.writeto(filesout[1],clobber=True)
    
def gettemp(path0,pathin,pathout,beta,rel_err):

    c=299792458. #speed of light
    
    filesin=[pathin+'/PSW.fits',pathin+'/PMW.fits',pathin+'/PLW.fits']
    filesout=[pathout+'/Temp_a1_map.fits',pathout+'/Temp_a2_map.fits',pathout+'/Temp_ar_map.fits',
              pathout+'/Temp-err_a1_map.fits',pathout+'/Temp-err_a2_map.fits',pathout+'/Temp-err_ar_map.fits']
    
    print 'Reading in maps...'    
    map_psw=pyfits.getdata(filesin[0])
    map_pmw=pyfits.getdata(filesin[1])
    map_plw=pyfits.getdata(filesin[2])
    
    hdr0_psw=pyfits.getheader(filesin[0])
    hdr0_pmw=pyfits.getheader(filesin[1])
    hdr0_plw=pyfits.getheader(filesin[2])
    
    hdr1_psw=pyfits.getheader(filesin[0],ext=1)
    hdr1_pmw=pyfits.getheader(filesin[1],ext=1)
    hdr1_plw=pyfits.getheader(filesin[2],ext=1)
    
    mapshape=map_psw.shape
    
    nx=mapshape[0]
    ny=mapshape[1]
    print nx,ny
    obs_in=zeros((nx,ny,3))
    obs_in[:,:,0]=map_psw
    obs_in[:,:,1]=map_pmw
    obs_in[:,:,2]=map_plw
    
    #obs_in=obs_in[400:500,400:500,:]

    print 'Making mask...'    
    obs_in=where(obs_in<0,nan,obs_in)
    obs_in=where(obs_in>1.e10,nan,obs_in)

    mask=ones((nx,ny))
    mask=where(isnan(obs_in[:,:,0]),0.,1.)
    mask=where(isnan(obs_in[:,:,1]),0.,mask)
    mask=where(isnan(obs_in[:,:,2]),0.,mask)
    
    afile=path0+'Code/temp_lookup.dat'
    ainput=reader(open(afile))
    a1in=[]
    a2in=[]
    arin_l1=[-999.]
    arin_g1=[]
    Tin=[]
    Tin_l1=[999.]
    Tin_g1=[]
    for row in ainput:
        if row[0][0]!='#':
            Tin.append(float(row[0]))
            a1in.append(float(row[1]))
            a2in.append(float(row[2]))
            if float(row[3]) <= 1:
                arin_l1.append(float(row[3]))
                Tin_l1.append(float(row[0]))
            else:
                arin_g1.append(float(row[3]))
                Tin_g1.append(float(row[0]))

    Tin_g1.append(min(Tin_l1))
    arin_g1.append(999.)
        
    Tin=array(Tin)
    a1in=array(a1in)
    a2in=array(a2in)
    Tin_l1=array(Tin_l1)
    arin_l1=array(arin_l1)
    Tin_l1[0]=max(Tin_g1)
    Tin_g1=array(Tin_g1)
    arin_g1=array(arin_g1)

    print 'a1:',min(a1in),max(a1in)
    print 'a2:',min(a2in),max(a2in)
    print 'ar<1:',min(arin_l1),max(arin_l1)
    print 'ar>1:',min(arin_g1),max(arin_g1)
    
    a1_int=interp1d(a1in,Tin,bounds_error=False,fill_value=-1.)
    a2_int=interp1d(a2in,Tin,bounds_error=False,fill_value=-1.)
    ar_int_l1=interp1d(arin_l1,Tin_l1,bounds_error=False,fill_value=-1.)
    ar_int_g1=interp1d(arin_g1,Tin_g1,bounds_error=False,fill_value=-1.)
    
    a1p=arange(-10,3,0.1)
    a2p=arange(-6,3,0.1)
    arp=arange(-4,4,0.1)
    
    T1p=a1_int(a1p)
    T2p=a2_int(a2p)
    
    arp_l1=arp[where(arp<=1.)]
    arp_g1=arp[where(arp>=1.)]
    Trp=where(arp < 1.,ar_int_l1(arp),ar_int_g1(arp))
    Trp_l1=ar_int_l1(arp_l1)
    Trp_g1=ar_int_g1(arp_g1)
    
    
    print 'Calculating alpha...'
    wl=array([250.,350.,500.])*1.e-6
    nu=c/wl
    
    alpha1=log(map_psw/map_pmw)/log(nu[0]/nu[1])
    alpha2=log(map_pmw/map_plw)/log(nu[1]/nu[2])
    alphar=alpha1/alpha2
    
    print 'alpha 1:',min(alpha1),max(alpha1)
    print 'alpha 2:',min(alpha2),max(alpha2)
    print 'alpha r:',min(alphar),max(alphar)
    
    Temp1=a1_int(alpha1)
    Temp2=a2_int(alpha2)
    Tempr=where(alphar < 1.,ar_int_l1(alphar),ar_int_g1(alphar))
    
    alpha_err = 2.*rel_err    
    Temp1_lo=a1_int(alpha1-alpha_err)
    Temp1_hi=a1_int(alpha1+alpha_err)
    Temp2_lo=a2_int(alpha2-alpha_err)
    Temp2_hi=a2_int(alpha2+alpha_err)
    Tempr_lo=where(alphar < 1.,ar_int_l1(alphar-alpha_err),ar_int_g1(alphar-alpha_err))
    Tempr_hi=where(alphar < 1.,ar_int_l1(alphar+alpha_err),ar_int_g1(alphar+alpha_err))
    
    Temp1_err = (Temp1_hi - Temp1_lo)/2.
    Temp2_err = (Temp2_hi - Temp2_lo)/2.
    Tempr_err = (Tempr_hi - Tempr_lo)/2.
    
    Temp1=where(mask==1,Temp1,0)
    Temp2=where(mask==1,Temp2,0)
    Tempr=where(mask==1,Tempr,0)

    Temp1_err=where(mask==1,Temp1_err,nan)
    Temp2_err=where(mask==1,Temp2_err,nan)
    Tempr_err=where(mask==1,Tempr_err,nan)

    
    prihdu=pyfits.PrimaryHDU()
    hdrt1=hdr1_psw.copy()
    hdrt1.update('UNITS','K',comment='Units of map')
    hdrt1.update('DATATYPE','T_grey_fit',comment='Best-fit greybody-temperature')
    hdrt1.update('FITPARAM','Alpha_250-300','Parameter used for fitting')
    hdrt1.update('CALIB','Proposed','Proposed pohtomoetry calibration scheme')
    hdrt1.update('BETA',beta,comment='Emmisivity spectral index assumed')
    hdrt1.update('TEMP_MIN',min(Tin),comment='Min of temperature range in fit')
    hdrt1.update('TEMP_MAX',max(Tin),comment='Max of temperature range in fit')
    hdrt1.update('REL_ERR',rel_err,comment='Relative error on input data')
    hdut1=pyfits.ImageHDU(Temp1,header=hdrt1)
    hdulist=pyfits.HDUList([prihdu,hdut1])
    hdulist.writeto(filesout[0],clobber=True)

    hdrt2=hdrt1.copy()
    hdrt2.update('FITPARAM','Alpha_350-500','Parameter used for fitting')
    hdut2=pyfits.ImageHDU(Temp2,header=hdrt2)
    hdulist=pyfits.HDUList([prihdu,hdut2])
    hdulist.writeto(filesout[1],clobber=True)

    hdrtr=hdrt1.copy()
    hdrtr.update('FITPARAM','Alpha_250-350/Alpha_350-500','Parameter used for fitting')
    hdutr=pyfits.ImageHDU(Tempr,header=hdrtr)
    hdulist=pyfits.HDUList([prihdu,hdutr])
    hdulist.writeto(filesout[2],clobber=True)
    
    hdrt1_err=hdrt1.copy()
    hdrt1_err.update('DATATYPE','T_grey_fit Error',comment='Error on best-fit greybody-temperature')
    hdut1_err=pyfits.ImageHDU(Tempr_err,header=hdrt1_err)
    hdulist=pyfits.HDUList([prihdu,hdut1_err])
    hdulist.writeto(filesout[3],clobber=True)

    hdrt2_err=hdrt2.copy()
    hdrt2.update('DATATYPE','T_grey_fit Error',comment='Error on best-fit greybody-temperature')
    hdut2_err=pyfits.ImageHDU(Tempr_err,header=hdrt2_err)
    hdulist=pyfits.HDUList([prihdu,hdut2_err])
    hdulist.writeto(filesout[4],clobber=True)

    hdrtr_err=hdrtr.copy()
    hdrtr.update('DATATYPE','T_grey_fit Error',comment='Error on best-fit greybody-temperature')
    hdutr_err=pyfits.ImageHDU(Tempr_err,header=hdrtr_err)
    hdulist=pyfits.HDUList([prihdu,hdutr_err])
    hdulist.writeto(filesout[5],clobber=True)
    
    return(Tempr,Tempr_err)
    
    #plot.figure(1)
    #plot.clf()
    #plot.plot(a1p,T1p,label='a1')
    #plot.plot(a2p,T2p,label='a2')
    #plot.plot(arp,Trp,label='ar')
    ##plot.plot(arp_l1,Trp_l1,'--')
    ##plot.plot(arp_g1,Trp_g1,'--')
    #plot.axhline(0.,color='k')
    #plot.legend(loc='upper left')
   # 
   # plot.figure(2)
    #plot.clf()
    #plot.plot(Tin,a1in,label='a1')
    #plot.plot(Tin,a2in,label='a2')
    #plot.plot(Tin_l1,arin_l1,label='ar<1')
    #plot.plot(Tin_g1,arin_g1,label='ar>1')
    #plot.ylim(-50,50)
    #plot.legend(loc='upper right')
    
    #plot.show()
    
if __name__=="__main__":
    maincode()