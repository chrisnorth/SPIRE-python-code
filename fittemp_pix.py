from numpy import exp,array,zeros,isnan,nan,ones,log10,min,max
from scipy import where
import matplotlib.pyplot as plot
import pyfits
import sys

def maincode():

    print 'Reading maps...'
    path='/data/Herschel/Calibration/Maps/Intensity maps/Proposed Photometry Scheme'
    pathn='/data/Herschel/Calibration/Maps/TMC'
    files=[path+'/PSW.fits',path+'/PMW.fits',path+'/PLW.fits']

    map_psw=pyfits.getdata(files[0])
    map_pmw=pyfits.getdata(files[1])
    map_plw=pyfits.getdata(files[2])
    
    hdr0_psw=pyfits.getheader(files[0])
    hdr0_pmw=pyfits.getheader(files[1])
    hdr0_plw=pyfits.getheader(files[2])
    
    hdr1_psw=pyfits.getheader(files[0],ext=1)
    hdr1_pmw=pyfits.getheader(files[1],ext=1)
    hdr1_plw=pyfits.getheader(files[2],ext=1)
    
    mapshape=map_psw.shape
    
    nx=mapshape[0]
    ny=mapshape[1]
    print nx,ny
    obs_in=zeros((nx,ny,3))
    obs_in[:,:,0]=map_psw
    obs_in[:,:,1]=map_pmw
    obs_in[:,:,2]=map_plw
    
    px=410
    py=460
    obs_pix=obs_in[px,py,:]
    obs_in=obs_in[px-20:px+20,py-20:py+20,:]

    print 'obs_pix:',obs_pix
    obs_in=where(obs_in<0,nan,obs_in)
    obs_in=where(obs_in>1.e10,nan,obs_in)

    rel_err=0.05
    err_in=rel_err*obs_in
    err_pix=rel_err*obs_pix    
    
    mask=ones((nx,ny))
    mask=where(isnan(obs_in[:,:,0]),0.,1.)
    mask=where(isnan(obs_in[:,:,1]),0.,mask)
    mask=where(isnan(obs_in[:,:,2]),0.,mask)
    #print obs_in.shape
    
    ##read maps
    
    #set chisq parameters
    Trange=(3,40) #temperature range to scan
    beta=2. #emissivity spectral index
    Tprec=1.0 #precision on Temperature

    c= 299792458.
    wls=array((250.,350.,500.)) #in microns
    freqs=c/(wls * 1.e-6) #in Hz
    #obsvals=array((1.05,1.45,1.65))
    
    ##do chi-squared test
    print 'Doing chi-squared test...'
    t_map,a_map,chisq_map = dochisq_map(freqs,obs_in,err_in,beta,Trange=Trange,Tprec=Tprec)
    
    t_arr,a_arr,chisq_arr = dochisq_pix(freqs,obs_pix,err_pix,beta,Trange=Trange,Tprec=Tprec)
    #print T_fit,a_fit
    
        
    t_map=where(mask == 1,t_map,nan)
    a_map=where(mask == 1,a_map,nan)
    chisq_map=where(mask == 1,chisq_map,nan)

    #nulist=arange    

    print 'Plotting Temperature map...'
    plot.figure(1,figsize=(6,8))
    plot.clf()
    plot.contourf(t_map,32)
    plot.axis('equal')
    plot.colorbar()
    plot.title('Greybody temperature (K)')
    plot.plot(20,20,'wo')
    #plot.savefig('/data/Herschel/Calibration/Maps/Temp_map_proposed.png',dpi=300)
    plot.show()
    
    print 'Plotting Magnitude map...'    
    plot.figure(2,figsize=(6,8))
    plot.clf()
    plot.contourf(a_map,32)
    plot.axis('equal')
    plot.colorbar()
    plot.title('Greybody amplitude')
    plot.plot(20,20,'wo')
    #plot.savefig('/data/Herschel/Calibration/Maps/Amp_map_proposed.png',dpi=300)
    plot.show()
    
    print 'Plotting Chi-sq map...'    
    plot.figure(3,figsize=(6,8))
    plot.clf()
    plot.contourf(chisq_map,32)
    plot.axis('equal')
    plot.colorbar()
    plot.title('Chi-squared fit')
    plot.plot(20,20,'wo')
    #plot.savefig('/data/Herschel/Calibration/Maps/Chi_sq_map_proposed.png',dpi=300)
    plot.show()
        
    print 'Chi-sq:',chisq_arr
    print 'Mag:',a_arr
    print 'Temp:',t_arr
    print min(chisq_arr)
    print 'Plotting chi-sq arr...'
    plot.figure(4,figsize=(6,8))
    plot.clf()
    plot.contourf(a_arr,t_arr,log10(chisq_arr),64)
    plot.colorbar()
    plot.ylabel('Temperature (K)')
    plot.xlabel('Magnitude (arb. units.)')
    plot.xscale('log')
    plot.plot(a_map[20,20],t_map[20,20],'wo')
    #plot.savefig('/data/Herschel/Calibration/Maps/PSW_map_proposed.png',dpi=300)
    plot.show()
    
    print 'Inputs:',obs_pix
    print 'Temp:',t_map[20,20]    
    print 'Mag:',a_map[20,20]
    print 'Chi-sq:',chisq_map[20,20]

def dochisq_map(freqs,obs_in,err_in,beta,Trange,Tprec):

    from numpy import zeros,arange,log,nan,min,max
    from scipy import where
    
    chi_init=999.    
    
    nx=obs_in.shape[0]
    ny=obs_in.shape[1]
    
    print Trange
    Tfit=arange(Trange[0],Trange[1]+Tprec,Tprec)
    nt=Tfit.size    
    nu0=freqs[2]

    na=10
    
    #t_arr=zeros((nt,na))
    #a_arr=zeros((nt,na))
    chisq_map=zeros((nx,ny))
    chisq_x=zeros((nx,ny))
    chisq_map[:,:]=chi_init #set initial value
    t_map=zeros((nx,ny))
    a_map=zeros((nx,ny))
    for t in range(nt):
        Temp=Tfit[t]
        #t_arr[:,:,t,:]=Temp
        obsmin=min(obs_in[where(obs_in == obs_in)])
        obsmax=max(obs_in[where(obs_in == obs_in)])
        print 'Temperature: %.1f'%(Temp)
        logamin=log(0.5*obsmin/greybody(1.,nu0,Temp,nu0,beta=beta))
        logamax=log(2.*obsmax/greybody(1.,nu0,Temp,nu0,beta=beta))
        logarng=logamax-logamin
        logalst=arange(na)*logarng/na + (logamin)
        alist=exp(logalst)
        #print alist
        for a in range(na):
            #a_arr[:,:,t,a]=alist[a]
            chisq_x[:,:]=0.
            for b in range(3):
                pred_val=greybody(alist[a],freqs[b],Temp,nu0,beta)
                chisq_x=where(isnan(obs_in[:,:,b]),chi_init,
                            chisq_x + (obs_in[:,:,b]-pred_val)**2/err_in[:,:,b]**2)
            #print 'chisq_map: ',min(chisq_map),max(chisq_map)
            #print 'chisq_x: ',min(chisq_x),max(chisq_x)
            #print 'temp: ',min(t_map),max(t_map)
            t_map=where(chisq_x < chisq_map,Temp,t_map)
            a_map=where(chisq_x < chisq_map,alist[a],a_map)
            chisq_map=where(chisq_x < chisq_map,chisq_x,chisq_map)
            
        print '%.1f%%'%(100.*(t+1)/nt)
                
    return(t_map,a_map,chisq_map)

def dochisq_pix(freqs,obs,errors,beta,Trange=(5,15),Tprec=0.5):
    
    from numpy import arange,log,exp,zeros,log10,min
    from scipy import where

    #print 'Doing chi-squared test...'
    
    Tfit=arange(Trange[0],Trange[1]+Tprec,Tprec)
    nt=Tfit.size    
    nu0=freqs[2]

    na=20
    
    a_arr=zeros((nt,na))
    t_arr=zeros((nt,na))
    chisq_arr=zeros((nt,na))
    for t in range(nt):
        Temp=Tfit[t]
        t_arr[t,:]=Temp
        aini=obs[0]/greybody(1.,nu0,Temp,nu0,beta=beta)
        
        #print 'A_ini:',aini
        logaini=log(aini)
        logalst=arange(na)*4./na + (logaini-2.)
        alist=exp(logalst)
        #print alist
        a_arr[t,:]=alist
        for a in range(na):
            for b in range(3):
                pred_val=greybody(alist[a],freqs[b],Temp,nu0,beta)
                chisq_arr[t,a]=chisq_arr[t,a] + ((obs[b]-pred_val)**2)/errors[b]**2
    
    #print 'min chi-square:',min(chisq_arr)
    T_fit=t_arr[where(chisq_arr == min(chisq_arr))]
    a_fit=a_arr[where(chisq_arr == min(chisq_arr))]
    #print 'best-fit tempt:',t_arr[where(chisq_arr == min(chisq_arr))]
    
    return t_arr,a_arr,chisq_arr

def greybody(amult,nu,T,nu0,beta):
    #produce a blackbody modified by nu**beta power law
    # units are W/m^2/sr/Hz (1e-26 Jy/sr)
    b=planck(nu,T)

    i_nu = amult * b * (nu/nu0)**beta
    
    #i_nu = i_nu * 1.e-20 # downgrade units
    return(i_nu)

def planck(nu,T):

    from numpy import expm1
    #return the planck function for frequency and temperature
    # units are W/m^2/sr/Hz
    h=6.626e-34
    c=299792458.
    kB=1.38e-23
    
    nom=1.e26 * 2.*h*nu**3 / c**2 ##units are janskys
    denom=expm1(h*nu/(kB*T))
    #print nom,denom
    planck = nom/denom
    
    return(planck)

def dBdT(nu,T):
    
    from numpy import sinh

    ##return the derivative of the Planck function

    h=6.626e-34
    c=299792458.
    kB=1.38e-23
    
    x=h*nu/(kB*T)
    arg=x/sinh(x/2.)

    fac=kB*nu**2/(2.*c**2)
    
    deriv=fac * arg**2
    
    return(deriv)

#def dXdA(A)

if __name__=="__main__":
    maincode()
