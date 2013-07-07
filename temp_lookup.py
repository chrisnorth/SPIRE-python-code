from numpy import arange,array,zeros,log

def maincode():
    
    import matplotlib.pyplot as plot
    
    c=299792458.
    
    #set temperature range and step
    Tmin=1.
    Tmax=50.
    Tstep=0.1
    
    Tarr=arange(Tmin,Tmax+Tstep,Tstep)
    nT=Tarr.size

    Iarr=zeros(6)
    SIarr=zeros((nT,6))
    
    wllist=array((70,100,160,250.,350.,500.))*1.e-6    
    nulist=c/wllist
    wl0=100.e-6
    nu0=c/wl0
    beta=2. #frequency-dependent emmisivity power-law
    
    fileout='/data/Herschel/Calibration/Code/temp_lookup.dat'
    
    plot.figure(1)
    plot.cla()
    
    for t in range(nT):
        Temp=Tarr[t]

        #calculate Intensity at all three wavelengths
        for i in range(6):
            Iarr[i]=greybody(1.,nulist[i],Temp,nu0,beta)
        #Iarr=Iarr/Iarr[0] #normalise to 250-micron intensity

        for i in range(5):
            SIarr[t,i]=(log(Iarr[i])-log(Iarr[i+1]))/(log(nulist[i])-log(nulist[i+1]))
        
        plot.plot(nulist/1.e9,Iarr)
        plot.xlabel('Frequency [GHz]')
        plot.ylabel('Intensity')
        plot.yscale('log')
        #print Temp,Temp/5.,int(Temp/5.)
        if Temp <= 10.:
            if int(t/10)==t/10.:
                plot.annotate(str(Temp),(nulist[0]/1.e9,Iarr[0]),verticalalignment='center')
        else:
            if int(t/100)==t/100.:
                plot.annotate(str(Temp),(nulist[0]/1.e9,Iarr[0]),verticalalignment='center')
    
    SIarr[:,5]=SIarr[:,3]/SIarr[:,4]

    f=open(fileout,'w')
    f.write('#Temp (K), alpha_250-350, alpha_350-500, alpha_250-350/alpha_350-500\n')
    
    for t in range(nT):
        f.write('%2.1f, %.4f, %.4f, %.4f\n'%(Tarr[t],SIarr[t,3],SIarr[t,4],SIarr[t,5]))
    
    f.close()
    
    print nulist
    print Tarr
    print SIarr
    
    plot.figure(2)
    plot.subplot(2,1,1)
    plot.cla()
    plot.plot(Tarr,SIarr[:,3],label=r'$\alpha_{250-350}$')
    plot.plot(Tarr,SIarr[:,4],label=r'$\alpha_{350-500}$')
    plot.ylabel(r'Spectral Index, $\alpha$')
    plot.axhline(0,color='k',linestyle='--')
    plot.legend(loc='lower right',ncol=2)
    plot.subplot(2,1,2)
    plot.cla()
    plot.axhline(1,color='k',linestyle=':')    
    plot.plot(Tarr,SIarr[:,5])
    plot.xlabel('Temperature [K]')
    plot.ylabel(r'$\alpha_{250-350} / \alpha_{350-500}$')
    plot.ylim(-5,5)
    plot.axhline(0,color='k',linestyle='--')
    
    plot.show()

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

if __name__=="__main__":
    maincode()