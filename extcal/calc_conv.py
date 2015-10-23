def calc_k4(alpha,RSRFarr,apfarr,nuarr,dnu,nu0):
    
#    import matplotlib
#    import matplotlib.pyplot as plot    
#    import scipy
    import scipy.integrate as integrate
    
    num_arg=RSRFarr * apfarr * nu0**alpha
    denom_arg=RSRFarr * apfarr * nuarr**alpha
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
#    plot.figure(4)
#    plot.clf()
#    plot.plot(nuarr,num_arg,'r-')  
#    plot.plot(nuarr,denom_arg,'g-')
    
    K4 = num_int / denom_int
    #print 'K4= %f' % K4
    return(K4)

def calc_k4e(alpha,RSRFarr,apfarr,nuarr,area_nu0,ind,dnu,nu0):

##assumes area is power-law area variation, index -2*ind
    
#    import matplotlib
#    import matplotlib.pyplot as plot    
#    import scipy
    import scipy.integrate as integrate
    
    num_arg=RSRFarr * apfarr
    denom_arg=area_nu0 * RSRFarr * apfarr * (nuarr/nu0)**(-2*ind) * (nuarr/nu0)**alpha
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
#    plot.figure(4)
#    plot.clf()
#    plot.plot(nuarr,num_arg,'r-')  
#    plot.plot(nuarr,denom_arg,'g-')
    
    K4e = num_int / denom_int
    #print 'K4= %f' % K4
    return(K4e)

def calc_k4pip(alpha,RSRFarr,nuarr,dnu,nu0):

    ##uses the old system to calculate the conversion parameters
    ##i.e. no aperture efficiency
    
#    import matplotlib
#    import matplotlib.pyplot as plot    
#    import scipy
    import scipy.integrate as integrate
    
    num_arg=RSRFarr * nu0**alpha
    denom_arg=RSRFarr * nuarr**alpha
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
#    plot.figure(4)
#    plot.clf()
#    plot.plot(nuarr,num_arg,'r-')  
#    plot.plot(nuarr,denom_arg,'g-')
    
    K4pip = num_int / denom_int
    #print 'K4= %f' % K4
    return(K4pip)

def calc_kc(alpha,RSRFarr,apfarr,nuarr,dnu,nu0,alpha0):
    
#    import matplotlib
#    import matplotlib.pyplot as plot    
#    import scipy
    import scipy.integrate as integrate
    
    num_arg=RSRFarr * apfarr * nuarr**alpha0 * nu0**(alpha-alpha0)
    denom_arg=RSRFarr * apfarr * nuarr**alpha
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
#    plot.figure(4)
#    plot.clf()
#    plot.plot(nuarr,num_arg,'r-')  
#    plot.plot(nuarr,denom_arg,'g-')
    
    Kc = num_int / denom_int
    #print 'K4= %f' % K4
    return(Kc)

def calc_kcpip(alpha,RSRFarr,nuarr,dnu,nu0,alpha0):
    
#    import matplotlib
#    import matplotlib.pyplot as plot    
#    import scipy
    import scipy.integrate as integrate
    
    num_arg=RSRFarr * nuarr**alpha0 * nu0**(alpha-alpha0)
    denom_arg=RSRFarr * nuarr**alpha
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
#    plot.figure(4)
#    plot.clf()
#    plot.plot(nuarr,num_arg,'r-')  
#    plot.plot(nuarr,denom_arg,'g-')
    
    Kc = num_int / denom_int
    #print 'K4= %f' % K4
    return(Kc)


def calc_k5(alpha,nuarr,rsrfarr,apfarr,areaarr,nu0):

    import scipy.integrate as integrate

    num_arg=nu0**alpha * rsrfarr * apfarr
    denom_arg=nuarr**alpha * rsrfarr * areaarr * apfarr
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
    K5 = num_int / denom_int
    
    return(K5)

def calc_kc2(alpha,nuarr,rsrfarr,apfarr,areaarr,nu0,alpha0=-1.):
    
    import scipy.integrate as integrate
    
    num_arg = nu0**(alpha - alpha0) * rsrfarr * apfarr * nuarr**alpha
    denom_arg = rsrfarr * apfarr * areaarr * nuarr**alpha
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
    Kc2 = num_int / denom_int
    
    return(Kc2)

def greybody(nu,T,beta,nu0,amult=1.):
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

def calc_k4temp(temp,beta,RSRFarr,apfarr,nuarr,dnu,nu0):
    
#    import matplotlib
#    import matplotlib.pyplot as plot    
#    import scipy
    import scipy.integrate as integrate

    num_arg=RSRFarr * apfarr * greybody(nu0,temp,beta,nu0)
    denom_arg=RSRFarr * apfarr * greybody(nuarr,temp,beta,nu0)
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
#    plot.figure(4)
#    plot.clf()
#    plot.plot(nuarr,num_arg,'r-')  
#    plot.plot(nuarr,denom_arg,'g-')
    
    K4t = num_int / denom_int
    #print 'K4= %f' % K4
    return(K4t)

def calc_k5temp(temp,beta,nuarr,rsrfarr,apfarr,areaarr,nu0):

    import scipy.integrate as integrate

    num_arg=greybody(nu0,temp,beta,nu0) * rsrfarr * apfarr
    denom_arg=greybody(nuarr,temp,beta,nu0) * rsrfarr * areaarr * apfarr
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
    K5t = num_int / denom_int
    
    return(K5t)

def calc_k4uni(alpha,nuarr,apfarr,nu0):
    import scipy.integrate as integrate

    num_arg = nu0**alpha * apfarr
    denom_arg = nuarr**alpha * apfarr
    
    num_int = integrate.trapz(num_arg,nuarr)
    denom_int = integrate.trapz(denom_arg,nuarr)
    
    K4uni = num_int / denom_int
    
    return(K4uni)

def calc_kcuni(alpha,nuarr,apfarr,nu0,alpha0=-1):
    import scipy.integrate as integrate
    
    num_arg = nu0**(alpha-alpha0) * nuarr**alpha0 *apfarr
    denom_arg = nuarr**alpha * apfarr
    
    num_int = integrate.trapz(num_arg,nuarr)
    denom_int = integrate.trapz(denom_arg,nuarr)
    
    Kcuni = num_int / denom_int
    
    return(Kcuni)

def calc_k5uni(alpha,nuarr,apfarr,areaarr,nu0):

    import scipy.integrate as integrate

    num_arg=nu0**alpha * apfarr
    denom_arg=nuarr**alpha * areaarr * apfarr
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
    K5uni = num_int / denom_int
    
    return(K5uni)

def calc_kc2uni(alpha,nuarr,apfarr,areaarr,nu0,alpha0):

    import scipy.integrate as integrate

    num_arg=nu0**(alpha-alpha0) * apfarr
    denom_arg=nuarr**alpha * areaarr * apfarr
    
    num_int=integrate.trapz(num_arg,nuarr)
    denom_int=integrate.trapz(denom_arg,nuarr)
    
    Kc2uni = num_int / denom_int
    
    return(Kc2uni)
