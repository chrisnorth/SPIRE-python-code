# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 15:16:08 2012

@author: chris
"""

from numpy import arange,concatenate,pi,zeros_like,min,max
from scipy import where
from scipy.integrate import trapz
import matplotlib.pyplot as plot

def maincode():
    
    ##################################################
    ## SET PARAMETERS
    ##################################################
    #physical parameters
    c=299792458.
    
    #dust model parameters
    T0=274 #temperature a 1au (K)
    wl0=160.e-6 #knee wavelength (m)
    n0=2.2e-8 #number density at 1au (cm-3)
    
    grainem=False    ##Set to true to calculate emmisivity based on grain structure.
    
    #grain size parameters
    ab=30.e-6 #characteristic grain radius (microns)
    k1=-2.3 #grain size power law exponent(a < ab)
    k2=-5.0 #grain size power law exponent(a >= ab)
    f0=5.e-4 #normalisation of grain distribution
        
    #list of wavelengths and frequencies
    wlarr = concatenate((arange(10.,100.,1.),arange(100.,1000.,10.)))
    wlarr = wlarr * 1.e-6
    nwl=len(wlarr)
    nuarr = c / wlarr
    print min(wlarr),max(wlarr)
    ####################################################
    ## Calculate emissivities
    ####################################################
    if grainem:
        #calcualte emissivities based on grains
        emarr=zeros_like(wlarr)
        for w in range(nwl):
            pfa=False
            if w==0:
                pfa=True
            emarr[w]=int_em(wlarr[w],n0,ab,f0,k1,k2,pfa=pfa)
    else:
        #simple emissivity calculation
        emarr=where(wlarr > wl0,(wlarr/wl0)**-2,1.)
        
        emarr = emarr * n0

    plot.figure(3)
    plot.clf()
    plot.plot(wlarr*1.e6,emarr)
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel(r'$\lambda$ ($\mu$m)')
    plot.ylabel('Emissivity')
    plot.draw()
    
    # Compute Planck function
    
    barr0=planck(wlarr,T0)
    
    Zarr0=barr0 * emarr
    
    plot.figure(4)
    plot.clf()
    plot.plot(wlarr*1.e6,Zarr0)
    plot.plot(wlarr*1.e6,planck(wlarr,T0))
    plot.xscale('log')
    plot.yscale('log')
    plot.xlabel(r'$\lambda$ ($\mu$m)')
    plot.ylabel('Intensity')
    plot.draw()
    
    #######################################################
    ## Generate structure of Zodiacal Light

def planck(wl,T):

    from numpy import expm1
    #return the planck function for frequency and temperature
    # units are W/m^2/sr/Hz
    h=6.626e-34
    c=299792458.
    kB=1.38e-23
    
    nom=2.*h*c**2 / wl**5 ##units are janskys
    denom=expm1(h*c/(wl*kB*T))
    #print nom,denom
    planck = nom/denom
    
    return(planck)
    
def f_a(a,ab=30.e-6,f0=5.e-4,k1=-2.3,k2=-5.):
    #calculate grain distribution function
    
    f_a = where(a <= ab , f0*(a/ab)**k1 , f0*(a/ab)**k2)
    
    return(f_a)

def q_a_wl(a,wl):
    # calculate dust emissivity (Q function)
    # for a single wavelength
    
    from numpy import pi

    q_a_wl = where(a < wl/(2.*pi),(wl/(2.*a*pi))**-2.,1.,)
    
    return(q_a_wl)

def int_em(wl,n0,ab,f0,k1,k2,pfa=False):
    #calculate emissivitiy of all grains
    
    ##list of grain sizes
    aarr = concatenate((arange(1.,10.,1.),arange(10.,100.,10),arange(100.,1000.,100.)))
    aarr = aarr * 1.e-6 #convert to microns

    f_aarr=f_a(aarr,ab=ab,f0=f0,k1=k1,k2=k2)

    if pfa:
        plot.figure(1)
        plot.clf()
        plot.plot(aarr*1.e6,f_aarr)
        plot.yscale('log')
        plot.xscale('log')
        plot.ylim(1e-10,1)
        plot.ylabel('f(a)')
        plot.xlabel(r'a ($\mu$m)')
        
        plot.draw()

    qarr=q_a_wl(aarr,wl)

    if pfa:
        plot.figure(2)
        plot.clf()
        plot.plot(aarr*1.e6,qarr)
        plot.yscale('log')
        plot.xscale('log')
        plot.draw()
    integrand = n0 * pi * aarr**2 * qarr * f_aarr

    int_em = trapz(integrand,aarr)

    return(int_em)

if __name__=="__main__":
    maincode()
    
plot.show()
    