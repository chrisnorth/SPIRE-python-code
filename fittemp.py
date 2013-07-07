from numpy import exp,array
import matplotlib.pyplot as plot

def maincode():

    print 'Reading maps...'
    ##read maps
    
    c= 299792458.
    wls=array((250.,350.,500.)) #in microns
    freqs=c/(wls * 1.e-6) #in Hz
    obsvals=array((1.05,1.45,1.65))
    
    #set chisq parameters
    Trange=(5,15) #temperature range to scan
    beta=2. #emissivity spectral index
    Tprec=0.1 #precision on Temperature
    errors=0.1*obsvals
    ##do chi-squares test
    Tfit,chisq = dochisq(freqs,obsvals,errors,beta,Trange=Trange,Tprec=Tprec)
    
def dochisq(freqs,obs,errors,beta,Trange=(5,15),Tprec=0.5):
    
    from numpy import arange,log,exp,zeros,log10,min
    from scipy import where

    print 'Doing chi-squared test...'
    
    Tfit=arange(Trange[0],Trange[1]+Tprec,Tprec)
    nt=Tfit.size    
    nu0=freqs[2]

    na=20
    
    plot.figure(1)
    plot.clf()
    plot.figure(2)
    plot.clf()
    #plot.show()
    a_arr=zeros((nt,na))
    t_arr=zeros((nt,na))
    chisq_arr=zeros((nt,na))
    for t in range(nt):
        Temp=Tfit[t]
        t_arr[t,:]=Temp
        aini=obs[0]/greybody(1.,nu0,Temp,nu0,beta=beta)
        nuplot=arange(1.e11,1.e13,1.e10)
        nnu=nuplot.size
        Tplot=zeros(nnu)
        for n in range(nnu):
            Tplot[n]=greybody(1.,nuplot[n],Temp,nu0,beta=0.)
        print min(nuplot),min(Tplot)
        print max(nuplot),max(Tplot)
        
        print 'A_ini:',aini
        logaini=log(aini)
        logalst=arange(na)*4./na + (logaini-2.)
        alist=exp(logalst)
        print alist
        a_arr[t,:]=alist
        for a in range(na):
            for b in range(3):
                pred_val=greybody(alist[a],freqs[b],Temp,nu0,beta)
                chisq_arr[t,a]=chisq_arr[t,a] + ((obs[b]-pred_val)**2)/errors[b]**2
        
        plot.figure(1)
        plot.plot(alist,chisq_arr[t,:],label='T=%.1f K'%(Temp))

        plot.figure(2)
        plot.plot(nuplot/1.e9,Tplot,label='T=%.1f K'%(Temp))
        
    
    print 'min chi-square:',min(chisq_arr)
    print 'best-fit tempt:',t_arr[where(chisq_arr == min(chisq_arr))]
    
    plot.figure(1)
    plot.ylabel(r'$\chi^2$')
    plot.xlabel('a')
    plot.xscale('log')
    plot.yscale('log')
    #plot.legend(loc='upper right')
    plot.show()
    
    plot.figure(2)
    plot.ylabel(r'$I_\nu$')
    plot.xlabel(r'$\nu$')
    #plot.legend(loc='upper right')
    plot.show()

    plot.figure(3)
    plot.clf()
    plot.contourf(a_arr,t_arr,log10(chisq_arr),16)
    plot.xscale('log')
    plot.colorbar()
    plot.show()
    
    return Tfit,chisq_arr

def greybody(amult,nu,T,nu0,beta):
        
    b=planck(nu,T)

    i_nu = amult * b * (nu/nu0)**beta
    
    #i_nu = i_nu * 1.e-20 # downgrade units
    return(i_nu)

def planck(nu,T):

    from numpy import expm1
    #return the planck function for frequency and temperature
    
    h=6.626e-34
    c=299792458.
    kB=1.38e-23
    
    nom=1.e26 * 2.*h*nu**2 / c**2 ##units are janskys
    denom=expm1(h*nu/(kB*T))
    #print nom,denom
    planck = nom/denom
    
    return(planck)
    
if __name__=="__main__":
    maincode()
