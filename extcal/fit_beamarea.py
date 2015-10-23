# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 16:33:20 2013

@author: chris
"""

from numpy import zeros,array,arange,min,max,polyfit
from scipy import where
import string
from csv import reader

import matplotlib.pyplot as plot

def maincode():
    
    ## set input parameters
    band=range(3)
    brad=600
    ind=0.85
    
    beta1=2.0
    
    beta2=1.75
    beta3=1.5
    
    ###set file suffixes
    
    fsuff='_theta_final_newBprofBS_Br%d_Ind%.2f'%(int(brad),ind)
    fsufft='_theta_newBprofBS_temp_beta%.1f_Br%d_Ind%.2f'%(beta1,int(brad),ind)

    fsufft2='_theta_newBprofBS_temp_beta%.1f_Br%d_Ind%.2f'%(beta2,int(brad),ind)
    fsufft3='_theta_newBprofBS_temp_beta%.1f_Br%d_Ind%.2f'%(beta3,int(brad),ind)
    
    ##set fitting parameters for alpha-fit
    afita=2
    
    ##set fitting parameters for temp-fit
    afitt=[3,2,1]
    nft=len(afitt)
    tlim=[[5,16],[14,31],[29,40]]
    tlim2=[[5,15],[15,30],[30,40]]

    ##set fitting parameters for beta-fit    
    afitb=2
    
    ##############################
    ##  Spectral index fitting  ##
    ##############################
    
    ##read in beam area against alpha
    alphabm=[]
    areaa_psw=[]
    areaa_pmw=[]
    areaa_plw=[]
    areapip=zeros(3)
    fileareaa='../Outputs/beamarea'+fsuff+'.csv'
    areaa_input=reader(open(fileareaa))
    for row in areaa_input:
        if string.find(row[0],'#') < 0:
            alphabm.append(float(row[0]))
            areaa_psw.append(float(row[1]))
            areaa_pmw.append(float(row[2]))
            areaa_plw.append(float(row[3]))
            if float(row[0])==-1:
                areapip[0]=float(row[1])
                areapip[1]=float(row[2])
                areapip[2]=float(row[3])
                #print areapip
    print 'Pipeline beam areas:',areapip
    alphabm=array(alphabm)
    nabm=alphabm.size
    area_eff=zeros((nabm,3))
    area_eff[:,0]=areaa_psw
    area_eff[:,1]=areaa_pmw
    area_eff[:,2]=areaa_plw
    
    print 'Fitting beam areas against spectral index...'
    fitpar=zeros((afita+1,3))
    areaa_fit=zeros((nabm,3))
    maxoff=zeros(3)
    print 'Fit order: %d'%afita
    for b in band:
        #print 'Band %d...'%(b)
        #for a in arange(nalph):
        #    area_eff[a,b]=measbeam_new(radarr,beam_scl[:,b],beam_fix[:,b],nuarr,rsrfarr[:,b],nueff[b],alpharr[a],ind=ind)
        ##fit power law
        #fitpar[0,b],fitpar[1,b],rval,pval,stderr=linregress(log10(alphabm-alpha_nep),log10(area_eff[:,b]/areas[b]))
        fitpar[:,b]=polyfit(alphabm,area_eff[:,b]/areapip[b],deg=afita)
        for d in range(0,afita+1):
            apow=afita-d
            #print d,apow,fitpar[d,b]
            areaa_fit[:,b]=areaa_fit[:,b] + fitpar[d,b]*alphabm[:]**(apow)
        areaa_fit[:,b]=areaa_fit[:,b]*areapip[b]
        maxoff[b]=max(100.*abs(areaa_fit[:,b]-area_eff[:,b])/area_eff[:,b])
        #areaa_fit[:,b]=areapip[b]*(fitpar[0,b]*alphabm**2 + fitpar[1,b]*alphabm + fitpar[2,b])
    print 'Max offset (%%): [%.3f, %.3f, %.3f]'%(maxoff[0],maxoff[1],maxoff[2])
    
    plotareaa(alphabm,area_eff,areapip,areaa_fit,fitpar,fsuff)

    ###########################
    ##  Temperature fitting  ##
    ###########################

    tempbm=[]
    areat_psw=[]
    areat_pmw=[]
    areat_plw=[]
    tmin=5
    fileareat='../Outputs/beamarea'+fsufft+'.csv'
    areat1_input=reader(open(fileareat))
    for row in areat1_input:
        if string.find(row[0],'#') < 0:
            if float(row[0])>tmin:
                tempbm.append(float(row[0]))
                areat_psw.append(float(row[1]))
                areat_pmw.append(float(row[2]))
                areat_plw.append(float(row[3]))
            
    tempbm=array(tempbm)
    ntbm=tempbm.size
    area_efft=zeros((ntbm,3))
    area_efft[:,0]=areat_psw
    area_efft[:,1]=areat_pmw
    area_efft[:,2]=areat_plw

                
    ##Fit beam areas
    print 'Fitting beam areas against temperature...'
    fitpart=[]
    areat_fit=[]
    maxoff=zeros((nft,3))
    for f in range(nft):
        print 'Fit part %d (%.0f-%.0fK) [order: %d]:'%(f+1,tlim[f][0],tlim[f][1],afitt[f])
        fitpart.append(zeros((afitt[f]+1,3)))
        #print tlim[f][0],tlim[f][1]
        infit=where(tempbm >= tlim[f][0],True,False)
        infit=where(tempbm > tlim[f][1],False,infit)
        infit2=where(tempbm >= tlim2[f][0],True,False)
        infit2=where(tempbm > tlim2[f][1],False,infit)
        areat_fit.append(zeros((ntbm,3)))
        for b in band:
            #print 'Band %d...'%(b)
            #print tempbm[where(infit)]
            #print area_efft[where(infit),b][0]
            fitpart[f][:,b]=polyfit(tempbm[where(infit)],area_efft[where(infit),b][0]/areapip[b],deg=afitt[f])
            for d in range(0,afitt[f]+1):
                tpow=afitt[f]-d
                #print d,tpow,fitpart[f][d,b]
                areat_fit[f][:,b]=areat_fit[f][:,b] + fitpart[f][d,b]*tempbm[:]**(tpow)
            
            areat_fit[f][:,b]=areat_fit[f][:,b]*areapip[b]
            maxoff[f,b]=max(100.*abs(areat_fit[f][where(infit2),b]-area_efft[infit2,b])/area_efft[infit2,b])
        print 'Max offsets (%%): [%.3f, %.3f, %.3f]'%(maxoff[f,0],maxoff[f,1],maxoff[f,2])
        #areaa_fit[:,b]=areapip[b]*(fitpar[0,b]*alphabm**2 + fitpar[1,b]*alphabm + fitpar[2,b])
        
    plotareat(tempbm,area_efft,areapip,areat_fit,fitpart,afitt,tlim2,fsuff)
    
    ####################
    ##  Beta fitting  ##
    ####################
#        
#    areat2_psw=[]
#    areat2_pmw=[]
#    areat2_plw=[]
#    fileareat2='../Outputs/beamarea'+fsufft2+'.csv'
#    #fileareat2=fileareat1
#    areat2_input=reader(open(fileareat2))
#    for row in areat2_input:
#        if string.find(row[0],'#') < 0:
#            if float(row[0])>tmin and float(row[0])<tmax:
#                areat2_psw.append(float(row[1]))
#                areat2_pmw.append(float(row[2]))
#                areat2_plw.append(float(row[3]))
#    
#    area_efft2=zeros((ntbm,3))
#    area_efft2[:,0]=areat2_psw
#    area_efft2[:,1]=areat2_pmw
#    area_efft2[:,2]=areat2_plw

    plot.show()

def plotareaa(alpharr,areaeff,areapip,areaa_fit,fitpar,fsuff):
      
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    #area0=array([426.,771.,1626.])
    #alpha_nep=array((1.26,1.39,1.45))
    
    xtickv=arange(-4,5,1)
    nt=xtickv.size
    xtickl=[]
    xtickbl=[]
    for n in arange(nt): 
        xtickl.append(str(xtickv[n]))
        xtickbl.append('')
    
    bandStr=[r'$250\,\mu$m',r'$350\,\mu$m',r'$500\,\mu$m']
    cols=['b','g','r']
    #styles=['-','--',':']
    plot.figure(1,figsize=(12,6))
    plot.clf()
    
    plot.subplot(1,2,1)
    for b in range(3):
        plot.plot(alpharr,areaeff[:,b]/areapip[b],c=cols[b],ls='-',label=bandStr[b],lw=3)
        #plot.plot(alpha_nep[b],areas[b],c=cols[b],ls='x')
    #plot.axhline(y=areas[0],color='b',linestyle='--',label=r'$\Omega_\mathrm{Nep}$')
    #plot.axhline(y=area0[0],color='b',linestyle=':',label=r'$\Omega_0$')
    #plot.axvline(x=alpha_nep[0],color='b',linestyle=':')
    plot.ylabel('Relative Beam Sold Angle')
    plot.xlabel(r'Spectral index, $\alpha$')
    
    plot.legend(loc='lower left',ncol=1,frameon=False)

    #plot.xticks(xtickv,xtickl)
    #plot.title('PSW')
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    afita=len(fitpar[:,0])
    #print afita
    for b in range(3):
        plot.plot(alpharr,areaa_fit[:,b]/areapip[b],c='k',ls='--',lw=1)
    
    plot.subplot(1,2,2)
    plot.xticks((0,1),('',''))
    plot.yticks((0,1),('',''))
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    for b in range(3):        
        lab=r'$\Omega/\Omega_\mathrm{pip}='
        apow=afita-1
        #print 0,apow,fitpar[0,b]
        for d in range(0,afita):
            apow=afita-d-1
            #print d,apow,fitpar[d,b]
            if fitpar[d,b]<0:
                lab=lab+'-'
            elif d>0:
                lab=lab+'+'
            lab=lab+r'%.3g'%(abs(fitpar[d,b]))
            if apow>0:
                lab=lab+r'\alpha'
            if apow>1:
                lab=lab+r'^%d'%(apow)
        lab=lab+r'$'
        #lab = lab+' ; $\Omega_\mathrm{pip}=%.1f$ arcsec$^2$'%(areapip[b])
        #plot.text(x0+xr*0.05,y0+yr*0.05*(1+b),r'%s: $\Omega/\Omega_0=%.3g\alpha^2 + %.3g\alpha + %.3g$'%(bandStr[b],fitpar[0,b],fitpar[1,b],fitpar[2,b]),color=cols[b])
        plot.text(x0+xr*0.05,y1-yr*0.05*(1+b),lab,color=cols[b])
    plot.legend(loc='upper right',ncol=1,frameon=False)

    figname='../Outputs/beamarea_fit_alpha'+str(afita-1)+fsuff
    
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)


def plotareat(temparr,area_efft,areapip,areat_fit,fitpart,afitt,tlim,fsuff):
      
    ############################################################
    #Plot beam area vs radius
    
    #c=299792458. #speed of light
    #wl0=array([250.,350.,500.])*1.e-6
    #nu0=(c/wl0)/1.e9 #in GHz
    #area0=array([426.,771.,1626.])
    #alpha_nep=array((1.26,1.39,1.45))
    
    xtickv=arange(-4,5,1)
    nt=xtickv.size
    xtickl=[]
    xtickbl=[]
    for n in arange(nt): 
        xtickl.append(str(xtickv[n]))
        xtickbl.append('')
    
    bandStr=[r'$250\,\mu$m',r'$350\,\mu$m',r'$500\,\mu$m']
    #cols=['b','g','r']
    plot.figure(2,figsize=(12,6))
    plot.clf()
    plot.subplot(1,2,1)
    plot.plot(temparr,area_efft[:,0]/areapip[0],'b-',lw=3,label=r'%s'%(bandStr[0]))
    plot.plot(temparr,area_efft[:,1]/areapip[1],'g-',lw=3,label=r'%s'%(bandStr[1]))
    plot.plot(temparr,area_efft[:,2]/areapip[2],'r-',lw=3,label=r'%s'%(bandStr[2]))
    
    plot.ylabel('Relative Beam Solid Angle')
    plot.xlabel(r'Temperature [K]')
    
    plot.legend(loc='upper right',ncol=1,frameon=False)

    #plot.xticks(xtickv,xtickl)
    #plot.title('PSW')
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    
    styles=('--',':','-.')
    #print afitt
    nft=len(afitt)
    for f in range(nft):
        #print afitt[f]
        infit=where(temparr >= tlim[f][0],True,False)
        infit=where(temparr > tlim[f][1],False,infit)
        #print min(temparr[infit]),max(temparr[infit])
        plot.plot(temparr[infit],areat_fit[f][infit,0]/areapip[0],c='k',ls=styles[f],lw=1)
        plot.plot(temparr[infit],areat_fit[f][infit,1]/areapip[1],c='k',ls=styles[f],lw=1)
        plot.plot(temparr[infit],areat_fit[f][infit,2]/areapip[2],c='k',ls=styles[f],lw=1)

    plot.subplot(1,2,2)
    plot.xticks((0,1),('',''))
    plot.yticks((0,1),('',''))
    (x0,x1)=plot.xlim()
    (y0,y1)=plot.ylim()
    xr=x1-x0
    yr=y1-y0
    cols=('b','g','r')
    for f in range(nft):
        for b in range(3):
            lab=r'$\Omega/\Omega_\mathrm{pip}='
            for d in range(0,afitt[f]+1):
                tpow=afitt[f]-d
                #print d,tpow,fitpart[f][d,b]
                if fitpart[f][d,b]<0:
                    lab=lab+'-'
                elif d>0:
                    lab=lab+'+'
                lab=lab+r'%.3g'%(abs(fitpart[f][d,b]))
                if tpow>0:
                    lab=lab+r'T'
                if tpow>1:
                    lab=lab+r'^%d'%(tpow)
            lab=lab+r'$'
            #lab = lab+' ; $\Omega_\mathrm{pip}=%.1f$ arcsec$^2$'%(areapip[b])
            #plot.text(x0+xr*0.05,y0+yr*0.05*(1+b),r'%s: $\Omega/\Omega_0=%.3g\alpha^2 + %.3g\alpha + %.3g$'%(bandStr[b],fitpar[0,b],fitpar[1,b],fitpar[2,b]),color=cols[b])
            plot.text(x0+xr*0.05,y1-yr*0.05*(2+b+4*f),lab,color=cols[b])
        lab=r'Fit %d (%.0f-%.0f K):'%(f+1,tlim[f][0],tlim[f][1])
        plot.text(x0+xr*0.05,y1-yr*0.05*(1+4*f),lab,color='k')
    
    plot.legend(loc='upper right',ncol=1,frameon=False)

    figname='../Outputs/beamarea_fit_temp_'
    for f in range(nft):#
        figname=figname+str(afitt[f]-1)+'_'
    figname=figname+fsuff
    plot.savefig(figname+'.png',transparent=False,dpi=300.)
    plot.savefig(figname+'.eps',transparent=False)

    
if __name__=="__main__":
    maincode()