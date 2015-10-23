from beam import getbeam,getbeamprofile,beam_azsm,beamarea_az
from numpy import zeros,arange,array
import matplotlib.pyplot as plot
from scipy import where
import argparse

def maincode():
    parser=argparse.ArgumentParser(description='Plot Herschel-SPIRE beams')
    parser.add_argument('--bmap',action='store_true',default=False,dest='bmap',help='Read measured & Theoretical beams from map')

    args=parser.parse_args()
    bmap=args.bmap

    band=range(3)
    
    brad=1000
    radarr=arange(brad)
    nrad=radarr.size
    beam_sm_m=zeros((nrad,3))
    beam_sm_m0=zeros((nrad,3))
    beam_sm_t=zeros((nrad,3))
    beam_sm_s=zeros((nrad,3))

    bgrid=1. #grid size in arcsec
    bwid=2001 #beam width in arcsec
    #bzlim=1.e-8
    #bzval=0.
    npx=bwid/bgrid

    beams_m=zeros((npx,npx,3))
    beams_m0=zeros((npx,npx,3))
    beams_t=zeros((npx,npx,3))
    areas_m=zeros(3)
    areas_t=zeros(3)
    
    if bmap==True:
        print 'Getting beams of type "M" and "T"...'
        for b in band:
            beamx,beamy,beams_m[:,:,b]=getbeam('M',b,bgrid,bwid)
            #beamx,beamy,beams_m0[:,:,b]=getbeam('M',b,bgrid,bwid,bzlim=1.e-7,bzval=1.e-7)
            #areas_m[b]=sum(beams_m[:,:,b])
            beamx,beamy,beams_t[:,:,b]=getbeam('T',b,bgrid,bwid)
            #areas_t[b]=sum(beams_t[:,:,b])
        
        #print 'Full measured beam areas: [%.2f , %.2f , %.2f] sq.arcsec' % (areas_m[0],areas_m[1],areas_m[2])
        #print 'Full theoretical beam areas: [%.2f , %.2f , %.2f] sq.arcsec' % (areas_t[0],areas_t[1],areas_t[2])

    #    plotbeammap(beamx,beamy,beams_m,beams_t)

        print 'Making azimuthally smoothed beam...'
        for b in band:
            (beam_sm_m[:,b],j1,j2,j3)=beam_azsm(beamx,beamy,beams_m[:,:,b],radarr,retall=True)
            #(beam_sm_m0[:,b],j1,j2,j3)=beam_azsm(beamx,beamy,beams_m0[:,:,b],radarr,retall=True)
            (beam_sm_t[:,b],j1,j2,j3)=beam_azsm(beamx,beamy,beams_t[:,:,b],radarr,retall=True,pess=True)
        beam_sm_m0=where(beam_sm_m < 1.e-8,1.e-8,beam_sm_m)
    else:
        print 'Getting beam profiles of type "M" and "T"...'
        for b in band:
            beam_sm_t[:,b]=getbeamprofile('T',b,radarr)
            beam_sm_m[:,b]=getbeamprofile('M',b,radarr)
            beam_sm_m0[:,b]=getbeamprofile('M',b,radarr,bzlim=1.e-8,bzval=1.e-8)
        
    print 'Getting beam profile of type "S"...'
    for b in band:
        beam_sm_s[:,b]=getbeamprofile('S',b,radarr)

    area_m=zeros((nrad,3))
    area_m0=zeros((nrad,3))
    area_t=zeros((nrad,3))
    area_s=zeros((nrad,3))
    for r in range(brad):
        for b in band:
            area_m[r,b]=beamarea_az(radarr,beam_sm_m[:,b],brad=radarr[r])
            area_m0[r,b]=beamarea_az(radarr,beam_sm_m0[:,b],brad=radarr[r])
            area_t[r,b]=beamarea_az(radarr,beam_sm_t[:,b],brad=radarr[r])
            area_s[r,b]=beamarea_az(radarr,beam_sm_s[:,b],brad=radarr[r])

    area0=array([426.,771.,1626.])

    print 'Measured:',area_m[where(radarr == 350.),:]
    print 'Theoretical:',area_t[where(radarr == 350.),:]
    print 'Combined:',area_s[where(radarr == 350.),:]

    #plotbeamprof(radarr,beam_sm_m,beam_sm_t,beam_sm_s)
    #plotbeamrad(radarr,area_m,area_t,area_s,area0)

    plotbeamprof2(radarr,beam_sm_m,beam_sm_t,beam_sm_s,area_m,area_t,area_s,area0)

    pos_trans=array(([250,301,289],[230,255,240],[225,260,254]))
    ##pos_trans: [band,position], position=[range min,range max, transition point]
    
    plottrans(radarr,beam_sm_m0,beam_sm_s,pos_trans)

    file_out='../Outputs/beamarea_rad.dat'
    file_out0='../Outputs/beamarea-lim1e-8_rad.dat'
    prof_out='../Outputs/beamprof_rad.dat'
    prof_out0='../Outputs/beamprof-lim1e-8_rad.dat'
    f=open(file_out,'w')
    f0=open(file_out0,'w')
    p=open(prof_out,'w')
    p0=open(prof_out0,'w')
        
    line='#rad , PSW , PMW , PLW \n'
    f.write(line)
    f0.write(line)
    for r in arange(brad):
        line='%.1f , %.4f , %.4f , %.4f \n'%(radarr[r],area_m[r,0],area_m[r,1],area_m[r,2])
        line0='%.1f , %.4f , %.4f , %.4f \n'%(radarr[r],area_m0[r,0],area_m0[r,1],area_m0[r,2])
        prof='%.1f , %.4e , %.4e , %.4e \n'%(radarr[r],beam_sm_m[r,0],beam_sm_m[r,1],beam_sm_m[r,2])
        prof0='%.1f , %.4e , %.4e , %.4e \n'%(radarr[r],beam_sm_m0[r,0],beam_sm_m0[r,1],beam_sm_m0[r,2])
        f.write(line)
        f0.write(line0)
        p.write(prof)
        p0.write(prof0)
    
    
    f.close()
    f0.close()
    p.close()
    p0.close()
        
    plot.show()

def plotbeammap(beamx,beamy,beams_m,beams_t,plimin=None):
    
    print 'Plotting beam maps...'
    
    from numpy import arange,floor,ceil,array,log10,min,max
    from scipy import where
    from pylab import NaN,inf
    
    print 'Plotting beams...'
    print 'WARNING: this can be slow'
    #plot beams
    plot.figure(1,figsize=(6,12))
    plot.clf()
    plot.hot()

    titles=('PSM Beam','PMW Beam','PLW Beam')
    
    for b in arange(3):
        print 'plot measured beam band %d'%b
        plot.subplot(3,2,2*b+1)
        plot.axis('equal')
        plot.gca().set_axis_bgcolor("#a5a5a5")
        bplot=log10(beams_m[:,:,b])
        #bplot=where(bplot == -pylab.inf,pylab.NaN,bplot) #replace -inf with NaN
        #bplot=where(bplot == pylab.inf,pylab.NaN,bplot) #replace inf with NaN
        #minb=min(bplot[where(bplot > -1.e30 )])
        #maxb=max(bplot[where(bplot < 1.e30 )])
        #print 'minb,maxb',minb,maxb
        print 'plimin',plimin
        if plimin == None:
            plim=array([-8,0])
        else: plim=array([plimin,0])
        print 'plim',plim
        (plim[0],plim[1])=(int(plim[0]),int(plim[1]))
        bplot=where(bplot < plim[0],plim[0]+1,bplot)
        bplot=where(bplot == NaN,plim[0]+1,bplot)
        bplot=where(bplot == -inf,plim[0]+1,bplot)
        #bplot=where(bplot != bplot,plim[0],bplot) #replace NaNs with min
        #bplot=where(bplot < plim[0],plim[0],bplot) #replace min vals with min
        #bplot=where(bplot > plim[1],plim[1],bplot) #replace max vals with max
        prange=plim[1]-plim[0]
        #ntick=8
        ncol=8
        #tlevs=arange(plim[0],plim[1]+prange/ntick,prange/ntick)
        clevs=arange(plim[0],plim[1]+prange/ncol,prange/ncol)
        #print 'Ticks:',tlevs
        print 'Colours:',clevs
        plot.contourf(beamx,beamy,bplot,levels=clevs)
        #plot.contourf(beamx,beamy,bplot)
        plot.clim(plim[0],plim[1])
        plot.xlim(min(beamx),max(beamx))
        plot.ylim(min(beamy),max(beamy))
        #plot.title(titles[b])
        if b == 0:
            plot.title('Measured Beam')

        if b == 2:
            plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))))
            plot.xlabel('x ["]')
        else:
            plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))),('','',''))

        plot.yticks((int(min(beamy[0,:])),0,int(max(beamy[0,:]))))
        plot.ylabel('y ["]')
        
        #cb=plot.colorbar(format='%.1f')#,ticks=tlevs)
        #cb.set_label('Log scale')



        print 'plot theoretical beam band %d'%b
        plot.subplot(3,2,2*b+2)
        plot.axis('equal')
        plot.gca().set_axis_bgcolor("#a5a5a5")
        bplot=log10(beams_t[:,:,b])
        #bplot=where(bplot == -pylab.inf,pylab.NaN,bplot) #replace -inf with NaN
        #bplot=where(bplot == pylab.inf,pylab.NaN,bplot) #replace inf with NaN
        #minb=min(bplot[where(bplot > -1.e30 )])
        #maxb=max(bplot[where(bplot < 1.e30 )])
        #print 'minb,maxb',minb,maxb
        print 'plimin',plimin
        if plimin == None:
            plim=array([-8,0])
        else: plim=array([plimin,0])
        print 'plim',plim
        (plim[0],plim[1])=(int(plim[0]),int(plim[1]))
        #bplot=where(bplot != bplot,plim[0],bplot) #replace NaNs with min
        #bplot=where(bplot < plim[0],plim[0],bplot) #replace min vals with min
        #bplot=where(bplot > plim[1],plim[1],bplot) #replace max vals with max
        prange=plim[1]-plim[0]
        #ntick=8
        ncol=8
        #tlevs=arange(plim[0],plim[1]+prange/ntick,prange/ntick)
        clevs=arange(plim[0],plim[1]+prange/ncol,prange/ncol)
        #print 'Ticks:',tlevs
        print 'Colours:',clevs
        plot.contourf(beamx,beamy,bplot,levels=clevs)
        #plot.contourf(beamx,beamy,bplot)
        plot.clim(plim[0],plim[1])
        plot.xlim(min(beamx),max(beamx))
        plot.ylim(min(beamy),max(beamy))

        if b == 0:
            plot.title('Theoretical Beam')
        
        if b == 2:
            plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))))
            plot.xlabel('x ["]')
        else:
            plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))),('','',''))
        
        plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:]))),('','',''),rotation=90)
        plot.ylabel(titles[b])

        cb=plot.colorbar(format='%.1f')#,ticks=tlevs)
        cb.set_label('Log scale')

        
#    print 'plot PMW'
#    plot.subplot(2,2,2)
#    plot.axis('equal')
#    plot.contourf(beamx,beamy,bplot[:,:,1]-max(bplot[:,:,1]),levels=levs)
#    plot.clim(plim[0],plim[1])
#    plot.title('PMW Beam')
#    plot.xticks((floor(min(beamx[:,0])),0,ceil(max(beamx[:,0]))))
#    #plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:])))),rotation=90)
#    plot.yticks((floor(min(beamy[0,:])),0,ceil(max(beamy[0,:]))),('','',''),rotation=90)
#    cb=plot.colorbar(format='%.1f',ticks=levs)
#    cb.set_label('Log scale')
#    plot.xlabel('x ["]')
#        
#    print 'plot PLW'
#    plot.subplot(2,2,3)
#    plot.axis('equal')
#    plot.contourf(beamx,beamy,bplot[:,:,2]-max(bplot[:,:,2]),levels=levs)
#    plot.clim(plim[0],plim[1])
#    plot.title('PLW Beam')
#    plot.xticks((floor((min(beamx[:,0]))),0,ceil(max(beamx[:,0]))))
#    plot.yticks((floor((min(beamy[0,:]))),0,ceil(max(beamy[0,:]))),rotation=90)
#    #plot.colorbar()
#    cb=plot.colorbar(format='%.1f',ticks=levs)
#    cb.set_label('Log scale')
#    #colorbar.set_label('Log scale')
#    plot.xlabel('x ["]')
#    plot.ylabel('y ["]')
        
    plot.savefig('../Outputs/beams_meas_theo.png',transparent=False,dpi=300.)
    #plot.savefig('../Outputs/beams'+fsuff+'.eps',transparent=False)
    
    
def plotbeamprof(radarr,beam_sm_m,beam_sm_t,beam_sm_s):
    
    print 'Plotting beam profiles...'
    
    plot.figure(2)
    plot.clf()
    lm,=plot.plot(radarr,beam_sm_m[:,0],'k:')    
    plot.plot(radarr,beam_sm_m[:,0],'b:')
    plot.plot(radarr,beam_sm_m[:,1],'g:')
    plot.plot(radarr,beam_sm_m[:,2],'r:')
    
    lt,=plot.plot(radarr,beam_sm_t[:,0],'k--')
    plot.plot(radarr,beam_sm_t[:,0],'b--')
    plot.plot(radarr,beam_sm_t[:,1],'g--')
    plot.plot(radarr,beam_sm_t[:,2],'r--')
    
    ls,=plot.plot(radarr,beam_sm_s[:,0],'k-')
    lpsw,=plot.plot(radarr,beam_sm_s[:,0],'b-')
    lpmw,=plot.plot(radarr,beam_sm_s[:,1],'g-')
    lplw,=plot.plot(radarr,beam_sm_s[:,2],'r-')
    
    plot.yscale('log')
    plot.xlabel('Radius [arcsec]')
    plot.ylabel('Beam response')
    plot.xlim(xmax=600)
    plot.ylim(1.e-8,1)
    
    plot.legend((lm,lt,ls,lpsw,lpmw,lplw),('Measured','Theoretical','Combined','PSW','PMW','PLW'),ncol=2,loc='upper right')    
    plot.gca().lines.remove(lm)
    plot.gca().lines.remove(lt)
    plot.gca().lines.remove(ls)
    
    plot.savefig('../Outputs/beam_profiles.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/beam_profiles.eps',transparent=False)
    
def plotbeamrad(radarr,area_m,area_t,area_s,area0):

    print 'Plotting beam areas...'
    
    plot.figure(3)
    plot.clf()
    
    #plot.subplot(2,1,1)
    plot.plot(radarr,area_m[:,0],'b:')
    plot.plot(radarr,area_m[:,1],'g:')
    plot.plot(radarr,area_m[:,2],'r:')

    plot.plot(radarr,area_t[:,0],'b--')
    plot.plot(radarr,area_t[:,1],'g--')
    plot.plot(radarr,area_t[:,2],'r--')

    plot.plot(radarr,area_s[:,0],'b-')
    plot.plot(radarr,area_s[:,1],'g-')
    plot.plot(radarr,area_s[:,2],'r-')
    
    plot.axhline(area0[0],linestyle='-.',color='b')
    plot.axhline(area0[1],linestyle='-.',color='g')
    plot.axhline(area0[2],linestyle='-.',color='r')
    
    plot.savefig('../Outputs/beam_radarea.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/beam_radarea.eps',transparent=False)


def plotbeamprof2(radarr,beam_sm_m,beam_sm_t,beam_sm_s,area_m,area_t,area_s,area0):
    
    import matplotlib.font_manager
    
    print 'Plotting beam profiles...'
    
    fontp=matplotlib.font_manager.FontProperties()
    fontp.set_size('small')            
    
    plot.figure(4,figsize=(6,16))
    plot.clf()
    plot.subplot(3,1,1)
    lm,=plot.plot(radarr,beam_sm_m[:,0],'k:')    
    plot.plot(radarr,beam_sm_m[:,0],'b:')
    plot.plot(radarr,beam_sm_m[:,1],'g:')
    plot.plot(radarr,beam_sm_m[:,2],'r:')
    
    lt,=plot.plot(radarr,beam_sm_t[:,0],'k--')
    plot.plot(radarr,beam_sm_t[:,0],'b--')
    plot.plot(radarr,beam_sm_t[:,1],'g--')
    plot.plot(radarr,beam_sm_t[:,2],'r--')
    
    ls,=plot.plot(radarr,beam_sm_s[:,0],'k-')
    lpsw,=plot.plot(radarr,beam_sm_s[:,0],'b-')
    lpmw,=plot.plot(radarr,beam_sm_s[:,1],'g-')
    lplw,=plot.plot(radarr,beam_sm_s[:,2],'r-')
    
    plot.yscale('log')
    #plot.xlabel('Radius [arcsec]')
    plot.ylabel(r'Beam response, $B(\theta)$')
    plot.xlim(xmax=600)
    plot.ylim(1.e-8,1)
    
    #plot.legend((lm,lt,ls,lpsw,lpmw,lplw),('Measured','Theoretical','Combined','PSW','PMW','PLW'),ncol=2,loc='upper right')
    plot.legend((lm,lpsw,lt,lpmw,ls,lplw),('Meas','PSW','Theo','PMW','Comb','PLW'),loc='upper right',ncol=3,columnspacing=1,prop=fontp)
    plot.gca().lines.remove(lm)
    plot.gca().lines.remove(lt)
    plot.gca().lines.remove(ls)
    
    plot.subplot(3,1,2)
    plot.plot(radarr,area_m[:,0],'b:')
    plot.plot(radarr,area_m[:,1],'g:')
    plot.plot(radarr,area_m[:,2],'r:')

    plot.plot(radarr,area_t[:,0],'b--')
    plot.plot(radarr,area_t[:,1],'g--')
    plot.plot(radarr,area_t[:,2],'r--')

    plot.plot(radarr,area_s[:,0],'b-')
    plot.plot(radarr,area_s[:,1],'g-')
    plot.plot(radarr,area_s[:,2],'r-')
    
    plot.axhline(area0[0],linestyle='-.',color='b')
    plot.axhline(area0[1],linestyle='-.',color='g')
    plot.axhline(area0[2],linestyle='-.',color='r')

    plot.ylabel(r'$\omega(\theta)$ [sq. arcsec]')
    plot.xlim(xmax=600)
    
    nr=radarr.size
    amax_m=area_m[499,:]
    amax_t=area_t[499,:]
    amax_s=area_s[499,:]

    a350_m=area_m[350,:]
    a350_t=area_t[350,:]
    a350_s=area_s[350,:]

    print 'Omega_0:',area0
    print ''
    print 'Measured (500''):',amax_m
    print 'Theoretical (500''):',amax_t
    print 'Combined (500''):',amax_s
    print ''
    print 'Measured (350''):',a350_m
    print 'Theoretical (350''):',a350_t
    print 'Combined (350''):',a350_s
    
    plot.subplot(3,1,3)
    plot.axhline(0.,color='k')
    
    #lm,=plot.plot(radarr,(area_m[:,0]/amax_m[0]),'k:',label='Meas')    
    plot.plot(radarr,(area_m[:,0]/amax_m[0]),'b:')
    plot.plot(radarr,(area_m[:,1]/amax_m[1]),'g:')
    plot.plot(radarr,(area_m[:,2]/amax_m[2]),'r:')

    #lt,=plot.plot(radarr,(area_t[:,0]/amax_m[0]),'k--',label='Theo')
    plot.plot(radarr,(area_t[:,0]/amax_m[0]),'b--')
    plot.plot(radarr,(area_t[:,1]/amax_m[1]),'g--')
    plot.plot(radarr,(area_t[:,2]/amax_m[2]),'r--')

    #ls,=plot.plot(radarr,(area_s[:,0]/amax_m[0]),'k-',label='Comb')
    #lpsw,=plot.plot(radarr,(area_s[:,0]/amax_m[0]),'b-',label='PSW')
    #lpmw,=plot.plot(radarr,(area_s[:,1]/amax_m[1]),'g-',label='PMW')
    #lplw,=plot.plot(radarr,(area_s[:,2]/amax_m[2]),'r-',label='PLW')
    plot.plot(radarr,(area_s[:,0]/amax_m[0]),'b-')
    plot.plot(radarr,(area_s[:,1]/amax_m[1]),'g-')
    plot.plot(radarr,(area_s[:,2]/amax_m[2]),'r-')
    
    #lc,=plot.plot(radarr,(area_s[:,0]/amax_m[0]),'k-.',label=r'$\Omega_0$')
    plot.axhline((area0[0]/amax_m[0]),linestyle='-.',color='k')
    plot.axhline((area0[0]/amax_m[0]),linestyle='-.',color='b')
    plot.axhline((area0[1]/amax_m[1]),linestyle='-.',color='g')
    plot.axhline((area0[2]/amax_m[2]),linestyle='-.',color='r')
    
    plot.ylim(0.9,1.05)
    plot.ylabel(r'$\omega(\theta)\,/\,\omega(\theta=500'')$')
    plot.xlabel(r'Radius, $\theta$ [arcsec]')
    plot.xlim(xmax=600)
    
    #plot.legend((lm,lt,ls),('Meas','Theo','Comb'),ncol=3,loc='upper right')
    #plot.figlegend((lm,lt,ls,lpsw,lpmw,lplw),('Meas','Theo','Comb','PSW','PMW','PLW'),loc='upper center',ncol=2)
    #plot.figlegend((lm,lpsw,lt,lpmw,ls,lplw,lc),('Meas','PSW','Theo','PMW','Comb','PLW',r'$\Omega_0$'),loc='upper center',ncol=4,columnspacing=1,prop=fontp)
    #plot.legend(loc='upper right')
    #plot.gca().lines.remove(lm)
    #plot.gca().lines.remove(lt)
    #plot.gca().lines.remove(ls)
    #plot.gca().lines.remove(lc)
    
    plot.savefig('../Outputs/beam_profiles_area.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/beam_profiles_area.eps',transparent=False)
    
def plottrans(radarr,beam_sm_m,beam_sm_s,pos_trans):
    
    from scipy.interpolate import interp1d
    plot.figure(5)
    plot.clf()
    

    lm,=plot.plot(radarr,beam_sm_m[:,0],'k--',label='Meas')
    plot.plot(radarr,beam_sm_m[:,0],'b--')
    plot.plot(radarr,beam_sm_m[:,1],'g--')
    plot.plot(radarr,beam_sm_m[:,2],'r--')
    
    lc,=plot.plot(radarr,beam_sm_s[:,0],'k-',label='Comb')
    plot.plot(radarr,beam_sm_s[:,0],'b-',label='PSW')
    plot.plot(radarr,beam_sm_s[:,1],'g-',label='PMW')
    plot.plot(radarr,beam_sm_s[:,2],'r-',label='PLW')

    b_trans=zeros((3,3))
    for b in range(3):
        spline=interp1d(radarr,beam_sm_s[:,b])
        b_trans[b,0]=spline(pos_trans[b,0])
        b_trans[b,1]=spline(pos_trans[b,1])
        b_trans[b,2]=spline(pos_trans[b,2])
        
    b_trans[:,0:2]=1.e-4        
    plot.plot(pos_trans[0,0],b_trans[0,0],'b>')
    plot.plot(pos_trans[1,0],b_trans[1,0],'g>')
    plot.plot(pos_trans[2,0],b_trans[2,0],'r>')
    
    plot.plot(pos_trans[0,1],b_trans[0,1],'b<')
    plot.plot(pos_trans[1,1],b_trans[1,1],'g<')
    plot.plot(pos_trans[2,1],b_trans[2,1],'r<')
   
    plot.axvline(pos_trans[0,0],color='b',linestyle=':')
    plot.axvline(pos_trans[0,1],color='b',linestyle=':')
    plot.axvline(pos_trans[1,0],color='g',linestyle=':')
    plot.axvline(pos_trans[1,1],color='g',linestyle=':')
    plot.axvline(pos_trans[2,0],color='r',linestyle=':')
    plot.axvline(pos_trans[2,1],color='r',linestyle=':')
    
    plot.plot(pos_trans[0,2],b_trans[0,2],'bD')
    plot.plot(pos_trans[1,2],b_trans[1,2],'gD')
    plot.plot(pos_trans[2,2],b_trans[2,2],'rD')

    plot.yscale('log')
    plot.xlim(200,350)
    plot.ylim(1.e-7,1.e-3)
    plot.ylabel(r'Beam response, $B(\theta)$')
    plot.xlabel(r'Radius, $\theta$ [arcsec]')
    
    plot.legend(loc='upper right',ncol=3)
    plot.gca().lines.remove(lc)
    plot.gca().lines.remove(lm)
    plot.savefig('../Outputs/beam_profiles_transition.png',transparent=False,dpi=300)
    plot.savefig('../Outputs/beam_profiles_transition.eps',transparent=False)
    
if __name__ == "__main__":
    maincode()