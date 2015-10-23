def getrsrf(rsrftype,band,nuarr,verbose=None):
    # Get RSRF

    from numpy import array,max,zeros
    from scipy import where
    from csv import reader
    from sys import exit
    from scipy.interpolate import interp1d

    file='../Inputs/SPIRE-Phot-RSRF.csv'
    rsrf_input=reader(open(file))
    wllist=[]
    rsrf_fromfile=[]
    for row in rsrf_input:
        wllist.append(float(row[0]))
        rsrf_fromfile.append(float(row[band+1]))

    c=299792458. #speed of light
    
    rsrfin=array(rsrf_fromfile)
    nuin=c/((array(wllist))*1.e-6)
    #dnuin=zeros(len(nuin)-1)
    #for i in range(len(dnuin)):
    #    dnuin[i]=nuin[i+1]-nuin[i]
    #    print nuin[i]/1.e9,dnuin[i]/1.e9,(wllist[i]-wllist[i+1])
    #print min(dnuin)/c,max(dnuin)/c
    #print nuin[0:10]/1.e9
#    print 'nuin range:  [%f,%f]' % (nuin.min(),nuin.max())
#    print 'nuarr range: [%f,%f]' % (nuarr.min(),nuarr.max())

    #nnu=nuarr.size
    intf=interp1d(nuin,rsrfin,bounds_error=False,fill_value=0)
    rsrf=intf(nuarr)

    if rsrftype == 'T':
        #Top Hat
        if verbose: print 'Generating Top Hat RSRF'
#        print type(rsrf)
        rsrf_max=max(rsrf)
#        print type(rsrf_max)
#        print 'max=%f' % rsrf_max
        rsrflim=0.5*rsrf_max

        rsrf=where(rsrf >= rsrflim,0.5,0.)

        #bandmin=array((490.,640.,1130.)) * 1.e9
        #bandmax=array((620.,1010.,1402.)) * 1.e9

        #rsrf=zeros(nnu)

        #gtl=sci.where((nuarr >= bandmin[band]),1,0)
        #ltl=sci.where((nuarr <= bandmax[band]),1,0)
        #inband=sci.where(gtl + ltl == 2)
        #rsrf[inband]=1
        
    elif rsrftype == 'M':
        #Get measured
        if verbose: print 'Generating measured RSRF'
    else:
        error='ERROR: Unknown RSRF type "%s". Must be one of M|T.' % rsrftype
        exit(error)
    
    return(rsrf)

def getapf(apftype,band,nuarr,verbose=None):

    from numpy import array
    from csv import reader
    from scipy.interpolate import interp1d
    from sys import exit

    if apftype == 'U':
        #Uniform aperture function
        if verbose: print 'Generating uniform aperture function'

        nnu=nuarr.size
        apfarr=array((nnu))
        apfarr[:]=1.
    elif apftype == 'R':
        #Real aperture function
        if verbose: print 'Generating real aperture function'
        file='../Inputs/app_eff.csv'
        apf_input=reader(open(file))
        nulist=[]
        apf_fromfile=[]
        for row in apf_input:
            nulist.append(float(row[0])) #in THz
            apf_fromfile.append(float(row[band+1]))
            
        apfin=array(apf_fromfile)
        nuin=(array(nulist))*1.e12
              
#        print 'nuin range:  [%f,%f]' % (nuin.min(),nuin.max())
#        print 'nuarr range: [%f,%f]' % (nuarr.min(),nuarr.max())
        
        nnu=nuarr.size
        intf=interp1d(nuin,apfin,bounds_error=False,fill_value=0)
        apfarr=intf(nuarr)

    else:
        error='ERROR: Unknown Aperture function type "%s". Must be one of U|R.' % apftype
        exit(error)

    return(apfarr)

def bandedge(nuarr,rsrfin,fact=1.,method='max',lim=None):
    
    from scipy import where
    from numpy import min,max,array,median
    
    #calculate points which have half max value
    lim=max(rsrfin)/2.
    if method=='fixed': lim=lim
    if fact < 1.:
        print 'Warning: using smaller bands'

    imax=nuarr.size
    inband=where(rsrfin > lim)
    ilim=array([min(inband),max(inband)])
    if fact != 1:
        nin=ilim[1]-ilim[0]
        ilim[0]=ilim[0]-nin*(fact-1.)/2
        if ilim[0] < 0: ilim[0]=0
        ilim[1]=ilim[1]+nin*(fact-1.)/2
        if ilim[1] > imax: ilim[1]=imax
    nulim=nuarr[ilim]

    if method=='median':
        #repeat using half median value of in-band RSRF as limit
        lim=median(rsrfin[ilim[0]:ilim[1]])/2.
        inband=where(rsrfin > lim)
        ilim=array([min(inband),max(inband)])
        if fact != 1:
            nin=ilim[1]-ilim[0]
            ilim[0]=ilim[0]-nin*(fact-1.)/2
            if ilim[0] < 0: ilim[0]=0
            ilim[1]=ilim[1]+nin*(fact-1.)/2
            if ilim[1] > imax: ilim[1]=imax
        nulim=nuarr[ilim]
    
    return(ilim,nulim)

def getapeffuni(fileApEff):
    
    from csv import reader
    from numpy import array
    import string
    from scipy.interpolate import interp1d
    
    fapEff=reader(open(fileApEff))
    nuArr=[]
    apEffArr=[]
    for row in fapEff:
        if string.find(row[0],'#') < 0:
            nuArr.append(float(row[0]))
            apEffArr.append(float(row[1]))
    nuArr=array(nuArr)
    apEffArr=array(apEffArr)
    
    apEffInt=interp1d(nuArr,apEffArr)
    
    return(apEffInt)

def readsvg(fileSVG,wnin=None,verbose=False,dofilt=True):
    
    import string
    from numpy import zeros,min,max,array
    from scipy.interpolate import interp1d
    from scipy import where
    
    f=open(fileSVG,'r')
    lines=f.readlines()
    nl=len(lines)
    
    detspecin=[]
    meanspecin=[]
    filt=[]
    l=-1
    nfiltered=0
    prevco=zeros((2,2))
    while l<nl:
        if string.find(lines[l],'pxmin=') != -1:
            pxmin=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'dxmin=') != -1:
            dxmin=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'pxmax=') != -1:
            pxmax=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'dxmax=') != -1:
            dxmax=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'pymin=') != -1:
            pymin=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'dymin=') != -1:
            dymin=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'pymax=') != -1:
            pymax=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'dymax=') != -1:
            dymax=float(string.split(lines[l],'=')[1])
        if string.find(lines[l],'filt') != -1:
            filt.append([float(string.split(lines[l],',')[1]),float(string.split(lines[l],',')[2])])
            if verbose:
                print 'filter: ',filt[-1]
        dloc=string.find(lines[l],'d="m')
        if dloc !=-1:
            ##readline
            coordstr=string.split(lines[l][dloc+5:-2],' ')
            ncoord=len(coordstr)
            coord=zeros((ncoord,2))
            #print string.split(coordstr[0],',')[0], string.split(coordstr[0],',')[1]
            coord[0,0]=float(string.split(coordstr[0],',')[0])
            coord[0,1]=float(string.split(coordstr[0],',')[1])
            for c in range(1,ncoord):
                #print coordstr[c]
                coord[c,0]=float(string.split(coordstr[c],',')[0])+coord[c-1,0]
                coord[c,1]=float(string.split(coordstr[c],',')[1])+coord[c-1,1]
            #print '--'
            if verbose:
                #print 'prevco: ',prevco,len(prevco)
                if coord[0,0] != prevco[0,0]:
                    print ncoord,' coordinates'
                    print 'start/end: ',coord[0,:],':',coord[ncoord-1,:],'(',coordstr[ncoord-1],')'                
                    print 'diff: %.5f,%.5f'%(coord[0,0]-prevco[1,0],coord[0,1]-prevco[1,1])
                    #print 'prevco: ',prevco
            prevco[0,:]=coord[0,:]
            prevco[1,:]=coord[-1,:]
            coord[:,0] = (coord[:,0] - pxmin)/(pxmax-pxmin) * (dxmax-dxmin) + dxmin
            coord[:,1] = (coord[:,1] - pymin)/(pymax-pymin) * (dymax-dymin) + dymin
            l=l+1
            foundcol=False
            while not foundcol:
                if string.find(lines[l],'stroke:#000000') != -1:
                    #found colour statement
                    foundcol=True
                    wnlist=coord[:,0]
                    if dofilt:
                        ##remove lines from filtered sections
                        filtcoord=zeros(ncoord)
                        for f in range(len(filt)):
                            filtcoord=where(coord[:,0]>filt[f][0],filtcoord+1,filtcoord)
                            filtcoord=where(coord[:,0]>filt[f][1],filtcoord-1,filtcoord)
                        nofilt=where(filtcoord == 0)
                        #print ncoord,filtcoord
                        #print 'no filt:',len(nofilt[0]),nofilt[0]
                        infilt=where(filtcoord > 0)
                        #print 'in filt:',len(infilt[0]),infilt[0]
                        #coord[infilt,1]=-1.
                        #coord=coord[nofilt[0],:]
                        meanspecin.append(coord[nofilt[0],:])
                    else:
                        meanspecin.append(coord)
                    if verbose:
                        print '(mean line)'
                        print '---'
                    
                elif string.find(lines[l],'stroke:#ff0000') != -1:
                    foundcol=True
                    detspecin.append(coord)
                    if verbose:
                        print '(detector line)'
                        print '---'
                if not foundcol:
                    l=l+1
        l=l+1
    
    ndet=len(detspecin)
    if wnin == None:
        wn=wnlist[1:-1]
    else:
        wn=wnin
    intmean=interp1d(meanspecin[0][:,0],meanspecin[0][:,1],bounds_error=True)
    meanspec=intmean(wn)
    if verbose:
        print 'wn:',min(wn),max(wn)
    detspec=zeros((len(wn),ndet))
    for d in range(ndet):
        spec=detspecin[d]
        intspec=interp1d(spec[:,0],spec[:,1],bounds_error=True)
        if verbose:
            print 'spec:',min(spec[:,0]),max(spec[:,0])
        detspec[:,d]=intspec(wn)

    return(wn,meanspec,detspec)
