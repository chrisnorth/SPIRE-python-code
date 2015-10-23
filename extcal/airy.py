def airyconv(xArr,yArr,airyArr,wid,radList):

    from scipy.special import j1,j0
    from numpy import pi,sum,ones,arange
    from scipy import where
    from scipy.signal import convolve2d
    from beam import beam_azsm

    nRad=len(radList)
    dd=radList[2]-radList[1]
    #print 'wid:',wid
    #print 'dd:',dd
    radListTot=arange(0,300,0.01)
    radListTot=where(radListTot == 0.,1.e-12,radListTot)
    airyList= ((2./pi) * (j1(pi*radListTot)/radListTot))**2
    
    airyTot=sum(airyList*2.*pi*radListTot*dd)
    #print 'airyTot:',airyTot
    npKer=wid / dd
    kernel=ones((npKer,npKer))
    #kernel = kernel/size(kernel)
    kernel = kernel * (dd**2.) / airyTot

    convSame=convolve2d(airyArr,kernel,mode="same")    
    convAz=beam_azsm(xArr,yArr,convSame,radList,nsamp=360)

    return(convAz)

def airyconv_MG(wid,radList):
    
    from scipy.special import j1,j0
    from numpy import pi,sum,ones,arange,zeros,sqrt
    from scipy import where

    nRad=len(radList)
    dd=radList[2]-radList[1]
    #print 'wid:',wid
    #print 'dd:',dd
    radListTot=arange(0,300,0.01)
    radListTot=where(radListTot == 0.,1.e-12,radListTot)
    airyList= ((2./pi) * (j1(pi*radListTot)/radListTot))**2
    
    airyTot=sum(airyList*2.*pi*radListTot*dd)

    conv1=zeros(nRad)
    conv2=zeros(nRad)
    npKer=wid/dd
    xList0=arange(-wid/2.,wid/2.+dd,dd)
    npKer=len(xList0)
    xArr0=zeros((npKer,npKer))
    yArr0=zeros((npKer,npKer))
    for x in range(npKer):
        xArr0[x,:]=xList0[x]
        yArr0[:,x]=xList0[x]
    
    for r in range(nRad):
        ##move along axis
        radArr=sqrt((xArr0+radList[r])**2 + yArr0**2)
        radArr=where(radArr == 0.,1.e-12,radArr) ##prevent div/0 errors
        bessArr = j1(pi*radArr)
        airyArr = ((2./pi) * (bessArr / radArr))**2
        conv1[r] = sum(airyArr)*dd**2./airyTot
    
        ##move at 45 degrees
        radArr=sqrt((xArr0+radList[r]/sqrt(2))**2 + (yArr0+radList[r]/sqrt(2))**2)
        radArr=where(radArr == 0.,1.e-12,radArr) ##prevent div/0 errors
        bessArr = j1(pi*radArr)
        airyArr = ((2./pi) * (bessArr / radArr))**2
        conv2[r] = sum(airyArr)*dd**2./airyTot
    
    convAz=(conv1 + conv2)/2.
    
    return(conv1,conv2,convAz)
    
def getairyconv(fileConv,norm=True):
    
    from csv import reader
    import string
    from numpy import array,zeros
    from scipy.interpolate import interp1d
    
    import matplotlib.pyplot as plot
    
    nuNorm=[]
    wlNorm=[]
    widNorm=[]
    fwhmList=[]
    apEffList=[]
    beamAreaList=[]
    
    fileConv='../Inputs/airy_conv.dat' ##hard-coded
    fConv=reader(open(fileConv,'r'))
    for row in fConv:
        if string.find(row[0],'#') < 0:
            nuNorm.append(float(row[0]))
            wlNorm.append(float(row[1]))
            widNorm.append(float(row[2]))
            fwhmList.append(float(row[3]))
            apEffList.append(float(row[4]))
            beamAreaList.append(float(row[5]))

    nuNorm=array(nuNorm)
    wlNorm=array(wlNorm)
    widNorm=array(widNorm)
    fwhmList=array(fwhmList)
    apEffList=array(apEffList)
    beamAreaList=array(beamAreaList)
    
    fwhmInt=interp1d(nuNorm,fwhmList,kind='cubic')
    apEffInt=interp1d(nuNorm,apEffList,kind='cubic')
    beamAreaInt=interp1d(nuNorm,beamAreaList,kind='cubic')

    if norm:
        fwhm0=fwhmInt(1.)
        apEff0=apEffInt(1.)
        beamArea0=beamAreaInt(1.)
        fwhmNorm=interp1d(nuNorm,fwhmList/fwhm0,kind='cubic')
        apEffNorm=interp1d(nuNorm,apEffList/apEff0,kind='cubic')
        beamAreaNorm=interp1d(nuNorm,beamAreaList/beamArea0,kind='cubic')

    if norm:
        return(fwhmNorm,apEffNorm,beamAreaNorm)
    else:
        return(fwhmInt,apEffInt,beamAreaInt)

    plot.figure(20)
    plot.clf()
    plot.plot(nuNorm,fwhmList,'ko')
    plot.plot(nuNorm,fwhmInt(nuNorm),'k-')
    plot.title('FWHM')
    
    plot.figure(21)
    plot.clf()
    plot.plot(nuNorm,apEffList,'ko')
    plot.plot(nuNorm,apEffInt(nuNorm),'k-')
    plot.title('Ap Eff')

    plot.figure(22)
    plot.clf()
    plot.plot(nuNorm,beamAreaList,'ko')
    plot.plot(nuNorm,beamAreaInt(nuNorm),'k-')
    plot.title('Ap Eff')
    plot.show()

def getairybeam(fileBeam,norm=True):
    
    from csv import reader
    import string
    from numpy import array,zeros
    from scipy.interpolate import interp2d

        
    ############################
    ## Read in beam file
    ############################
    fileBeam='../Inputs/airy_beams.dat' ##hard-coded
    fB=open(fileBeam,'r')
    fBlines=fB.readlines()
    nRad=len(fBlines)-1
    nWid=len(string.split(fBlines[0],','))-1
    fB.close()
    #nRad=len(fBeam)-1
    #nWid=len(fBeam[0])-1
    radList=zeros(nRad)
    widList=zeros(nWid)
    beamArr=zeros((nRad,nWid))
    
    row=string.split(fBlines[0][0:-1],',')
    for w in range(nWid):
        widList[w]=float(row[w+1])
    
    for r in range(nRad):
        row=string.split(fBlines[r+1],',')
        radList[r]=float(row[0])
        for w in range(nWid):
            beamArr[r,w]=float(row[w+1])
    
    beamInt=interp2d(radList,widList,beamArr,copy=True,kind='cubic')
    
    if norm:
        return(beamInt)
    else:
        return(beamInt)
        
def calcAiryTot():
    
    from numpy import pi,arange
    from scipy import where
    from scipy.special import j1

    dd=0.01    
    radListTot=arange(0.,300.,dd)
    radListTot=where(radListTot == 0.,1.e-12,radListTot)
    nRadTot=len(radListTot)
    airyList= ((2./pi) * (j1(pi*radListTot)/radListTot))**2
    airyTot=sum(airyList*2.*pi*radListTot*dd)
    
def makeAiryArr(maxRad,dd):
    
    from numpy import pi,sqrt,arctan2,zeros,zeros_like
    from scipy import where
    from scipy.special import j1
    
    pixCtr = maxRad / dd

    ##make position/radius Arrays
    npix = int(2.*maxRad/dd + 1)
    xArr = zeros((npix,npix))
    yArr = zeros_like(xArr)
    
    for x in range(npix):
        xArr[x,:] = -maxRad + x*dd
        yArr[:,x] = -maxRad + x*dd
        
    radArr = sqrt(xArr**2 + yArr**2)
    azArr = arctan2(yArr,xArr)
        
    ##remove zero
    radArr = where(radArr == 0.,1.e-12,radArr)
        
        ###radius is the angle in units of lambda/d

    ##Make Airy function
    bessArr = j1(pi*radArr)
    airyArr = ((2./pi) * (bessArr / radArr))**2
    
    return(xArr,yArr,airyArr)