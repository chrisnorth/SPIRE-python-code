import string
from numpy import array
from beam import arcsec2sr
from csv import reader 

##Make latex tables for Observers' Manual

######################################################################
######   MAKE POWER LAW TABLES
######################################################################


file_pathin='../Outputs/'
file_pathout='../CalProducts/HIPE11/'
file_suff='_theta_final_newBprofBS_Br600_Ind0.85_ESA4'
filein_Kc=file_pathin+'KcolP'+file_suff+'.csv'
filein_Kc2=file_pathin+'KcolEuni'+file_suff+'.csv'
filein_area=file_pathin+'beamarea'+file_suff+'.csv'

fileout_KcolP=file_pathout+'KcolP_alpha_Cal.csv'
fileout_KcolE=file_pathout+'KcolE_alpha_Cal.csv'
fileout_area=file_pathout+'beamarea_alpha_Cal.csv'
fileout_arearel=file_pathout+'beamarea-rel_alpha_Cal.csv'

##read beam area
fareain=reader(open(filein_area))
alpha=[]
area_arcsec=[]
area_sr=[]
area_rel=[]
for row in fareain:
    if string.find(row[0],'#')<0:
        print row
        alpha.append(float(row[0]))
        areas=array((float(row[1]),float(row[2]),float(row[3])))
        area_arcsec.append(areas)
        area_sr.append(arcsec2sr(areas)*1.e8)
        if float(row[0])==-1:
            area_pip_arcsec=areas

for a in range(len(alpha)):
    area_rel.append(array([area_arcsec[a][0]/area_pip_arcsec[0],area_arcsec[a][1]/area_pip_arcsec[1],area_arcsec[a][2]/area_pip_arcsec[2]]))
    
fKcin=reader(open(filein_Kc))
KcPSW=[]
KcPMW=[]
KcPLW=[]
for row in fKcin:
    if string.find(row[0],'#')<0:
        print row
        KcPSW.append(float(row[1]))
        KcPMW.append(float(row[2]))
        KcPLW.append(float(row[3]))

KcPSW=array(KcPSW)
KcPMW=array(KcPMW)
KcPLW=array(KcPLW)

fKc2in=reader(open(filein_Kc2))
Kc2PSW=[]
Kc2PMW=[]
Kc2PLW=[]
for row in fKc2in:
    if string.find(row[0],'#')<0:
        print row
        Kc2PSW.append(float(row[1]))
        Kc2PMW.append(float(row[2]))
        Kc2PLW.append(float(row[3]))
        
KextPSW=array(Kc2PSW)
KextPMW=array(Kc2PMW)
KextPLW=array(Kc2PLW)

nalph=len(alpha)
print nalph
fKcPout=open(fileout_KcolP,'w')
fKcEout=open(fileout_KcolE,'w')
fareaout=open(fileout_area,'w')
farearelout=open(fileout_arearel,'w')


#write headers (csv files)
line1='#Colour correction factor (point source) with spectral index\n'
line2='#alpha, PSW, PMW, PLW\n'
line3='#Double, Double, Double, Double\n'
fKcPout.write(line1+line2+line3)
fKcEout.write(line1+line2+line3)
line1='#Effective beam solid angle with spectral index\n'
line2='#alpha, PSW ("), PMW ("), PLW (")\n'
line3='#Double, Double, Double, Double\n'
fareaout.write(line1+line2+line3)
line1='#Relative beam solid angle with spectral index\n'
line2='#alpha, PSW, PMW, PLW\n'
fareaout.write(line1+line2+line3)
for a in range(nalph):

    ##Kc
    line='%.1f , %.4f , %.4f , %.4f \n'%(alpha[a],KcPSW[a],KcPMW[a],KcPLW[a])
    fKcPout.write(line)
    line='%.1f , %.4f , %.4f , %.4f \n'%(alpha[a],KextPSW[a],KextPMW[a],KextPLW[a])
    fKcEout.write(line)
    
    #area file
    line='%.1f , %.4f , %.4f , %4f \n'%(alpha[a],area_arcsec[a][0],area_arcsec[a][1],area_arcsec[a][2])
    fareaout.write(line)    

    #area-rel file
    line='%.1f , %.4f , %.4f , %4f \n'%(alpha[a],area_rel[a][0],area_rel[a][1],area_rel[a][2])
    farearelout.write(line)    

fKcPout.close()
fKcEout.close()
fareaout.close()
farearelout.close()

######################################################################
######   MAKE GREY BODY TABLES
######################################################################

beta1=1.5
beta2=1.75
beta3=2.0
file_sufft1='_theta_newBprofBS_temp_beta%.1f_Br600_Ind0.85_ESA4'%(beta1)
file_sufft2='_theta_newBprofBS_temp_beta%.1f_Br600_Ind0.85_ESA4'%(beta2)
file_sufft3='_theta_newBprofBS_temp_beta%.1f_Br600_Ind0.85_ESA4'%(beta3)

filein_Kc_t1=file_pathin+'KcolP-t'+file_sufft1+'.csv'
filein_Kc2_t1=file_pathin+'KcolEuni-t'+file_sufft1+'.csv'
filein_area_t1=file_pathin+'beamarea'+file_sufft1+'.csv'

filein_Kc_t2=file_pathin+'KcolP-t'+file_sufft2+'.csv'
filein_Kc2_t2=file_pathin+'KcolEuni-t'+file_sufft2+'.csv'
filein_area_t2=file_pathin+'beamarea'+file_sufft2+'.csv'

filein_Kc_t3=file_pathin+'KcolP-t'+file_sufft3+'.csv'
filein_Kc2_t3=file_pathin+'KcolEuni-t'+file_sufft3+'.csv'
filein_area_t3=file_pathin+'beamarea'+file_sufft3+'.csv'

fileout_KcP_t=file_pathout+'KcolP_temp_Cal.csv'
fileout_KcE_t=file_pathout+'KcolE_temp_Cal.csv'
fileout_area_t=file_pathout+'beamarea_temp_Cal.csv'
fileout_arearel_t=file_pathout+'beamarea-rel_temp_Cal.csv'

##read beam area (beta1)
fareain_t1=reader(open(filein_area_t1))
temp=[]
area_arcsec_t1=[]
area_sr_t1=[]
area_rel_t1=[]
for row in fareain_t1:
    if string.find(row[0],'#')<0:
        #print row
        temp.append(float(row[0]))
        areas=array((float(row[1]),float(row[2]),float(row[3])))
        area_arcsec_t1.append(areas)
        area_sr_t1.append(arcsec2sr(areas)*1.e8)
        area_rel_t1.append(areas/area_pip_arcsec)
        
##read beam area (beta2)
fareain_t2=reader(open(filein_area_t2))
area_arcsec_t2=[]
area_sr_t2=[]
area_rel_t2=[]
for row in fareain_t2:
    if string.find(row[0],'#')<0:
        print row
        #temp.append(float(row[0]))
        areas=array((float(row[1]),float(row[2]),float(row[3])))
        area_arcsec_t2.append(areas)
        area_sr_t2.append(arcsec2sr(areas)*1.e8)
        area_rel_t2.append(areas/area_pip_arcsec)
        
##read beam area (beta3)
fareain_t3=reader(open(filein_area_t3))
area_arcsec_t3=[]
area_sr_t3=[]
area_rel_t3=[]
for row in fareain_t3:
    if string.find(row[0],'#')<0:
        print row
        #temp.append(float(row[0]))
        areas=array((float(row[1]),float(row[2]),float(row[3])))
        area_arcsec_t3.append(areas)
        area_sr_t3.append(arcsec2sr(areas)*1.e8)
        area_rel_t3.append(areas/area_pip_arcsec)
        

##read KcolP (beta1)
fKcin_t1=reader(open(filein_Kc_t1))
KcPSW_t1=[]
KcPMW_t1=[]
KcPLW_t1=[]
for row in fKcin_t1:
    if string.find(row[0],'#')<0:
        #print row
        KcPSW_t1.append(float(row[1]))
        KcPMW_t1.append(float(row[2]))
        KcPLW_t1.append(float(row[3]))

KcPSW_t1=array(KcPSW_t1)
KcPMW_t1=array(KcPMW_t1)
KcPLW_t1=array(KcPLW_t1)

##read KcolP (beta2)
fKcin_t2=reader(open(filein_Kc_t2))
KcPSW_t2=[]
KcPMW_t2=[]
KcPLW_t2=[]
for row in fKcin_t2:
    if string.find(row[0],'#')<0:
        #print row
        KcPSW_t2.append(float(row[1]))
        KcPMW_t2.append(float(row[2]))
        KcPLW_t2.append(float(row[3]))

KcPSW_t2=array(KcPSW_t2)
KcPMW_t2=array(KcPMW_t2)
KcPLW_t2=array(KcPLW_t2)

##read KcolP (beta3)
fKcin_t3=reader(open(filein_Kc_t3))
KcPSW_t3=[]
KcPMW_t3=[]
KcPLW_t3=[]
for row in fKcin_t3:
    if string.find(row[0],'#')<0:
        #print row
        KcPSW_t3.append(float(row[1]))
        KcPMW_t3.append(float(row[2]))
        KcPLW_t3.append(float(row[3]))

KcPSW_t3=array(KcPSW_t3)
KcPMW_t3=array(KcPMW_t3)
KcPLW_t3=array(KcPLW_t3)

##read KcolE (beta1)
fKc2in_t1=reader(open(filein_Kc2_t1))
Kc2PSW_t1=[]
Kc2PMW_t1=[]
Kc2PLW_t1=[]
for row in fKc2in_t1:
    if string.find(row[0],'#')<0:
        #print row
        Kc2PSW_t1.append(float(row[1]))
        Kc2PMW_t1.append(float(row[2]))
        Kc2PLW_t1.append(float(row[3]))

KextPSW_t1=array(Kc2PSW_t1)
KextPMW_t1=array(Kc2PMW_t1)
KextPLW_t1=array(Kc2PLW_t1)

##read KcolE (beta2)
fKc2in_t2=reader(open(filein_Kc2_t2))
Kc2PSW_t2=[]
Kc2PMW_t2=[]
Kc2PLW_t2=[]
for row in fKc2in_t2:
    if string.find(row[0],'#')<0:
        #print row
        Kc2PSW_t2.append(float(row[1]))
        Kc2PMW_t2.append(float(row[2]))
        Kc2PLW_t2.append(float(row[3]))
        
KextPSW_t2=array(Kc2PSW_t2)
KextPMW_t2=array(Kc2PMW_t2)
KextPLW_t2=array(Kc2PLW_t2)

##read KcolE (beta3)
fKc2in_t3=reader(open(filein_Kc2_t3))
Kc2PSW_t3=[]
Kc2PMW_t3=[]
Kc2PLW_t3=[]
for row in fKc2in_t3:
    if string.find(row[0],'#')<0:
        #print row
        Kc2PSW_t3.append(float(row[1]))
        Kc2PMW_t3.append(float(row[2]))
        Kc2PLW_t3.append(float(row[3]))
        
KextPSW_t3=array(Kc2PSW_t3)
KextPMW_t3=array(Kc2PMW_t3)
KextPLW_t3=array(Kc2PLW_t3)

ntemp=len(temp)
print ntemp
fKcPout_t=open(fileout_KcP_t,'w')
fKcEout_t=open(fileout_KcE_t,'w')
fareaout_t=open(fileout_area_t,'w')
farearelout_t=open(fileout_arearel_t,'w')

#write headers (csv files)
line1='#Colour correction parameter (point source) with temperature and emissivity index\n'
line2='#Beta: , %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n'%(beta1,beta1,beta1,beta2,beta2,beta2,beta3,beta3,beta3)
line3='#Temp (K), PSW, PMW, PLW, PSW, PMW, PLW, PSW, PMW, PLW\n'
line4='#Double, Double, Double, Double, Double, Double, Double, Double, Double, Double\n'
fKcPout_t.write(line1+line2+line3+line4)

line1='#Colour correction parameter (extended source source) with temperature and emissivity index\n'
fKcEout_t.write(line1+line2+line3+line4)

line1='#Effective beam area with temperature and emissivity index\n'
line3='#Temp (K), PSW ("), PMW ("), PLW ("), PSW ("), PMW ("), PLW ("), PSW ("), PMW ("), PLW (")\n'
fareaout_t.write(line1+line2+line3)

line1='#Relative beam area with temperature and emissivity index\n'
line3='#Temp (K), PSW, PMW, PLW, PSW, PMW, PLW, PSW, PMW, PLW\n'
farearelout_t.write(line1+line2+line3)


for t in range(ntemp):
    ## print every 1K for 3 K < T < 40 K
    if (temp[t] >= 3. and temp[t] <= 40. and temp[t]/int(temp[t])==1.):
        print temp[t]
        ##KcP file
        line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(temp[t],KcPSW_t1[t],KcPMW_t1[t],KcPLW_t1[t],KcPSW_t2[t],KcPMW_t2[t],KcPLW_t2[t],KcPSW_t3[t],KcPMW_t3[t],KcPLW_t3[t])
        fKcPout_t.write(line)
        
        ##KcE file
        line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(temp[t],KextPSW_t1[t],KextPMW_t1[t],KextPLW_t1[t],KextPSW_t2[t],KextPMW_t2[t],KextPLW_t2[t],KextPSW_t3[t],KextPMW_t3[t],KextPLW_t3[t])
        fKcEout_t.write(line)
        
        #area file
        line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(temp[t],area_arcsec_t1[t][0],area_arcsec_t1[t][1],area_arcsec_t1[t][2],area_arcsec_t2[t][0],area_arcsec_t2[t][1],area_arcsec_t2[t][2],area_arcsec_t3[t][0],area_arcsec_t3[t][1],area_arcsec_t3[t][2])
        fareaout_t.write(line)    
    
        #area-rel file
        line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(temp[t],area_rel_t1[t][0],area_rel_t1[t][1],area_rel_t1[t][2],area_rel_t2[t][0],area_rel_t2[t][1],area_rel_t2[t][2],area_rel_t3[t][0],area_rel_t3[t][1],area_rel_t3[t][2])
        farearelout_t.write(line)    

fKcPout_t.close()
fKcEout_t.close()
fareaout_t.close()
farearelout_t.close()
