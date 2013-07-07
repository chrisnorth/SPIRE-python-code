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

fileout_Kc=file_pathout+'Kc_alpha_OM.csv'
fileout_area=file_pathout+'beamarea_alpha_OM.csv'
fileout_arearel=file_pathout+'beamarea-rel_alpha_OM.csv'
filetex_Kc=file_pathout+'Kc_alpha_OM.tex'
filetex_area=file_pathout+'beamarea_alpha_OM.tex'
filetex_arearel=file_pathout+'beamarea-rel_alpha_OM.tex'
filetwiki_Kc=file_pathout+'Kc_alpha_OM_twiki.dat'
filetwiki_arearel=file_pathout+'beamarea-rel_alpha_OM_twiki.dat'

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

#print area_rel
#print area_pip_arcsec 
#print shape(area_rel)

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
fKcout=open(fileout_Kc,'w')
fareaout=open(fileout_area,'w')
farearelout=open(fileout_arearel,'w')

fKctex=open(filetex_Kc,'w')
fareatex=open(filetex_area,'w')
fareareltex=open(filetex_arearel,'w')

fKctwiki=open(filetwiki_Kc,'w')
fareareltwiki=open(filetwiki_arearel,'w')


#write headers (csv files)
line='#alpha, PSW (psrc), PMW (psrc), PLW (psrc), PSW (extd), PMW (extd), PLW (extd)\n'
fKcout.write(line)
line='#alpha, PSW ("), PMW ("), PLW ("), PSW (sr), PMW (sr), PLW (sr)\n'
fareaout.write(line)
line='#alpha, PSW, PMW, PLW,\n'
farearelout.write(line)

#write headers (tex files)
line1='\\begin{tabular}{c|lll|lll|}\n'
line2=' & \\multicolumn{3}{c|}{Point Source (psrc)} & \\multicolumn{3}{c|}{Extended source (extd)} \\\\\n'
line3='$\\alpha$ & PSW & PMW & PLW & PSW & PMw & PLW \\\\\n'
line4='\\hline\n'
fKctex.write(line1+line2+line3+line4)

##change for area
line1='\\begin{tabular}{c|lll|lll|}\n'
line2=' & \\multicolumn{3}{c|}{Relative beam solid angle & \\multicolumn{3}{c|}{Relative beam solid angle \\\\\n'
fareatex.write(line1+line2+line3+line4)

##change for area-re;
line1='\\begin{tabular}{c|lll|}\n'
line2=' & \\multicolumn{3}{c|}{Relative beam solid angle \\\\\n'
fareareltex.write(line1+line2+line3+line4)

#write headers (twiki files)
line1='| | psrc | psrc | psrc | extd | extd | extd |\n'
line2='| alpha | PSW | PMW | PLW | PSW | PMW | PLW |\n'
fKctwiki.write(line1+line2)
line1='| alpha | PSW | PMW | PLW |\n'
fareareltwiki.write(line1)

for a in range(nalph):

    ##Kc
    line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(alpha[a],KcPSW[a],KcPMW[a],KcPLW[a],KextPSW[a],KextPMW[a],KextPLW[a])
    fKcout.write(line)
    line='%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n'%(alpha[a],KcPSW[a],KcPMW[a],KcPLW[a],KextPSW[a],KextPMW[a],KextPLW[a])
    fKctex.write(line)
    line='| %.1f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |\n'%(alpha[a],KcPSW[a],KcPMW[a],KcPLW[a],KextPSW[a],KextPMW[a],KextPLW[a])
    fKctwiki.write(line)
    
    #area file
    line='%.1f , %d , %d , %d , %.4f , %.4f , %.4f \n'%(alpha[a],area_arcsec[a][0],area_arcsec[a][1],area_arcsec[a][2],area_sr[a][0],area_sr[a][1],area_sr[a][2])
    fareaout.write(line)    
    line='%.1f , %.4f , %.4f , %.4f \n'%(alpha[a],area_rel[a][0],area_rel[a][1],area_rel[a][2])
    farearelout.write(line)    
    line='%.1f & %d & %d & %d & %.4f & %.4f & %.4f \\\\ \n'%(alpha[a],area_arcsec[a][0],area_arcsec[a][1],area_arcsec[a][2],area_sr[a][0],area_sr[a][1],area_sr[a][2])
    fareatex.write(line)
    line='%.1f & %.4f & %.4f & %.4f \\\\ \n'%(alpha[a],area_rel[a][0],area_rel[a][1],area_rel[a][2])
    fareareltex.write(line)
    line='| %.1f | %.4f | %.4f | %.4f |\n'%(alpha[a],area_rel[a][0],area_rel[a][1],area_rel[a][2])
    fareareltwiki.write(line)
    
line='\\end{tabular}\n'
fKctex.write(line)
fareatex.write(line)

fKcout.close()
fareaout.close()
farearelout.close()
fKctex.close()
fareatex.close()
fareareltex.close()
fKctwiki.close()
fareareltwiki.close()

######################################################################
######   MAKE GREY BODY TABLES
######################################################################

beta1=1.5
beta2=2.0
file_sufft1='_theta_newBprofBS_temp_beta%.1f_Br600_Ind0.85_ESA4'%(beta1)
file_sufft2='_theta_newBprofBS_temp_beta%.1f_Br600_Ind0.85_ESA4'%(beta2)

filein_Kc_t1=file_pathin+'KcolP-t'+file_sufft1+'.csv'
filein_Kc2_t1=file_pathin+'KcolEuni-t'+file_sufft1+'.csv'
filein_area_t1=file_pathin+'beamarea'+file_sufft1+'.csv'

filein_Kc_t2=file_pathin+'KcolP-t'+file_sufft2+'.csv'
filein_Kc2_t2=file_pathin+'KcolEuni-t'+file_sufft2+'.csv'
filein_area_t2=file_pathin+'beamarea'+file_sufft2+'.csv'

fileout_KcP_t=file_pathout+'KcP_temp_OM.csv'
fileout_KcE_t=file_pathout+'KcE_temp_OM.csv'
fileout_arearel_t=file_pathout+'beamarea-rel_temp_OM.csv'
fileout_area_t1=file_pathout+'beamarea_temp_beta%.1f_OM.csv'%beta1
fileout_area_t2=file_pathout+'beamarea_temp_beta%.1f_OM.csv'%beta2
filetex_KcP_t=file_pathout+'KcP_temp_OM.tex'
filetex_KcE_t=file_pathout+'KcE_temp_OM.tex'
filetex_area_t=file_pathout+'beamarea_temp_OM.tex'
filetex_arearel_t=file_pathout+'beamarea-rel_temp_OM.tex'
filetex_area_t1=file_pathout+'beamarea_temp_beta%.1f_OM.tex'%beta1
filetex_area_t2=file_pathout+'beamarea_temp_beta%.1f_OM.tex'%beta2

filetwiki_KcP_t=file_pathout+'KcP_temp_OM_twiki.dat'
filetwiki_KcE_t=file_pathout+'KcE_temp_OM_twiki.dat'
filetwiki_arearel_t=file_pathout+'beamarea-rel_temp_OM_twiki.dat'

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
        area_rel_t1.append(areas/area_pip_arcsec)
        area_sr_t1.append(arcsec2sr(areas)*1.e8)
        
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
        area_rel_t2.append(areas/area_pip_arcsec)
        area_sr_t2.append(arcsec2sr(areas)*1.e8)

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

ntemp=len(temp)
print ntemp
fKcPout_t=open(fileout_KcP_t,'w')
fKcEout_t=open(fileout_KcE_t,'w')
fareaout_t1=open(fileout_area_t1,'w')
fareaout_t2=open(fileout_area_t2,'w')
farearelout_t=open(fileout_arearel_t,'w')

fKcPtex_t=open(filetex_KcP_t,'w')
fKcEtex_t=open(filetex_KcE_t,'w')
fareareltex_t=open(filetex_arearel_t,'w')
fareatex_t1=open(filetex_area_t1,'w')
fareatex_t2=open(filetex_area_t2,'w')

fKcPtwiki_t=open(filetwiki_KcP_t,'w')
fKcEtwiki_t=open(filetwiki_KcE_t,'w')
fareareltwiki_t=open(filetwiki_arearel_t,'w')

#write headers (csv files)
line1='#Beta: , %.1f, %.1f, %.1f, %.1f, %.1f, %.1f\n'%(beta1,beta1,beta1,beta2,beta2,beta2)
line2='#Temp (K), PSW, PMW, PLW, PSW, PMW, PLW\n'
fKcPout_t.write(line1+line2)
fKcEout_t.write(line1+line2)
line2='#Temp (K), PSW ("), PMW ("), PLW ("), PSW (sr), PMW (sr), PLW (sr)\n'
fareaout_t1.write(line1+line2)
fareaout_t2.write(line1+line2)

#write headers (tex files)
line1='\\begin{tabular}{c|lll|lll|}\n'
line2=' & \\multicolumn{3}{c|}{Beam solid angle (arcsec$^2$)} & \\multicolumn{3}{c|}{Beam solid angle (sr)} \\\\\n'
line3='$\Temp (K)$ & PSW & PMW & PLW & PSW & PMw & PLW \\\\\n'
line4='\\hline\n'
fKcPtex_t.write(line1+line2+line3+line4)
fKcEtex_t.write(line1+line2+line3+line4)
line2=' & \\multicolumn{3}{c|}{Beam solid angle (arcsec$^2$)} & \\multicolumn{3}{c|}{Beam solid angle (sr)} \\\\\n'
fareatex_t1.write(line1+line2+line3+line4)
fareatex_t2.write(line1+line2+line3+line4)

line2=' & \\multicolumn{3}{c|}{$\beta = %.1f$} & \\multicolumn{3}{c|}{$\beta = %1.f$} \\\\\n'%(beta1,beta2)
line3='$\Temp (K)$ & PSW & PMW & PLW & PSW & PMw & PLW \\\\\n'
fareareltex_t.write(line1+line2+line3+line4)

#write headers (twiki files)
line1='| Beta > | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f |\n'%(beta1,beta1,beta1,beta2,beta2,beta2)
line2='| Temp (K) | PSW | PMW | PLW | PSW | PMW | PLW |\n'
fKcPtwiki_t.write(line1+line2)
fKcEtwiki_t.write(line1+line2)
fareareltwiki_t.write(line1+line2)

for t in range(ntemp):
    ## print every 1K for 3 K < T < 10 K and every 5 K for 10 K < T < 40 K    
    if (temp[t] >= 3.0 and temp[t] <=10. and temp[t]/int(temp[t])==1.) or (temp[t] >= 10.0 and temp[t] <=40. and (temp[t]/5.)/int(temp[t]/5.)==1.):
        print temp[t],area_arcsec_t1[t],area_arcsec_t2[t]
        ##KcP file
        line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(temp[t],KcPSW_t1[t],KcPMW_t1[t],KcPLW_t1[t],KcPSW_t2[t],KcPMW_t2[t],KcPLW_t2[t])
        fKcPout_t.write(line)
        line='%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n'%(temp[t],KcPSW_t1[t],KcPMW_t1[t],KcPLW_t1[t],KcPSW_t2[t],KcPMW_t2[t],KcPLW_t2[t])
        fKcPtex_t.write(line)
        line='| %.1f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |\n'%(temp[t],KcPSW_t1[t],KcPMW_t1[t],KcPLW_t1[t],KcPSW_t2[t],KcPMW_t2[t],KcPLW_t2[t])
        fKcPtwiki_t.write(line)
        
        ##KcE file
        line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(temp[t],KextPSW_t1[t],KextPMW_t1[t],KextPLW_t1[t],KextPSW_t2[t],KextPMW_t2[t],KextPLW_t2[t])
        fKcEout_t.write(line)
        line='%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n'%(temp[t],KextPSW_t1[t],KextPMW_t1[t],KextPLW_t1[t],KextPSW_t2[t],KextPMW_t2[t],KextPLW_t2[t])
        fKcEtex_t.write(line)
        line='| %.1f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |\n'%(temp[t],KextPSW_t1[t],KextPMW_t1[t],KextPLW_t1[t],KextPSW_t2[t],KextPMW_t2[t],KextPLW_t2[t])
        fKcEtwiki_t.write(line)
        
        #area file
        line='%.1f , %.4f , %.4f , %.4f , %.4f , %.4f , %.4f \n'%(temp[t],area_rel_t1[t][0],area_rel_t1[t][1],area_rel_t1[t][2],area_rel_t2[t][0],area_rel_t2[t][1],area_rel_t2[t][2])
        farearelout_t.write(line)    
        line='%.1f , %d , %d , %d , %.4f , %.4f , %.4f \n'%(temp[t],area_arcsec_t1[t][0],area_arcsec_t1[t][1],area_arcsec_t1[t][2],area_sr_t1[t][0],area_sr_t1[t][1],area_sr_t1[t][2])
        fareaout_t1.write(line)    
        line='%.1f , %d , %d , %d , %.4f , %.4f , %.4f \n'%(temp[t],area_arcsec_t2[t][0],area_arcsec_t2[t][1],area_arcsec_t2[t][2],area_sr_t2[t][0],area_sr_t2[t][1],area_sr_t2[t][2])
        fareaout_t2.write(line)    
        line='%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\\\ \n'%(temp[t],area_rel_t1[t][0],area_rel_t1[t][1],area_rel_t1[t][2],area_rel_t2[t][0],area_rel_t2[t][1],area_rel_t2[t][2])
        fareareltex_t.write(line)
        line='%.1f & %d & %d & %d & %.4f & %.4f & %.4f \\\\ \n'%(temp[t],area_arcsec_t1[t][0],area_arcsec_t1[t][1],area_arcsec_t1[t][2],area_sr_t1[t][0],area_sr_t1[t][1],area_sr_t1[t][2])
        fareatex_t1.write(line)
        line='%.1f & %d & %d & %d & %.4f & %.4f & %.4f \\\\ \n'%(temp[t],area_arcsec_t2[t][0],area_arcsec_t2[t][1],area_arcsec_t2[t][2],area_sr_t2[t][0],area_sr_t2[t][1],area_sr_t2[t][2])
        fareatex_t2.write(line)

        line='| %.1f | %.4f | %.4f | %.4f | %.4f | %.4f | %.4f |\n'%(temp[t],area_rel_t1[t][0],area_rel_t1[t][1],area_rel_t1[t][2],area_rel_t2[t][0],area_rel_t2[t][1],area_rel_t2[t][2])
        fareareltwiki_t.write(line)

line='\\end{tabular}\n'
fKcPtex_t.write(line)
fKcEtex_t.write(line)
fareatex_t1.write(line)
fareatex_t2.write(line)

fKcPout_t.close()
fKcEout_t.close()
fareaout_t1.close()
farearelout_t.close()
fareaout_t2.close()
fKcPtex_t.close()
fKcEtex_t.close()
fareareltex_t.close()
fareatex_t1.close()
fareatex_t2.close()
fKcPtwiki_t.close()
fKcEtwiki_t.close()
fareareltwiki_t.close()
