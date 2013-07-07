# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:43:57 2012

@author: chris
"""

from csv import reader
import string
from numpy import array

betas=[1.5]
for beta in betas:
    filepip='../Outputs/Kc2_theta_final_newBprofBS_Br600_Ind0.85_ESA4.csv'
    file1 = '../Outputs/Kct_theta_newBprofBS_temp_beta%.1f_Br600_Ind0.85_ESA4.csv'%(beta)
    file2 = '../Outputs/Kct2_theta_newBprofBS_temp_beta%.1f_Br600_Ind0.85_ESA4.csv'%(beta)

    fpip=reader(open(filepip))
    f1=reader(open(file1))
    f2=reader(open(file2))

    for row in fpip:
        if string.find(row[0],'#')<0:
            if float(row[0]) == -1:
                kpipE=array((float(row[1]),float(row[2]),float(row[3])))
    
    temp=[]
    temp2=[]
    kct=[]
    kct2=[]
    for row in f1:
        if string.find(row[0],'#')<0:
            temp.append(float(row[0]))
            kct_t=array((float(row[1]),float(row[2]),float(row[3])))
            kct.append(kct_t)
    ntemp=len(temp)
        
    for row in f2:
        if string.find(row[0],'#')<0:
            temp2.append(float(row[0]))
            kct2_t=array((float(row[1]),float(row[2]),float(row[3])))
            kct2.append(kct2_t/kpipE)
    if len(temp2) != ntemp:
        print 'different numbers of temperatures!'
            
    for t in range(ntemp):
        if ((int(temp[t]) == temp[t]) and (temp[t] < 10)):
            #only print every 1K for T<=10
            line=r'%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\'%(temp[t],kct[t][0],kct[t][1],kct[t][2],kct2[t][0],kct2[t][1],kct2[t][2])
            print line
        if int(temp[t]/5.) == temp[t]/5. and temp[t] >= 10:
            #only print every 5K for T>10
            line=r'%.1f & %.4f & %.4f & %.4f & %.4f & %.4f & %.4f \\'%(temp[t],kct[t][0],kct[t][1],kct[t][2],kct2[t][0],kct2[t][1],kct2[t][2])
            print line