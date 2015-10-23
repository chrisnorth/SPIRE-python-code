# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 16:08:29 2011

@author: chris
"""
from numpy import array,zeros,arange
from csv import reader
from matplotlib import pyplot as plot
def maincode():
    
    file1='../Outputs/Kc2b_theta_newBprof_Br500_Ind0.65.csv'
    file2='../Outputs/Kc2b_theta_newBprof_Br500_Ind1.00.csv'
    
    f1=reader(open(file1))
    f2=reader(open(file2))
    
    alist=[]
    f1_1=[]
    f1_2=[]
    f1_3=[]
    for row in f1:
        if 
        alist.append(float(row[0]))
        f1_1.append(float(row[1]))
        f1_2.append(float(row[1]))
        f1_3.append(float(row[1]))
    
    f2_1=[]
    f2_2=[]
    f2_3=[]
    for row in f2:
        f2_1.append(float(row[1]))
        f2_2.append(float(row[1]))
        f2_3.append(float(row[1]))
    
    alpharr=array(alist)
    
    res1=array([[f1_1],[f1_2],[f1_3]])
    res2=array([[f2_1],[f2_2],[f2_3]])
    
    print res1
    
if __name__=="__main__":
    maincode()