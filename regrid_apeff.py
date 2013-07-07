# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 11:03:08 2013

@author: chris
"""

from numpy import array,zeros
from csv import reader
from sys import exit
from scipy.interpolate import interp1d

filersrf='../Inputs/SPIRE-Phot-RSRF_nu.csv'
rsrf_input=reader(open(filersrf))
nursrf=[]
for row in rsrf_input:
    try:
        f=float(row[0])
        nursrf.append(float(row[0])) #GHz
    except:
        print row
nursrf=array(nursrf)

fileapef='../Inputs/app_eff.csv'
apef_input=reader(open(fileapef))
nuapef=[]
apeff=[]
for row in apef_input:
    nuapef.append(float(row[0])*1.e3) #THz->GHz
    apeff.append([float(row[1]),float(row[2]),float(row[3])])

apeff=array(apeff)
apeff_regrid=zeros((len(nursrf),4))
apeff_regrid[:,0]=nursrf
for b in range(3):
    intap=interp1d(nuapef,apeff[:,b],bounds_error=False,fill_value=0)
    apeff_regrid[:,b+1]=intap(nursrf)

fileout='../Inputs/SPIRE-Phot-ApEff_nu.csv'
fout=open(fileout,'w')
apeff_rows=[]
fout.write("Frequency,PSW,PMW,PLW\n")
fout.write("Double,Double,Double,Double\n")
fout.write("GHz,,,\n")
fout.write(",,,\n")

for r in range(len(apeff_regrid)):
    fout.write("%.12g,%.4g,%.4g,%.4g\n"%(apeff_regrid[r,0],apeff_regrid[r,1],apeff_regrid[r,2],apeff_regrid[r,3]))

fout.close()