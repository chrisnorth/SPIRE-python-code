# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 17:43:57 2012

@author: chris
"""

from csv import reader
import string
from numpy import array,pi
from beam import arcsec2sr
file='../Outputs/beamarea_theta_newBprofBS_Br700_Ind0.85.csv'

f=reader(open(file))

alpha=[]
area_arcsec=[]
area_sr=[]
for row in f:
    if string.find(row[0],'#')<0:
        alpha.append(float(row[0]))
        areas=array((float(row[1]),float(row[2]),float(row[3])))
        area_arcsec.append(areas)
        area_sr.append(arcsec2sr(areas)*1.e8)
        
nalph=len(alpha)

for a in range(nalph):
    line=r'%.1f & %d & %d & %d & %.4f & %.4f & %.4f \\'%(alpha[a],area_arcsec[a][0],area_arcsec[a][1],area_arcsec[a][2],area_sr[a][0],area_sr[a][1],area_sr[a][2])
    print line
