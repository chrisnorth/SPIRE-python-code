# -*- coding: utf-8 -*-
"""
Created on Tue Jul 17 22:40:52 2012

@author: chris
"""

fileBeam='../Inputs/airy_beams.dat'
fBeam=open(fileBeam,'w')
line='#Rad'
for w in range(nWid):
    line='%s,%.3f'%(line,wid[w])
line=line+'\n'
fBeam.write(line)

for r in range(nRad):
    line='%.3f'%(radList[r])
    for w in range(nWid):
        line='%s,%.5g'%(line,convArr[r,w])
    line=line+'\n'
    fBeam.write(line)
fBeam.close()