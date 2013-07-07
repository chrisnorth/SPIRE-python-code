# -*- coding: utf-8 -*-
"""
Created on Thu Sep 27 18:30:38 2012

@author: chris
"""

##test extended source calibration
from numpy import array
from csv import reader
import string

from beam import arcsec2sr

##calculate initial parameters
wl0=250.*1.e-6
I0=100.e6 #Source surface brightness in MJy/sr
alphaS=-3.

wl=array([250.,350.,500])*1.e-6

##calculate monochromatic surface brightness
I_S0=(wl0/wl)**alphaS * I0
print 'True Monochromatic Surface Brightness (I_S0,MJy/sr): %.5g , %.5g , %.5g'%(I_S0[0]/1e6,I_S0[1]/1e6,I_S0[2]/1e6)

##read beam areas in arcsec^2
file_area='../Outputs/beamarea_theta_newBprofBS_Br700_Ind0.85.csv'
farea=reader(open(file_area,'r'))
for row in farea:
    if string.find(row[0],'#')<0:
        if float(row[0]) == -1.:
            area0_arcsec=array((float(row[1]),float(row[2]),float(row[3])))
        if float(row[0]) == alphaS:
            area_arcsec=array((float(row[1]),float(row[2]),float(row[3])))

##set neptune beam area
areaNep_arcsec=array((450.63, 798.13, 1668.45))
##calculate beam area in sr
area_sr=arcsec2sr(area_arcsec)
area0_sr=arcsec2sr(area0_arcsec)
areaNep_sr=arcsec2sr(areaNep_arcsec)
print 'True beam areas (arcsec^2): %.5g , %.5g , %.5g'%(area_arcsec[0],area_arcsec[1],area_arcsec[2])
print 'True beam areas (sr): %.5g , %.5g , %.5g'%(area_sr[0],area_sr[1],area_sr[2])
print 'a=-1 beam areas (arcsec^2): %.5g , %.5g , %.5g'%(area0_arcsec[0],area0_arcsec[1],area0_arcsec[2])
print 'a=-1 beam areas (sr): %.5g , %.5g , %.5g'%(area0_sr[0],area0_sr[1],area0_sr[2])
area0_err=100.*(area0_arcsec - area_arcsec)/area_arcsec
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(area0_err[0],area0_err[1],area0_err[2])
print 'Neptune beam areas (arcsec^2): %.5g , %.5g , %.5g'%(areaNep_arcsec[0],areaNep_arcsec[1],areaNep_arcsec[2])
print 'Neptune beam areas (sr): %.5g , %.5g , %.5g'%(areaNep_sr[0],areaNep_sr[1],areaNep_sr[2])
areaNep_err=100.*(areaNep_arcsec - area_arcsec)/area_arcsec
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(areaNep_err[0],areaNep_err[1],areaNep_err[2])
##calculate monochromatic flux density
S_S0 = I_S0 * area_sr
print 'True Monochromatic Flux Density (S_S0,Jy/beam): %.5g , %.5g , %.5g'%(S_S0[0],S_S0[1],S_S0[2])

##read K_5
file_K5='../Outputs/K5_theta_newBprofBS_Br700_Ind0.85.csv'
fk5=reader(open(file_K5,'r'))
for row in fk5:
    if string.find(row[0],'#')<0:
        if float(row[0]) == alphaS:
            K_5=array((float(row[1]),float(row[2]),float(row[3])))

print 'Extended conversion (K_5): %.5g , %.5g , %.5g'%(K_5[0],K_5[1],K_5[2])

##read K_CE
file_Kc2='../Outputs/Kc2_theta_newBprofBS_Br700_Ind0.85.csv'
fkc2=reader(open(file_Kc2,'r'))
for row in fkc2:
    if string.find(row[0],'#')<0:
        if float(row[0]) == alphaS:
            K_CE=array((float(row[1]),float(row[2]),float(row[3])))
        if float(row[0]) == -1:
            K_CE0=array((float(row[1]),float(row[2]),float(row[3])))

##calculate measured flux density
S_S = I_S0 / K_5
print 'Broadband RSRF-weighted Flux Density (S_S,Jy/beam): %.5g , %.5g , %.5g'%(S_S[0],S_S[1],S_S[2])

##read K_CP
file_Kc='../Outputs/K4_theta_newBprofBS_Br700_Ind0.85.csv'
fkc=reader(open(file_Kc,'r'))
for row in fkc:
    if string.find(row[0],'#')<0:
        if float(row[0]) == -1:
            K4=array((float(row[1]),float(row[2]),float(row[3])))
print 'Pipeline conversion (K4): %.5g , %.5g , %.5g'%(K4[0],K4[1],K4[2])

Spip = S_S * K4
print 'Pipeline flux density (Spip): %.5g , %.5g , %.5g'%(Spip[0],Spip[1],Spip[2])
Spip_err = 100.*(Spip - S_S0)/S_S0
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(Spip_err[0],Spip_err[1],Spip_err[2])

##read K_CP
file_Kc='../Outputs/Kc_theta_newBprofBS_Br700_Ind0.85.csv'
fkc=reader(open(file_Kc,'r'))
for row in fkc:
    if string.find(row[0],'#')<0:
        if float(row[0]) == alphaS:
            K_CP=array((float(row[1]),float(row[2]),float(row[3])))
print 'Point source colour correction (K_CP): %.5g , %.5g , %.5g'%(K_CP[0],K_CP[1],K_CP[2])

##calculate point source flux density
S0_P = Spip * K_CP
print 'Point source approx flux density (S0_P/Jy/beam): %.5g , %.5g , %.5g'%(S0_P[0],S0_P[1],S0_P[2])
S0_P_err = 100.*(S0_P - S_S0)/S_S0
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(S0_P_err[0],S0_P_err[1],S0_P_err[2])

I0_P = S0_P / area_sr
print 'Point source approx surface brightness (I0_P,MJy/sr): %.5g , %.5g , %.5g'%(I0_P[0]/1e6,I0_P[1]/1e6,I0_P[2]/1e6)
I0_P_err = 100.*(I0_P - I_S0)/I_S0
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(I0_P_err[0],I0_P_err[1],I0_P_err[2])

I0_P_a0 = S0_P / area0_sr
print 'Point source approx surface brightness (I0_P_a0, a=-1, MJy/sr): %.5g , %.5g , %.5g'%(I0_P_a0[0]/1e6,I0_P_a0[1]/1e6,I0_P_a0[2]/1e6)
I0_P_a0_err = 100.*(I0_P_a0 - I_S0)/I_S0
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(I0_P_a0_err[0],I0_P_a0_err[1],I0_P_a0_err[2])

print 'Extended colour correction (K_CE, MJy/sr per Jy): %.5g , %.5g , %.5g'%(K_CE[0]/1e6,K_CE[1]/1e6,K_CE[2]/1e6)
I0_E = Spip * K_CE
print 'Extended surface brightness (I0_E, MJy/sr): %.5g , %.5g , %.5g'%(I0_E[0]/1e6,I0_E[1]/1e6,I0_E[2]/1e6)
I0_E_err = 100.*(I0_E - I_S0)/I_S0
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(I0_E_err[0],I0_E_err[1],I0_E_err[2])

print 'Extended colour correction (K_CE_0, MJy/sr per Jy,a=-1): %.5g , %.5g , %.5g'%(K_CE0[0]/1e6,K_CE0[1]/1e6,K_CE0[2]/1e6)
I0_E0 = Spip * K_CE0

print 'Extended surface brightness (I0_E, MJy/sr,a=-1): %.5g , %.5g , %.5g'%(I0_E0[0]/1e6,I0_E0[1]/1e6,I0_E0[2]/1e6)
I0_E0_err = 100.*(I0_E0 - I_S0)/I_S0
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(I0_E0_err[0],I0_E0_err[1],I0_E0_err[2])

S_Anep = I0_E * areaNep_sr
print 'Flux density with Neptune beam (S_Anep, Jy/beam): %.5g , %.5g , %.5g'%(S_Anep[0],S_Anep[1],S_Anep[2])
S_Anep_err = 100.*(S_Anep - S_S0)/S_S0
print '  Relative Error (%%): %.5g , %.5g , %.5g'%(S_Anep_err[0],S_Anep_err[1],S_Anep_err[2])