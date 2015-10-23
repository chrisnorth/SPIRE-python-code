# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 13:50:28 2012

@author: chris
"""

###Plot beam map

import aplpy
#import pyfits
import matplotlib.pyplot as plot

beamdir='../Inputs/spire_beams_measured/'
bands=['psw','pmw','plw']

files=[]
for b in range(3):
    files.append('%s%s_beam_1arcsec.fits'%(beamdir,bands[b]))

for b in range(1):
    fig=aplpy.FITSFigure(files[b])    
    fig.show_grayscale()
    
plot.show()
