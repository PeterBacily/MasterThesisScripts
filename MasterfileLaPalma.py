from __future__ import division
import pyfits as pf
import glob
import airmass
# import datareduc
import matplotlib.pyplot as plt
import numpy as np

filelist = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')

file  =  pf.open(filelist[8])
days_P_er =  ((file[0].header['BJD'] - 2454391.)/6.829)* 0.001
er1 = days_P_er/6.829
er2 = (days_P_er+0.5)/6.829
print er1
print er2

print ' -------------------'


days_P_er =  ((2457473.333 - 2454391.)/6.829)* 0.001
er1 = days_P_er/6.829
er2 = (days_P_er+0.5)/6.829
print er1
print er2

