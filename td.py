from __future__ import division
import matplotlib.pyplot as plt
import glob
import pyfits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
from scipy.optimize import *
from scipy.stats import chi2
from PyAstronomy import pyasl
filelist_lapalma = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
filelist2 =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
header = pf.open(filelist_lapalma[0])[0].header
asdf = airmass.fitfraunlp(filelist_lapalma[4])

# print header['RA'],header['DEC']
# def suck():
#     # header = pf.open(file)[0].header
#     JD = 2457306.585396
#     DEC = (-1.94)* (2*math.pi/360)
#     RA = 85.18* (2*math.pi/360)
#     LAT = 28.7636* (2*math.pi/360)
#     LON = -17.8947* (2*math.pi/360)
#     D = JD - 2451545.0
#     GMST = ((18.697374558 + 24.06570982441908*D)%24) * (360/24)* (2*math.pi/360)
#     sinalt = math.sin(DEC)*math.sin(LAT) + math.cos(DEC)*math.cos(LAT)*math.cos((GMST + LON - RA))
#     alt = np.arcsin(sinalt)*360/(2*math.pi)
#     airmass = 1/sinalt
#     return airmass,alt,JD
# am,alt,jd = suck()
# print alt