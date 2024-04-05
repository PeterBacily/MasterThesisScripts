# __author__ = 'PeterBacily'
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

c_light = 299792.458

filelist_lapalma = glob.glob('D:\Peter\Master Thesis\Data\LaPalmaData' + '\*.fits')
def bjd(file):
    header = pf.open(file)[0].header
    HJD =  float(header['BJD'])
    return HJD


def aphase(filetime):
    # datafile = pf.open(file)
    period = 6.829
    jdstart = Time(2454391., format='jd').jd
    # header =  datafile[0].header
    # filetime = Time(header['MJD-OBS'], format='jd').jd
    phase = ((filetime- jdstart)% period )/ period
    return phase

# file = pf.open(filelist_lapalma[0])
# header = file[0].header
# for i, item in enumerate(header):
#     print i, item
# startdate = header['DATE-OBS']
# HJD =  header['BJD']
# exptime = header['EXPTIME']
# SNR = header['SNR50']
# ALT = header['TELALT']
# AM = 1/np.sin(ALT)
# phase = aphase(HJD)
# BC = header['BVCOR']
# fwl = airmass.fitfraun(file)
# velshift = c_light*(fwl-5895.92)/5895.92
# print header['TELALT']
# print header['DATE-OBS'],header['DATE-END'],header['DATE-AVG'],header['EXPTIME']


newlist = sorted(filelist_lapalma, key=lambda x: bjd(x), reverse=False)
print r'\begin{tabular}{ l|| r| r| r| r|r|r|r|r }'
print r'# & Start date & Mid HJD& T_{\textrm{exp}}  & SNR & airmass & Alt & phase & v_{\baryc}& $ v_{IS}$\\'
print r' & 2015 & + 2457000 & (s)& & & (deg)& p = 6.829 d& (km/s) & (km/s)\\'

for i,file in enumerate(newlist):
    header = pf.open(file)[0].header
    startdate = airmass.timeanddatelp(file)
    HJD =  float(header['BJD'])
    exptime = header['EXPTIME']
    SNR = header['SNR50']
    ALT = float(header['TELALT'])
    AM = 1/np.sin(2*np.pi*ALT/360)
    phase = round(aphase(HJD),3)
    BC = float(header['BVCOR'])
    fwl = airmass.fitfraunlp(file)
    velshift = c_light*(fwl-5895.92)/5895.92
    print str(i), '&' ,str(startdate),'&' ,str(round(HJD-2457000,3)), '&' ,str(int(round(exptime))), '&' ,str(int(np.round(SNR))), '&' ,str(np.round(AM,1)),'&' ,str(np.round(ALT,1)), '&' ,str(phase),'&' ,str(int(np.round(BC))),'&', str(int(np.round(velshift))),r'\\'
print r'\end{tabular}'