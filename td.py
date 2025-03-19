from __future__ import division
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import Path_check
import os
import pickle
import calendar
import numpy as np
# import airmass
from scipy.optimize import *
from scipy.stats import chi2
from PyAstronomy import pyasl

import airmass

folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)
[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
datafile_folder_omar = str(converted_Data_folder)+r'\dataset_omar\\'
fl_dataset_omar = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\Dataset_Omar\fits\*.fits')
converted_datafiles_omar = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\'+r'*.txt')
datafiles = []
for filepath in converted_datafiles_omar:
    a = open(filepath, 'rb')
    b = pickle.load(a)
    datafiles.append(b)
    a.close()
i=0
for file in datafiles:
    print(i, file.filename, file.header['DATE-OBS'])
    print(file.HJD)
    print('-------------')
    i+=1
# testfile=pf.open(fl_dataset_omar[0])
# header=testfile[0].header
# time = header['DATE-OBS']
# bjd1=airmass.bjd_lapalma_from_date_zet_ori(time)
# bjd2 =header['BJD']
# dif = (bjd1-bjd2)*24*60*60
# a=[1,2,3,4]
# b=str(a)
# print(b[0])
# testfile.close()
# print(bjd1)
# print(bjd2)
# print(dif)
# rng = np.random.default_rng()
#
# A = 2.  # amplitude
# c = 2.  # offset
#
# prot = 4
# w0 = 2*np.pi/prot  # rad/sec
# nin = 150
# nout = 1002
# x = rng.uniform(1, 20, nin)
# w = (2*np.pi)/np.linspace(1,10,1000)
# t= (2*np.pi)/w
# print(np.linspace(1,10,1000),t)
# y = A * np.cos(w0*x) + c
#
# # w = np.linspace(0.25, 10, nout)
# from scipy.signal import lombscargle
# pgram_power = lombscargle(x, y, w, normalize=False)
# plt.plot(t,pgram_power)
# plt.show()
# plt.close()
#
# a = np.linspace(1,10,10)
# b=1/a
# print(b)

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