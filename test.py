# __author__ = 'PeterBacily'
from __future__ import division
import warnings
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
import open_masterfiles
import scipy.stats as ss
from scipy.optimize import *
from PyAstronomy import pyasl

# file = r"D:\peter\Master_Thesis\Datareduction\Data\LaPalmaData\zet Ori2220151012.fits"
# data = pf.open(file)
apo_files=open_masterfiles.apo_demetra_orders()[3]
hd = apo_files.header
date = apo_files.header['JD-MID']
print(date)
# norm_wls = airmass.velocity_to_wl([404,455],5875.621)
# print(norm_wls)
# ll_lapalma = [[r'Ha', 6562.819, 6551, 6552, 6578, 6579], [r'Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hy', 4340.472, 4322, 4324, 4357, 4360], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8], ['He_II', 4541.6, 4498, 4499, 4580, 4581], ['He_II', 4685.804, 4679, 4680, 4690, 4691], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1]]
#
#
# def my_sin(x,  amplitude, phase, offset):
#     return amplitude*np.sin(2*np.pi*x  + phase)  + offset
#
# def my_line(x,a,b):
#     return a*x+b
#
# def flat_line(x,a):
#     return a
#
# def chisq(exp,obs):
#     obs1 = np.array(obs)
#     exp1 = np.array(exp)
#     sigma = np.avg(obs)
#
# filelist = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
# file1 = filelist[0]
# # datafolder = r'C:\peter\School\Master Scriptie\Data\EW\lp/'
# # datafile_list = glob.glob(datafolder + '*.npy')
# # file1 =np.load(datafile_list[0])
# # ew = file1[0]
# # print ew
# # phase = file1[1]
# # error = file1[2]
# # lineused = file1[3][0]
# # # print lineused
# # linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
# # p1=[0.005,0.5, 1]
# # fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
# # p2 = [0,1]
# # fit2 = curve_fit(my_line, phase, ew, p0=p2,  sigma=error, absolute_sigma=True)
# # p3 =[1]
# # fit3 = curve_fit(flat_line , phase, ew, p0=p3,  sigma=error, absolute_sigma=True)
# # # t= np.linspace(0,1,50)
# # data_fit_sin = my_sin(np.array(phase), *fit1[0])
# # data_fit_line = my_line(np.array(phase), *fit2[0])
# # data_fit_flat = flat_line(np.array(phase), *fit3[0])
# #
# # chi2_sin,p_sin = ss.chisquare(ew,data_fit_sin,ddof=3)
# # chi2_line,p_line = ss.chisquare(ew,data_fit_line,ddof=2)
# # chi2_flat,p_flat = ss.chisquare(ew,data_fit_flat,ddof=1)
#
#
# N = result.nobs
# SSR = result.ssr
# s2 = SSR / N
# L = ( 1.0/np.sqrt(2*np.pi*s2) ) ** N * np.exp( -SSR/(s2*2.0) )
# print 'ln(L) =', np.log( L )