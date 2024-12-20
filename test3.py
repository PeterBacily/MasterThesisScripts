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
import Path_check
import os
import pickle

a=[1,2,3]
print(a[0:])
print(type(a[0:-2]),type(a[-1]))
# a =[[1,1,1],[10,10,10],[100,100,100]]
# print(np.average(a,axis=0))
# b=[a[0],a[1]]
exit()
# a=6
# b=[7,9]
#
# print(np.min(a),np.max(a))
# print(np.min(b),np.max(b))
# import seaborn as sns
folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)
#
[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
data_full_night_all= str(converted_Data_folder)+r'\demetra\with_orders\full_night\\'
data_merc = str(converted_Data_folder)+r'\mercator\ll_apo_vcor_2\\'
data_audela = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\AudeLA\all\\'
filelist_merc=open_masterfiles.mercator(data_merc)
filelist_apo=open_masterfiles.apo_demetra_orders(data_full_night_all)
filelist_audela = open_masterfiles.apo(data_audela)
data_individual = str(converted_Data_folder)+r'\demetra\with_orders\Individual\\'
# # for file in filelist_merc:
# #     print(file.header['CDELT1'])
# #     nf_ha = file.line6562.normalizationflux
# #     snr_ha = 1 / np.std(nf_ha)
# #     print(snr_ha)
# # print('------')
# # for file in filelist_apo:
# #     wlarr= file.line6562_order.wl
# #     wldist = []
# #     for i in range (len(wlarr)-5):
# #         wld = wlarr[i+1]-wlarr[i]
# #         wldist.append(wld)
# #     print(wldist)
data_individual_list = open_masterfiles.apo_demetra_orders(path = data_individual,manual_filelist=None,sort_data_files='on')
data_full_night_all_list = open_masterfiles.apo_demetra_orders(path = data_full_night_all,manual_filelist=None,sort_data_files='on')
datafile_merc=filelist_merc[0]
wl=datafile_merc.wl_rebin2
flux=datafile_merc.flux_rebin2
plt.plot(wl,flux)
# plt.show()
# plt.close()
wlpiece = [5335, 5345]
order=airmass.find_order(wlpiece,data_individual_list[0])
demwl, demflux=order.wl_rebin, order.flux_rebin
# plt.plot(demwl,demflux)
plt.show()
plt.close()
# i=0
# orderfile_comp=[data_individual_list[3],data_full_night_all_list[7]]
# # print(orderfile_comp)
# # for file in orderfile_comp:
# #     wave = file.line6562_order_rebin.wl
# #     flux = file.line6562_order_rebin.flux
# #     plt.plot(wave,flux+i)
# #     i += 0.1
# # mf = filelist_merc[6]
# # wave = mf.line6562_rebin.wl
# # flux = mf.line6562_rebin.flux
# # plt.plot(wave,flux+i)
# # plt.show()
# # plt.close()
#
# from collections import defaultdict
#
#
# groups = defaultdict(list)
#
# for obj in data_individual_list:
#     # print(obj.time_and_date)
#     groups[obj.time_and_date[0:5]].append(obj)
#
# new_list = groups.values()
# # day=new_list[0]
# # for day in new_list:
# full_data = []
# for day in list(new_list):
#     day_data = []
#     for i in range(18):
#         wl_rebin = np.array(day[0].orders[i].wl_original)[10:-10]
#         day_order_data = []
#         for k in range(len(day)):
#             wl1 = np.array(day[k].orders[i].wl_original)
#             flux1=np.array(day[k].orders[i].flux_original)
#             flux_rebin = airmass.rebin_spec(wl1,flux1,wl_rebin)
#             snr_ha = airmass.snr_ha(day[k], return_only_snr=True)
#             # print(snr_ha)
#             day_order_data.append([wl1,flux1,wl_rebin,flux_rebin,snr_ha])
#             # plt.plot(wl_rebin,flux_rebin)
#
#             # print('wlarraystepmin=',min(np.diff(wl1)))
#             # print(len(wl1))
#             # print(eq)
#             # if not np.allclose(wl1,wl2):
#             #     print('not close')
#             # print('------------')
#             # print(wl1)
#             # print(wl2)
#             # print('------------')
#             # flux =observation.orders[i].flux_original
#         # plt.show()
#         # plt.close()
#         day_data.append(day_order_data)
#     full_data.append(day_data)
#
# testday = full_data[0]
# testorder = testday[12]
# for spec in testorder:
#     plt.plot(spec[2],spec[3])
#     print(len(spec[3]))
#
# plt.show()
# plt.close()
    # print('cd',file.header['CDELT1'])
    # nf_ha = file.line6562.normalizationflux
    # snr_ha = 1 / np.std(nf_ha)
    # print(snr_ha)
# print('-----')
# for file in filelist_audela:
#     # print(file.header['CDELT1'])
#     nf_ha = file.line6562.normalizationflux
#     snr_ha = 1 / np.std(nf_ha)
#     print(snr_ha)
#
# file_apo=filelist_audela[7]
# file_merc=filelist_merc[6]
# wl_apo,flux_apo = file_apo.wl_original,file_apo.flux_original
# [a,b,c,d]=file_apo.line6562.normalization_boundaries_wl
# normwave_apo = np.hstack((wl_apo[(wl_apo > a) & (wl_apo < b)], wl_apo[(wl_apo > c) & (wl_apo < d)]))
# normflux_apo = np.hstack((flux_apo[(wl_apo > a) & (wl_apo < b)], flux_apo[(wl_apo > c) & (wl_apo < d)]))
# flux_left = flux_apo[(wl_apo > a) & (wl_apo < b)]
# flux_right =  flux_apo[(wl_apo > c) & (wl_apo < d)]
# n_fl = flux_left/np.average(flux_left)
# n_fr = flux_right/np.average(flux_right)
# normflux_apo_2  = np.hstack((n_fl,n_fr))
# snr2 = 1/np.std(normflux_apo_2)
# slope,height = np.polyfit(normwave_apo,normflux_apo,1)
#
#     # print 'slope and height are', slope, height
# fit = np.poly1d([slope,height])
# print(np.poly1d(fit))
#
#
# nnf=[]
# for k, nwl in enumerate(normwave_apo):
#     nnf.append(normflux_apo[k] / fit(nwl))
# # plt.plot(normflux_apo)
# # plt.plot(nnf)
# # plt.show()
# # plt.close()
# print(nnf)
# print(np.std(nnf))
# snr = 1/np.std(nnf)
# print('snr=',snr)
# print('snr2=',snr2)
# print(a,b,c,d)
# # print(normwave_apo,normflux_apo)
# # wl_merc,flux_merc = file_merc.wl_original,file_merc.flux_original
# plt.plot(wl_apo,flux_apo)
#
# # plt.plot(wl_merc,flux_merc)
# plt.show()
# plt.close()
# file = r"D:\peter\Master_Thesis\Datareduction\Data\LaPalmaData\zet Ori2220151012.fits"
# data = pf.open(file)
# apo_file=open_masterfiles.apo_demetra_orders(data_full_night_all)[5]
# apo_lines = ['line6562_order']
# linedata = getattr(apo_file, 'line6562_order')
# flux = linedata.flux
# wave = linedata.wl
# [a,b,c,d]=linedata.normalization_boundaries_wl
# normwave = np.hstack((wave[(wave > a) & (wave < b)], wave[(wave > c) & (wave < d)]))
# normflux = np.hstack((flux[(wave > a) & (wave < b)], flux[(wave > c) & (wave < d)]))
# print(normflux)
# nf= linedata.normalizationflux
# print(1/np.std(nf))

# hd = apo_files.header
# date = apo_files.header['JD-MID']
# print(date)
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