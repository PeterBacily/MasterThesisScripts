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
import specutils
import pickle


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
data_full_night_all= str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin02\combined\high_snr\\'
data_merc = str(converted_Data_folder)+r'\mercator\ll_apo_vcor_2\\'
data_audela = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\AudeLA\all\\'
filelist_merc=open_masterfiles.mercator(data_merc)
filelist_apo=open_masterfiles.apo_demetra_orders(data_full_night_all)
filelist_audela = open_masterfiles.apo(data_audela)
data_individual = str(converted_Data_folder)+r'\demetra\with_orders\Individual\\'

binsize = '01'
di = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin' + binsize + r'\single_obs\\'
fn = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin' + binsize + r'\combined\high_snr\\'
# data_individual_list = open_masterfiles.apo_demetra_orders(path=di, manual_filelist=None, sort_data_files='on')
data_full_night_all_list = open_masterfiles.apo_demetra_orders(path=fn, manual_filelist=None, sort_data_files='on')
tf = data_full_night_all_list[0]


ha_hb_linelist = ['line6562','line4861']

def degrade_spectrum(wl,flux,spectral_resolution=10000, desired_snr=110):
    deg_wl = wl
    noise_array=[]
    deg_flux, fwhm = pyasl.instrBroadGaussFast(wl, flux, spectral_resolution,
                                         edgeHandling="firstlast", fullout=True, maxsig=5.0)

    for f in deg_flux:
        noise_exp = f/desired_snr
        noise = np.random.normal(0, noise_exp)
        noise_array.append(noise)
    deg_flux_with_noise = np.array(deg_flux)+np.array(noise_array)
    return deg_wl,deg_flux_with_noise

def degrade_spectrum2(wl,flux,spectral_resolution=10000, desired_snr=70):
    deg_wl = wl
    noise_array=[]




    for f in flux:
        noise_exp = f/desired_snr

    #     # print(f, noise_exp)
    #     noise= np.random.poisson(lam=noise_exp, size=None)
        noise = np.random.normal(0, noise_exp)
        noise_array.append(noise)
    flux_with_noise = np.array(flux)+np.array(noise_array)
    # test_array = np.column_stack((deg_flux, deg_flux_with_noise,noise_array))
    deg_flux_with_noise, fwhm = pyasl.instrBroadGaussFast(wl, flux_with_noise, spectral_resolution,
                                         edgeHandling="firstlast", fullout=True, maxsig=5.0)
    return deg_wl,deg_flux_with_noise



def slice_and_norm(wl,flux,start,end,rebin=None):
    slice_flux = flux[(wl > start) & (wl < end)]
    slice_wl = wl[(wl > start) & (wl < end)]
    slice_flux_norm = slice_flux / np.average(slice_flux)
    if rebin == None:
        return slice_wl,slice_flux_norm
    else:
        slice_wl_rebinned,slice_flux_rebinned = airmass.rebin2(slice_wl,slice_flux_norm,step=rebin)
        return slice_wl_rebinned,slice_flux_rebinned

dark_path = r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\Zet_Ori_Data_Altair_Response\20160228\20160304-210922-DARK-600s-1.fit'
flat_path = r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\Zet_Ori_Data_Altair_Response\20160228\20160228-224418-FLAT-85s-1.fit'

# flat_file = pf.open(dark_path)
# print(flat_file.info())
# flatdata = flat_file[0].data
# plt.imshow(flatdata,clim=(185, 230))
# # plt.hist(flatdata.ravel(), bins=range(256), fc='k', ec='k')
# plt.show()
# plt.close()
# exit()

bd = [5170, 5210]
# bd =[6549, 6550.7, 6576.0, 6578.0]
merc_file= filelist_merc[0]
merc_wl=merc_file.wl_rebin
merc_flux=merc_file.flux_rebin
list_of_orders = tf.orders
order_file = airmass.find_order(bd,tf)
print(order_file.order_number_demetra)

apo_wl = order_file.wl_rebin
apo_flux = order_file.flux_rebin
apo_wl_slice,apo_flux_slice = slice_and_norm(apo_wl,apo_flux,bd[0],bd[1])
# print(airmass.velocity_to_wl([-1174,-985],6562.819))
# wl_apo = tf.line4861_order_rebin.wl
# flux_apo = tf.line4861_order_rebin.flux
# start_wl = wl_apo[0]
# end_wl= wl_apo[-1]

# print(merc_wl[4]-merc_wl[3],merc_wl[-9]-merc_wl[-10])

# exit()
sr = 10000
d_snr = 35
deg_wl,deg_flux = airmass.degrade_spectrum_noise_first(merc_wl,merc_flux,spectral_resolution=sr, desired_snr=d_snr,pre_rebin = 0.1)
# deg_wl,deg_flux = airmass.degrade_spectrum(merc_wl,merc_flux,spectral_resolution=10000, desired_snr=120,pre_rebin = 0.05)
# deg_wl_rebin,deg_flux_rebin = airmass.rebin2(deg_wl,deg_flux,0.1)
# wl_apo_sn,flux_apo_sn = slice_and_norm(wl_apo,flux_apo,start_wl,end_wl,rebin=None)
# wl_merc_sn,flux_merc_sn = slice_and_norm(merc_wl,merc_flux,start_wl,end_wl,rebin=None)
md_snr_ha,md_snr_straight = airmass.SNR_merc_degen(deg_wl,deg_flux)
apo_snr_ha,apo_snr_straight = airmass.SNR_apo_orders(tf)
print(md_snr_straight,apo_snr_straight)
wl_deg_sn, flux_deg_sn = slice_and_norm(deg_wl,deg_flux,bd[0],bd[1],rebin=None)
plt.plot(wl_deg_sn,flux_deg_sn,label='Mercator degenerated, R='+str(sr)+'\nSNR given='+str(d_snr)+' SNR measured='+str(np.round(md_snr_straight)))
plt.plot(apo_wl_slice,apo_flux_slice,label='APO SNR measured='+str(np.round(apo_snr_straight)))
plt.title('Continuum slice of APO and degenrated mercator spectra\n rebinned to 0.1 Å')
plt.xlabel('Wavelength (Å)')
plt.ylabel('Relative flux')
plt.legend()
plt.show()
plt.close()

exit()

snr_merc_original = airmass.snr_2(merc_wl,merc_flux,boundaries=bd,rebin=False,rebin_size=0.1,separate=False)
snr_merc_degrade_no_rebin = airmass.snr_2(deg_wl,deg_flux,boundaries=bd,rebin=False,rebin_size=0.1,separate=False)


print('orig',snr_merc_original)
print('deg, no rebin',snr_merc_degrade_no_rebin)
plt.plot(merc_wl,merc_flux)
plt.plot(deg_wl,deg_flux)
plt.legend()
plt.show()
plt.close()
# ha_order = list(filter(lambda x: x.order_number_demetra == '34', list_of_orders))[0]
exit()
# wl=ha_order.wl_rebin
# flux = ha_order.flux_rebin
plt.plot(wl_apo_sn,flux_apo_sn,label='apo')
# plt.plot(wl_merc_sn,flux_merc_sn,label = 'merc')
plt.plot(wl_deg_sn,flux_deg_sn,label = 'naughty,bad wl')
plt.legend()
plt.show()
plt.close()


# datafile_merc=filelist_merc[0]
# wl=datafile_merc.wl_original
# flux=datafile_merc.flux_original
# binsizes = [j-i for i, j in zip(wl[:-1], wl[1:])]
# print(binsizes)
# print(np.average(binsizes[3:-2]))
# bs_ratio = []
# for i in range(len(binsizes)):
#     bs_ratio.append(wl[i]/binsizes[i])
# print(sorted(bs_ratio)[:3],sorted((bs_ratio)[-3:]))
# print(np.average(bs_ratio[3:-2]))
# axes = {}
# fig = plt.figure()
# axes['ax0'] = fig.add_subplot(212)
# for i in range(3):
#     subplotnum = 330+i+1
#     axes[f'ax{i+1}'] = fig.add_subplot(subplotnum)
#
#
#
# plt.show()
# plt.close()
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

# def returnfunction(x):
#     return x,x**2
# class HolyGrail:
#
#     def __init__(self,steps = [1,2,3]):
#         self.start = 'start_at_init'
#         setattr(self,'a',returnfunction(2)[0])
#         setattr(self, 'b', returnfunction(2)[1])
#     # function definition in question:
#     # TypeError: 'str' object is not callable
#
# testobj = HolyGrail()
# print(testobj.a,testobj.b)

def test_binsize_raw():
    binsize='05'
    di = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin' + binsize + r'\single_obs\\'
    fn = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin' + binsize + r'\combined\\'
    data_individual_list = open_masterfiles.apo_demetra_orders(path=di, manual_filelist=None, sort_data_files='on')
    data_full_night_all_list = open_masterfiles.apo_demetra_orders(path=fn, manual_filelist=None, sort_data_files='on')
    tf = data_full_night_all_list[0]
    linedata = tf.line6562_order
    wl = linedata.wl
    binsizes = [j-i for i, j in zip(wl[:-1], wl[1:])]
    print(binsizes)
    print(np.average(binsizes))
exit()
rebinsizes = ['025','01','02','05']
for binsize in rebinsizes:
    di=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin'+binsize+r'\single_obs\\'
    fn =r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin'+binsize+r'\combined\\'
    data_individual_list = open_masterfiles.apo_demetra_orders(path = di,manual_filelist=None,sort_data_files='on')
    data_full_night_all_list = open_masterfiles.apo_demetra_orders(path = fn,manual_filelist=None,sort_data_files='on')
    # # datafile_merc=filelist_merc[0]
    tf = data_full_night_all_list[0]
    # tf = data_individual_list[0]
    linedata = tf.line6562_order_rebin
    wl=linedata.wl
    flux=linedata.flux
    print(tf.mark)
    plt.plot(wl,flux)
    plt.show()
    plt.close()
# # wl=datafile_merc.wl_rebin2
# # flux=datafile_merc.flux_rebin2
# # wlpiece = [5335, 5345]
# order=airmass.find_order(wlpiece,data_individual_list[0])
# demwl, demflux=order.wl_rebin, order.flux_rebin
# plt.plot(demwl,demflux)
# plt.show()
# plt.close()
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