from __future__ import division
import matplotlib.pyplot as plt
import glob
import pyfits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
import scipy.stats as ss
from scipy.optimize import *
from PyAstronomy import pyasl
import itertools


c_light = 299792.458
# ['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],
linelist_apo = [ ['H_beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I', 4713.1457, 4707.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2],  ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
linelist = [4861.01,4471.47,6678.15,4921.93,4647.42,5015.68,5411.52,4650.85,4541.52,5592.25]
ll2= [['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0],['H_gamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0],['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6],['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0],['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914],['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471],['He_II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],['C_IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
filelist_lapalma = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
filelist_apo =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
nall = [['Na doublet',5892.95,5888.5,5889,5894.2,5894.75]]
ll3 = [[r'H$\beta$', 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I', 5875.621, 5865, 5869.1, 5879.9, 5882]]


def equivalent_width(wl,flux,linecenter,snr):
    # print wl, linecenter
    velo,vsini = airmass.wl_to_velocity(wl,linecenter)
    # print vsini
    vlim=175
    wl_linepart = wl[(velo>-vlim)&(velo<vlim)]
    dwl = wl_linepart[-1]-wl_linepart[0]
    v_linepart = velo[(velo>-vlim)&(velo<vlim)]
    # print velo
    flux_linepart = flux[(velo>-vlim)&(velo<vlim)]
    # print flux_linepart
    F_avg = np.average(flux_linepart)
    ew =dwl*(1-F_avg)
    # print F_avg,dwl,ew,snr
    er = np.sqrt(1+(1/F_avg))*(dwl-ew)/snr
    # print er
    return er, ew

def aphase(filetime):
    # datafile = pf.open(file)
    period = 6.829
    jdstart = Time(2452734.2, format='jd').jd
    # header =  datafile[0].header
    # filetime = Time(header['MJD-OBS'], format='jd').jd
    phase = ((filetime- jdstart)% period )/ period
    return phase

def my_sin2(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*x  + phase)  + offset




# samenightlist = []
# for key, group in itertools.groupby(filelist_apo, lambda x: airmass.date(x)):
#     # print key, list(igroup)
#     samenightlist.append(list(group))
#
# print samenightlist[0]
# print samenightlist[1]
#
# for line in linelist_apo:
#     ews = []
#     ebs = []
#     phases = []
#     for night in samenightlist:
#         # print len(night)
#         ews1n = []
#         eb1n=[]
#         phases1n = []
#         for file in night:
#             wl_rebin,flux_rebin = airmass.reduce_spectrum(file,radial_velocity=18.5)
#
#             wln, fluxn, _ = airmass.normalize(wl_rebin,flux_rebin,line[2],line[3],line[4],line[5],line[2],line[5])
#             snr = airmass.snr(file)
#             er,ew = equivalent_width(wln,fluxn,line[1],snr)
#             # print ew
#             b_cor,HJD = airmass.barcor(file)
#             phase = aphase(HJD)
#             ews1n.append(ew)
#             eb1n.append(er)
#             phases1n.append(phase)
#         ph = np.average(phases1n)
#         eqw = np.average(ews1n)
#         eb = np.sqrt(np.sum(np.array(eb1n)**2))/len(eb1n)
#         ews.append(eqw)
#         ebs.append(eb)
#         phases.append(ph)
#     print line
#     print ews
#     print ebs
#     print phases
#     dat = np.array([ews,ebs,phases,line])
#     np.save(r'C:\peter\School\Master Scriptie\Data\EW\samenight\\'+line[0]+'.npy',dat)
#-------------------------------------------------------------------------------------------------
for line in linelist_apo:
    dat = np.load(r'C:\peter\School\Master Scriptie\Data\EW\samenight\\'+line[0]+'.npy')
    ew = dat[0]
    eb = dat[1]
    phases = dat[2]
    fit = curve_fit(my_sin2, phases, ew, sigma=eb, absolute_sigma=True)
    t= np.linspace(0,1,50)
    data_fit = my_sin2(t, *fit[0])
    plt.title(line[0]+str(int(line[1])))
    plt.xlabel('Phase (P=6.29d)')
    plt.ylabel(r'Equivalent width ($\AA$)')
    plt.errorbar(phases,ew,yerr=eb,fmt='o')
    plt.plot(t,data_fit)
    plt.savefig(r'C:\peter\School\Master Scriptie\figures\EWs\eShel\samenight\\'+line[0]+str(int(line[1]))+'.pdf')


# ------LAPALMA------------------------------------------------------------------------------------------------
# all_ews = []
# all_phases = []
# ebs = []
# for file in filelist_lapalma:
#     datafile = pf.open(file)
#     header = datafile[0].header
#     snr = header['SNR50']
#     phase =  aphase(header['BJD'])
#     naxis1 = header['NAXIS1']
#     crval1 = header['CRVAL1']
#     cdelt1 = header['CDELT1']
#     flux = datafile[0].data
#     wl = np.exp(np.arange(naxis1)*cdelt1 + crval1 )
#     ews = []
#     errs = []
#     linesused = []
#     for line in [ll2[2]]:
#         print line[0],line[1]
#         wln, fluxn, _ = airmass.normalize(np.array(wl),np.array(flux),line[2],line[3],line[4],line[5],line[2],line[5])
#         er,ew = equivalent_width(wln,fluxn,line[1],snr)
#         ews.append(ew)
#         errs.append(er)
#         linesused.append(line)
#     eb = np.sqrt(np.sum(np.array(errs)**2))/len(errs)
#     ebs.append(eb)
#     equi_width = np.average(ews)
#     all_ews.append(equi_width)
#     all_phases.append(phase)
# all_ews_n = np.array(all_ews)/np.average(all_ews)
# dat = np.array([all_ews_n,all_phases,ebs,linesused])
#
# np.save('C:\peter\School\Master Scriptie\Data\EW\lapalmahe5875.npy',dat)
# plt.scatter(all_phases,all_ews_n)
# plt.show()
# # ------------APO-----------------------------------------------------------
# #
# all_ews = []
# all_phases = []
# ebs = []
# for file in filelist_apo:
#     snr = airmass.snr(file)
#     ews = []
#     errs = []
#     b_cor,HJD = airmass.barcor(file)
#     phase = aphase(HJD)
#     wl_rebin,flux_rebin = airmass.reduce_spectrum(file,radial_velocity=18.5)
#     linesused = []
#     for line in [ll3[1]]:
#         print line[0],line[1]
#         wln, fluxn, _ = airmass.normalize(wl_rebin,flux_rebin,line[2],line[3],line[4],line[5],line[2],line[5])
#         # print fluxn
#         er,ew = equivalent_width(wln,fluxn,line[1],snr)
#         ews.append(ew)
#         errs.append(er)
#         linesused.append(line)
#     eb = np.sqrt(np.sum(np.array(errs)**2))/len(errs)
#     ebs.append(eb)
#     # print ebs
#     equi_width = np.average(ews)
#     all_ews.append(equi_width)
#     all_phases.append(phase)
# all_ews_n = np.array(all_ews)/np.average(all_ews)
# dat = np.array([all_ews_n,all_phases,ebs,linesused])
#
# np.save('C:\peter\School\Master Scriptie\Data\EW\eshelhe5875-5.npy',dat)
# plt.scatter(all_phases,all_ews_n)
# plt.show()
# # ------------Fitting-----------------------------------------------------------
#
#
# def my_sin(x, freq, amplitude, phase, offset):
#     return np.sin(x * freq + phase) * amplitude + offset
# def my_sin2(x,  amplitude, phase, offset):
#     return amplitude*np.sin(2*np.pi*x  + phase)  + offset
#
# # dat = np.load('C:\peter\School\Master Scriptie\Data\EW\eshel4.npy')
# # ew = dat[0]
# # eb = dat[2]
# # phases = dat[1]
# # dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lapalmahe5875.npy')
# # dat2 = np.load('C:\peter\School\Master Scriptie\Data\EW\eshelhe5875-5.npy')
# dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lapalmahb.npy')
# dat2 = np.load('C:\peter\School\Master Scriptie\Data\EW\eshelhb-5.npy')
#
# linesused = dat1[3]
# linenames = ''
# for line in linesused:
#     linenames+=' ' + line[0]
#     linenames+= str(int(round(line[1])))
# # ew = np.concatenate((dat1[0],dat2[0]))
# # eb = np.concatenate((dat1[2],dat2[2]))
# ew1 = dat1[0]
# ew2 = dat2[0]
# eb1 = dat1[2]
# eb2 = dat2[2]
# phases1 = dat1[1]
# phases2 = dat2[1]
# # print ew
# print len(ew2),len(eb2),len(phases2)
# p0=[4, 0.005,0.5, 1]
# p1=[0.005,0.5, 1]
# fit2 = curve_fit(my_sin2, phases2, ew2,  p0=p1, sigma=eb2, absolute_sigma=True)
# # print fit[0][0]/(2*np.pi)
# fit1 = curve_fit(my_sin2, phases1, ew1, p0=p1,  sigma=eb1, absolute_sigma=True)
# print fit1[0]
# print fit2[0]
#
# t= np.linspace(0,1,50)
# data_fit = my_sin2(t, *fit1[0])
# data_fit2 = my_sin2(t, fit1[0][0],fit1[0][1]%(2*np.pi),fit1[0][2])
# print 'phaselapalma = ' , (fit1[0][1]%(2*np.pi))/(2*np.pi)
# print 'errlaplma = ', fit1[1]
# print 'phaseeshel = ' , (fit2[0][1]%(2*np.pi))/(2*np.pi)
# print 'erreshel = ', fit2[1]
# # print fit2[0][1]%(2*np.pi)
# # print (2*np.pi)+fit2[0][1]
# # plt.plot(t,data_fit,c='g')
# # plt.plot(t,data_fit2,c='r')
# # plt.show()
# # ---------------plotting----------
# t= np.linspace(0,1,50)
# data_fit = my_sin2(t, *fit1[0])
# data_fit2 = my_sin2(t, *fit2[0])
#
# f, (ax1, ax2) = plt.subplots(2, sharex=True)
# plt.suptitle(r'H$\beta$',size=18 )
# # plt.suptitle(r'He I 5875',size=18 )
# ax1.set_title('HERMES',size=12)
# ax2.set_title('eShel',size=12)
# ax1.errorbar(phases1, ew1, yerr=eb1, fmt='o', label = 'HERMES', c='r')
# ax2.errorbar(phases2, ew2, yerr=eb2, fmt='o', label = 'Data', c='b')
# ax2.plot(t,data_fit2, label='Fit', c='b')
# ax1.plot(t,data_fit, label='Fit HERMES',c='r')
# # ax2.set_xlabel('Phase (6.829 d)',size='16')
# ax1.set_ylim([0.85,1.1])
# # plt.ylabel('Equivalent width $\AA$',size='14')
# plt.xlim([0,1])
#
# f.add_subplot(111, frameon=False)
# # hide tick and tick label of the big axes
# plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
# plt.xlabel('Phase (6.829 d)',size='16')
# # plt.ylabel("common Y")
# plt.ylabel('Equivalent width $\AA$',size='14')
# # ax1.legend()
# # ax2.legend()
# plt.legend()
# # plt.savefig(r'C:\peter\School\Master Scriptie\figures\EWs\together\He5875.pdf')
# plt.savefig(r'C:\peter\School\Master Scriptie\figures\EWs\together\Hb.pdf')
# plt.show()
# # ---------------EW-EW-Comparison-------------------------------------------------
# #
# dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lapalmahb.npy')
# dat2 = np.load('C:\peter\School\Master Scriptie\Data\EW\lapalmahe5875.npy')
# # dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\eshelhb-5.npy')
# # dat2 = np.load('C:\peter\School\Master Scriptie\Data\EW\eshelhe5875-5.npy')
#
# ew1 = np.array(dat1[0])
# phase1 = np.array(dat1[1])
# phase2 = np.array(dat2[1])
# ew2 = np.array(dat2[0])
# eb1 = np.array(dat1[2])
# eb2 = np.array(dat2[2])
# # plt.errorbar(ew1, ew2, xerr=eb2,yerr=eb2, fmt='o')
# # plt.title(r'EW comparision of He I 5876 and H$\beta$')
# # plt.xlabel(r'EW H$\beta$ ($\AA$)')
# # plt.ylabel(r'EW He I 5875 ($\AA$)')
# # plt.show()
#
# pearson_cor = ss.pearsonr(ew1,ew2)
#
# print 'pearson cor = ', pearson_cor
#
# def func(x,a,b):
#     return a*x+b
# #
# # w = np.random.rand(36)
# # w1 = 1/(eb1**2)
# # w2 = 1/(eb2**2)
# # w2/=np.sum(w1)
# # w2/=np.sum(w2)
# #
# #     # Actual weighting
# # print type(ew1)
# # x = ew1*(1/eb1**2)
# # y =ew2*(1/eb2**2)
# # x = ew1*w1
# # y =ew2*w2
# # print x.shape, y.shape
# # print np.corrcoef(ew1,ew2)
# # print np.cov(ew1,ew2)
# # print np.cov(ew1,ew2,ddof=1)*sum(ew1**2)
# # print np.corrcoef(x,y)
# pars, corr = curve_fit(func, ew1, ew2, p0=[1, 0], sigma=eb1,absolute_sigma=True)
# # print 'fit corr = ', corr
# plt.errorbar(ew1, ew2, xerr=eb2,yerr=eb2, fmt='o')
# # for i in range(len(ew1)):
# #     plt.text(ew1[i], ew2[i], str(i+1), color="red", fontsize=12)
# plt.plot([0.8,1.2],func(np.array([0.8,1.2]), *pars),c = 'red')
# plt.title(r'EW comparision of He I 5876 and H$\beta $ (HERMES)')
# plt.xlim([0.85,1.1])
# plt.xlim([0.85,1.1])
# plt.xlabel(r'EW H$\beta$ ($\AA$)')
# plt.ylabel(r'EW He I 5875 ($\AA$)')
# plt.savefig(r'C:\peter\School\Master Scriptie\figures\EWs\lapalma\ewcomp.pdf')
# plt.show()