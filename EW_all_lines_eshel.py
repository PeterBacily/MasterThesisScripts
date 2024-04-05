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
import warnings


c_light = 299792.458
# ['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],
linelist_apo = [ ['Hbeta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0],['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2],  ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
linelist = [4861.01,4471.47,6678.15,4921.93,4647.42,5015.68,5411.52,4650.85,4541.52,5592.25]
ll2= [['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0],['H_gamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0],['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6],['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0],['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914],['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471],['He_II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],['C_IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
filelist_lapalma = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
filelist_apo =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
nall = [['Na doublet',5892.95,5888.5,5889,5894.2,5894.75]]
ll3 = [[r'H$\beta$', 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I', 5875.621, 5865, 5869.1, 5879.9, 5882]]
# ll_laplma = [[r'Halpha', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'Hbeta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hgamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I', 5875.621, 5860.6241285010492, 5861.5995122694139, 5860.6241285010492, 5870.5989810949914], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4681.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5806.33334373, 5808.314546]]
ll_lapalma = [[r'Ha', 6562.819, 6551, 6552, 6578, 6579], [r'Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hy', 4340.472, 4322, 4324, 4357, 4360], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8], ['He_II', 4541.6, 4498, 4499, 4580, 4581], ['He_II', 4685.804, 4679, 4680, 4690, 4691], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1]]
ll_extra = [['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0]]
print filelist_apo
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
    # print F_avg
    ew =dwl*(1-F_avg)
    # print F_avg,dwl,ew,snr
    er = np.sqrt(1+(1/F_avg))*(dwl-ew)/snr
    # print er
    return er, ew

def aphase(filetime):
    # datafile = pf.open(file)
    period = 6.829621
    jdstart = 2454391.
    erphase = int((filetime-jdstart)/period)*0.001
    unc1 = erphase/period
    unc2 = erphase+0.5/period
    # header =  datafile[0].header
    # filetime = Time(header['MJD-OBS'], format='jd').jd
    phase = ((filetime- jdstart)% period )/ period
    return phase
# ------LAPALMA------------------------------------------------------------------------------------------------
#
# for line in ll_lapalma:
#     all_ews = []
#     all_phases = []
#     ebs = []
#     ews = []
#     errs = []
#     for file in filelist_lapalma:
#         datafile = pf.open(file)
#         header = datafile[0].header
#         snr,snr_poisson = airmass.real_snr(file,line[1],SG=False)
#         # snr,snr_poisson = airmass.real_snr(file,line[1],SG=True)
#         phase =  aphase(header['BJD'])
#         naxis1 = header['NAXIS1']
#         crval1 = header['CRVAL1']
#         cdelt1 = header['CDELT1']
#         flux = datafile[0].data
#         wl = np.exp(np.arange(naxis1)*cdelt1 + crval1 )
#
#         linesused = [line]
#
#         wln, fluxn, _ = airmass.normalize(np.array(wl),np.array(flux),line[2],line[3],line[4],line[5],line[2],line[5])
#         print line
#         # print len(wln), len(fluxn)
#         # er,ew = equivalent_width(wln,fluxn,line[1],snr)
#         er,ew = equivalent_width(wln,fluxn,line[1],snr)
#         # print er,ew
#         ews.append(ew)
#         errs.append(er)
#         # linesused.append(line)
#         # print len(errs)
#         # eb = np.sqrt(np.sum(np.array(errs)**2))/len(errs)
#         # ebs.append(eb)
#         # equi_width = np.average(ews)
#         # all_ews.append(equi_width)
#         # all_ews.append(equi_width)
#         all_phases.append(phase)
#         # print all_phases, ews,errs
#     # all_ews_n = np.array(all_ews)/np.average(all_ews)
#     all_ews_n = np.array(ews)/np.average(ews)
#     # print 'marker'
#     # print all_ews_n
#     # dat = np.array([all_ews_n,all_phases,ebs,linesused])
#     dat = np.array([all_ews_n,all_phases,errs,linesused])
#     print len(all_ews_n),len(all_phases),len(errs)
#     # print dat
#     chisq, red_chisq, AIC_1,AIC_2,AIC_3,llhr = airmass.EW_stats(dat)
#     dat2 = np.array([all_ews_n,all_phases,errs,linesused,chisq, red_chisq, AIC_1,AIC_2,AIC_3,llhr])
#     linename = line[0]+str(int(line[1]))+'.npy'
#     np.save('C:\peter\School\Master Scriptie\Data\EW\lp\\'+linename,dat2)

# # ------------APO-----------------------------------------------------------
#
# print filelist_apo
for line in linelist_apo:
    all_ews = []
    all_phases = []
    ebs = []
    ews = []
    errs = []
    for file in filelist_apo:
        b_cor, HJD = airmass.barcor(file)
        phase = aphase(HJD)
        wl_rebin, flux_rebin = airmass.reduce_spectrum(file, radial_velocity=18.5)
        linesused = []
        print line[2]
        snr = airmass.snr(file,line[2],line[3])
        # print line[0],line[1]
        wln, fluxn, _ = airmass.normalize(wl_rebin,flux_rebin,line[2],line[3],line[4],line[5],line[2],line[5])
        # print fluxn
        er,ew = equivalent_width(wln,fluxn,line[1],snr)
        ews.append(ew)
        errs.append(er)
        all_phases.append(phase)
        linesused.append(line)
    # eb = np.sqrt(np.sum(np.array(errs)**2))/len(errs)
    # ebs.append(eb)
    # print ebs
    all_ews_n = np.array(ews)
    # all_ews_n = np.array(ews) / np.average(ews)
    # equi_width = np.average(ews)
    # all_ews.append(equi_width)

# all_ews_n = np.array(all_ews)/np.average(all_ews)
#     dat = np.array([all_ews_n,all_phases,ebs,linesused])
    dat = np.array([all_ews_n,all_phases,errs,linesused])
    linename = line[0] + str(int(line[1]))
    chisq, red_chisq, AIC_1, AIC_2, AIC_3, llhr = airmass.EW_stats(dat)
    dat2 = np.array([all_ews_n, all_phases, errs, linesused, chisq, red_chisq, AIC_1, AIC_2, AIC_3, llhr])
    np.save('C:\peter\School\Master Scriptie\Data\EW\\apo\\'+linename+'_nn.npy',dat2)
# plt.scatter(all_phases,all_ews_n)
# plt.show()
# # ------------Fitting-----------------------------------------------------------
#
#
# def my_sin(x, freq, amplitude, phase, offset):
#     return np.sin(x * freq + phase) * amplitude + offset
def my_sin2(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*x  + phase)  + offset
#
# # dat = np.load('C:\peter\School\Master Scriptie\Data\EW\eshel4.npy')
# # ew = dat[0]
# # eb = dat[2]
# # phases = dat[1]
# for line in [ll_lapalma[0]]:
#     linename = line[0]+str(int(line[1]))
#     linename2 = linename[:-4]+' '+linename[-4:]
#     dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lp\\'+linename+'.npy')
#
#     ew1 = dat1[0]
#
#     eb1 = dat1[2]
#
#     phases1 = dat1[1]
#     p1=[0.005,0.5, 1]
#
#     fit1 = curve_fit(my_sin2, phases1, ew1, p0=p1,  sigma=eb1, absolute_sigma=True)
#     t= np.linspace(0,1,50)
#
#
#     t= np.linspace(0,1,50)
#     data_fit = my_sin2(t, *fit1[0])
#     # data_fit2 = my_sin2(t, *fit2[0])
#     #
#     f, (ax1) = plt.subplots(1, sharex=True)
#
#     plt.suptitle(linename2,size=18 )
#     ax1.errorbar(phases1, ew1, yerr=eb1, fmt='o', label = 'HERMES', c='r')
#
#     ax1.plot(t,data_fit, label='Fit HERMES',c='r')
#
#     # ax1.set_ylim([0.85,1.1])
#     plt.ylabel('Equivalent width $\AA$',size='14')
#     plt.xlim([0,1])
#     plt.xlabel('Phase (6.829 d)',size='16')
#
#     plt.savefig(r'C:\peter\School\Master Scriptie\figures\EWs\lapalma\\'+linename+'.pdf')
# plt.show()
# # ---------------EW-EW-Comparison-------------------------------------------------
#
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

