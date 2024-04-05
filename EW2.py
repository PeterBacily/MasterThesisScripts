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
from PyAstronomy import pyasl


c_light = 299792.458
# ['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],
linelist_apo = [ ['H_beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2],  ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
linelist = [4861.01,4471.47,6678.15,4921.93,4647.42,5015.68,5411.52,4650.85,4541.52,5592.25]
ll2= [['H_beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0],['H_gamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0],['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0],['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914],['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471],['He_II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],['C_IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
filelist_lapalma = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
filelist_apo =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
nall = [['Na doublet',5892.95,5888.5,5889,5894.2,5894.75]]
def equivalent_width(wl,flux,linecenter,snr):
    # print wl, linecenter
    velo,vsini = airmass.wl_to_velocity(wl,linecenter)
    # print vsini
    wl_linepart = wl[(wl>5889.2)&(wl<5891)]
    dwl = wl_linepart[-1]-wl_linepart[0]
    v_linepart = velo[(wl>5889.5)&(wl<5891)]
    # print velo
    flux_linepart = flux[(wl>5889.5)&(wl<5891)]
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
# ------LAPALMA------------------------------------------------------------------------------------------------
# all_ews = []
# all_phases = []
# ebs = []
for file in filelist_lapalma:
    datafile = pf.open(file)
    header = datafile[0].header
    # for item in header:
    #      print item
    date = header['DATE-OBS']
    JD = header['BJD']
    print date, JD
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
#     for line in nall:
#         print line[0],line[1]
#         wln, fluxn, _ = airmass.normalize(np.array(wl),np.array(flux),line[2],line[3],line[4],line[5],line[2],line[5])
#         er,ew = equivalent_width(wln,fluxn,line[1],snr)
#         ews.append(ew)
#         errs.append(er)
#         linesused.append(line)
#     eb = np.sqrt(np.sum(np.array(errs)**2))/len(errs)
#     # ebs.append(eb)
#     equi_width = np.average(ews)
#     all_ews.append(equi_width)
#     all_phases.append(phase)
# all_ews_n = np.array(all_ews)/np.average(all_ews)
# dat = np.array([all_ews_n,all_phases,ebs,linesused])
# #
# np.save('C:\peter\School\Master Scriptie\Data\EW\lapalmana5892.npy',dat)
# plt.scatter(all_phases,all_ews_n)
# plt.show()
# ------------APO-----------------------------------------------------------
#
# all_ews = []
# all_phases = []
# ebs = []
# for file in filelist_apo[1:]:
#     snr = airmass.snr(file)
#     ews = []
#     errs = []
#     phase = airmass.phase(file)
#     wl_rebin,flux_rebin = airmass.reduce_spectrum(file,radial_velocity=18.5)
#     linesused = []
#     for line in nall:
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
# np.save('C:\peter\School\Master Scriptie\Data\EW\eshelna.npy',dat)
# plt.scatter(all_phases,all_ews_n)
# plt.show()
# ------------Fitting-----------------------------------------------------------


# def my_sin(x, freq, amplitude, phase, offset):
#     return np.sin(x * freq + phase) * amplitude + offset
# def my_sin2(x,  amplitude, phase, offset):
#     return np.sin(x * 2*np.pi + phase) * amplitude + offset
#
# # dat = np.load('C:\peter\School\Master Scriptie\Data\EW\eshel4.npy')
# # ew = dat[0]
# # eb = dat[2]
# # phases = dat[1]
# dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lapalmana.npy')
# dat2 = np.load('C:\peter\School\Master Scriptie\Data\EW\eshelna.npy')
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
#
# p0=[4, 0.005,0.5, 1]
# p1=[0.005,0.5, 1]
# # fit = curve_fit(my_sin, phases, ew,  p0=p0, sigma=eb, absolute_sigma=True)
# # print fit[0][0]/(2*np.pi)
# fit2 = curve_fit(my_sin2, phases1, ew1, p0=p1,  sigma=eb1, absolute_sigma=True)
#
# # # ---------------plotting----------
# t= np.linspace(0,1,50)
# # data_fit = my_sin(t, *fit[0])
# data_fit2 = my_sin2(t, *fit2[0])
#
# plt.title(linenames )
# plt.errorbar(phases1, ew1, yerr=eb1, fmt='o', label = 'HERMES', c='r')
# plt.errorbar(phases2, ew2, yerr=eb2, fmt='o', label = 'eShel', c='b')
# # plt.plot(t,data_fit, label='best fit')
# plt.plot(t,data_fit2, label='Expected')
# plt.xlabel('Phase (6.83 d)')
# plt.ylim([0.8,1.2])
# plt.ylabel('Equivalent width $\AA$')
# plt.xlim([0,1])
# plt.legend()
# plt.show()
