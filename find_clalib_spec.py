# __author__ = 'PeterBacily'
from __future__ import division
import warnings
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
from SavitzkyGolay import savitzky_golay

fl_eshel_clean_folder = r'D:\Peter\Master Thesis\Data\eShelData\data\clean'
filelist = glob.glob(fl_eshel_clean_folder+r'\*.fit')
filelist_lapalma = glob.glob('D:\Peter\Master Thesis\Data\LaPalmaData/*.fits')
filelist2 =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
ll_lapalma = [[r'Ha', 6562.819, 6551, 6552, 6578, 6579], [r'Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hy', 4340.472, 4322, 4324, 4357, 4360], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8], ['He_II', 4541.6, 4498, 4499, 4580, 4581], ['He_II', 4685.804, 4679, 4680, 4690, 4691], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1]]
workfile = 'C:/peter/School/Master Scriptie/Data/EW/lp/lp_snr_ll'
ll_bad_ones = [[r'Ha', 6562.819, 6551, 6552, 6578, 6579], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8], ['He_II', 4685.804, 4679, 4680, 4690, 4691], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1]]

# ll_n = [['He', 6678.15, 6656, 6660, 6690, 6695]]
ll_n =[['He', 4921.93, 4910, 4913, 4928.2, 4931.5]]

# start = 6500
# stop = 6600
plt.style.use('classic')
testfile = filelist_lapalma[0]
testfile_eshel = filelist[8]
def real_snr(file,linecenter,continuumpiece = None):
    linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, [6615.0, 6618.3]], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, [4877.0, 4883.0]], ['Hy', 4340.472, 4322, 4324, 4357, 4360, [4298.0, 4302.0]], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, [4031.0, 4038.0]], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, [4458.0, 4463.0]], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, [4734.0, 4740.0]], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, [5819.0, 5833.0]], ['He_II', 4541.6, 4498, 4499, 4580, 4581, [4579.2, 4587.6]], ['He_II', 4685.804, 4679, 4680, 4690, 4691, [4718.0, 4740.0]], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, [5370.0, 5380.0]], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, [5564.3, 5575.0]], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, [5831.2, 5833.1]]]
    sl = sorted(linelist, key=lambda x: abs(x[1]-linecenter))
    print sl
    if continuumpiece==None:
        closest_line = sl[0]
        start = closest_line[6][0]
        stop = closest_line[6][1]
    else:
        start = continuumpiece[0]
        stop = continuumpiece[1]
    dif = np.min(np.absolute(np.array([start,stop])-linecenter))
    if dif > 100:
        warnings.warn('Warning: SNR calculated from continuum piece '+str(dif)+' Angstrom from linecenter, SNR might be inaccurate')

def normalize_fluxarray(flux):
    wl = range(len(flux))
    fitparams = np.polyfit(wl,flux,1)
    x1 = fitparams[0]
    x0 = fitparams[1]
    # fitted_flux = np.array(wl)*x1 +x0
    fitted_flux = savitzky_golay(flux,25,4)
    normalized_flux = flux/fitted_flux
    return normalized_flux

def find_snr_continuum_piece_eshel(file,ll):
    wl, flux = airmass.extractdata(35, file)
    snr_ll = []
    # for line in ll_lapalma:
    for line in ll:
        dlambda = 50
        a,b,c,d = line[2],line[3],line[4],line[5]
        linewave, fluxarray, nnf = airmass.normalize(np.array(wl), np.array(flux), a,b,c,d,line[2]-dlambda,line[5]+dlambda)
        sat = 'n'
        while sat == 'n':
            calibflux = np.hstack((fluxarray[(linewave>a)&(linewave<b)],fluxarray[(linewave>c)&(linewave<d)]))
            f,(ax1,ax2,ax3) = plt.subplots(3)
            ax1.plot(linewave,fluxarray)
            ax1.axvline(a)
            ax1.axvline(b)
            ax1.axvline(c)
            ax1.axvline(d)
            ax2.plot(calibflux)
            plt.show()
            plt.close()
        # # print type(linewave),type(fluxarray)

        # # print calibflux
        # avg = np.average(calibflux)
        # std = np.std(calibflux)
        # # print avg, std
        # snr = avg/std
        # print '######################'
        # print 'line: ', line[0],line[1]
        # print 'Oude SNR:', snroud
        # print 'nieuwe snr:  SNR: ', snr
        # print '######################'
        ##########################

        ##########################
        # start = line[1]-dlambda
        # stop = line[1]+dlambda

        # flux_line = flux[(wl>start) & (wl<stop)]
        # wl_line = wl[(wl>start) & (wl<stop)]

            # plt.plot(wl_line,flux_line)
            # plt.axvline(line[1])
            # plt.show()
            # plt.close()
            v1 = float(raw_input("start: "))
            plt.plot(linewave,fluxarray)
            plt.axvline(v1)
            plt.show()
            plt.close()
            v2 = float(raw_input("stop: "))
            plt.close()
            print linewave
            snrflux = fluxarray[(linewave>v1) & (linewave<v2)]
            norm_snrflux = normalize_fluxarray(snrflux)
            f,(ax1,ax2) = plt.subplots(2)
            ax1.plot(linewave,fluxarray)
            ax1.axvline(v1)
            ax1.axvline(v2)
            ax2.plot(norm_snrflux)
            plt.show()
            plt.close()
            sat = raw_input('satisfied?')
            newline = line
            print newline
            newline.append([v1,v2])
        snr_ll.append(newline)
        print newline
    #     snr_ll.append(newline)
    print snr_ll
    # np.savetxt(r'C:\peter\School\Master Scriptie\Data\snr_linelist.txt',snr_ll)


def find_snr_continuum_piece(file,ll):
    datafile = pf.open(file)
    header = datafile[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    snroud = header['SNR50']
    rawflux = np.array(datafile[0].data)
    rawwl = np.exp(np.arange(naxis1)*cdelt1 + crval1 - 18.5/299792.458)
    wl, flux = airmass.remove_nan(rawwl,rawflux)
    snr_ll = []
    # for line in ll_lapalma:
    for line in ll:
        dlambda = 50
        a,b,c,d = line[2],line[3],line[4],line[5]
        linewave, fluxarray, nnf = airmass.normalize(np.array(wl), np.array(flux), a,b,c,d,line[2]-dlambda,line[5]+dlambda)
        sat = 'n'
        while sat == 'n':
            calibflux = np.hstack((fluxarray[(linewave>a)&(linewave<b)],fluxarray[(linewave>c)&(linewave<d)]))
            f,(ax1,ax2,ax3) = plt.subplots(3)
            ax1.plot(linewave,fluxarray)
            ax1.axvline(a)
            ax1.axvline(b)
            ax1.axvline(c)
            ax1.axvline(d)
            ax2.plot(calibflux)
            plt.show()
            plt.close()
        # # print type(linewave),type(fluxarray)

        # # print calibflux
        # avg = np.average(calibflux)
        # std = np.std(calibflux)
        # # print avg, std
        # snr = avg/std
        # print '######################'
        # print 'line: ', line[0],line[1]
        # print 'Oude SNR:', snroud
        # print 'nieuwe snr:  SNR: ', snr
        # print '######################'
        ##########################

        ##########################
        # start = line[1]-dlambda
        # stop = line[1]+dlambda

        # flux_line = flux[(wl>start) & (wl<stop)]
        # wl_line = wl[(wl>start) & (wl<stop)]

            # plt.plot(wl_line,flux_line)
            # plt.axvline(line[1])
            # plt.show()
            # plt.close()
            v1 = float(raw_input("start: "))
            plt.plot(linewave,fluxarray)
            plt.axvline(v1)
            plt.show()
            plt.close()
            v2 = float(raw_input("stop: "))
            plt.close()
            print linewave
            snrflux = fluxarray[(linewave>v1) & (linewave<v2)]
            norm_snrflux = normalize_fluxarray(snrflux)
            f,(ax1,ax2) = plt.subplots(2)
            ax1.plot(linewave,fluxarray)
            ax1.axvline(v1)
            ax1.axvline(v2)
            ax2.plot(norm_snrflux)
            plt.show()
            plt.close()
            sat = raw_input('satisfied?')
            newline = line
            print newline
            newline.append([v1,v2])
        snr_ll.append(newline)
        print newline
    #     snr_ll.append(newline)
    print snr_ll
    # np.savetxt(r'C:\peter\School\Master Scriptie\Data\snr_linelist.txt',snr_ll)


def check_TVS(file, filelist):
    linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, [6615.0, 6618.3]], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, [4877.0, 4883.0]], ['Hy', 4340.472, 4322, 4324, 4357, 4360, [4298.0, 4302.0]], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, [4031.0, 4038.0]], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, [4458.0, 4463.0]], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, [4734.0, 4740.0]], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, [5819.0, 5833.0]], ['He_II', 4541.6, 4498, 4499, 4580, 4581, [4579.2, 4587.6]], ['He_II', 4685.804, 4679, 4680, 4690, 4691, [4718.0, 4740.0]], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, [5370.0, 5380.0]], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, [5564.3, 5575.0]], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, [5831.2, 5833.1]]]
    datafile = pf.open(file)
    header = datafile[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    snroud = header['SNR50']
    rawflux = np.array(datafile[0].data)
    rawwl = np.exp(np.arange(naxis1)*cdelt1 + crval1 - 18.5/299792.458)
    wl, flux = airmass.remove_nan(rawwl,rawflux)
    snr_ll = []
    for line in ll_lapalma:
        dlambda = 50
        a,b,c,d = line[2],line[3],line[4],line[5]

        linewave, fluxarray, nnf = airmass.normalize(np.array(wl), np.array(flux), a,b,c,d,line[2]-50,line[5]+50)
        sat = 'n'
        while sat == 'n':
            calibflux = np.hstack((fluxarray[(linewave>a)&(linewave<b)],fluxarray[(linewave>c)&(linewave<d)]))
            f,(ax1,ax2,ax3) = plt.subplots(3)
            ax1.plot(linewave,fluxarray)
            ax1.axvline(a)
            ax1.axvline(b)
            ax1.axvline(c)
            ax1.axvline(d)
            ax2.plot(calibflux)
            plt.show()
            plt.close()

# find_snr_continuum_piece_eshel(testfile_eshel,ll_n)
# find_snr_continuum_piece(testfile,ll_n)

snr_l= [['-','_',5375,5385,5385,5385]]
find_snr_continuum_piece_eshel(testfile_eshel,snr_l)