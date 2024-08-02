from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.style
# import matplotlib as mpl
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import os
import pickle
import numpy as np
import airmass
from scipy.optimize import *
from scipy.stats import chi2
import Datafile_class
import zipfile
from pathlib import Path
import shutil
from PyAstronomy import pyasl
import open_masterfiles
matplotlib.style.use('classic')
vsini =127
c_light = 299792.458
fl_lowsnr = glob.glob('C:\peter\School\Master Scriptie\Data\eShelData\data\lowsnr/*.fit')
# ['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],
linelist_apo = [ ['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2],  ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
linelist = [4861.01,4471.47,6678.15,4921.93,4647.42,5015.68,5411.52,4650.85,4541.52,5592.25]
ll2= [['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0],['H_gamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0],['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6],['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0],['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914],['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471],['He_II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],['C_IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
filelist_lapalma = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
filelist_apo =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/lowsnr*.fit')
nall = [['Na doublet',5892.95,5888.5,5889,5894.2,5894.75]]
ll3 = [[r'H$\beta$', 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I', 5875.621, 5865, 5869.1, 5879.9, 5882]]
# ll_laplma = [[r'Halpha', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'Hbeta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hgamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I', 5875.621, 5860.6241285010492, 5861.5995122694139, 5860.6241285010492, 5870.5989810949914], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4681.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5806.33334373, 5808.314546]]
ll_lapalma = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
ll = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll_lapalma2 = [[r'H$\alpha$', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'H$\beta$', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'H$\gamma$', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4546.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4690.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
# print len(ll_lapalma2)
# print len(ll_lapalma)


for line in ll_lapalma:
    if line[0] == 'He_I':
        line.append('He I '+ str(int(np.round(line[1]))))
    if line[0] == 'He_II':
        line.append('He II '+ str(int(np.round(line[1]))))


def line(x, x0, a1, b1, tau1):
    return a1 * np.exp(-tau1 * np.exp(-((x - x0) / (2 * b1)) ** 2))

def Lapalma_vrad(datafile_folder, plot_save_folder, linelist):
    filelist = glob.glob(datafile_folder+'\*.fits')
    for line in linelist:
        # print line[6]
        swl = line[2]-10
        ewl = line[5]+10
        vs,lfs =airmass.overplot_LaPalma(filelist,line,swl,ewl, v_rad=18.5)
        # print vs[0], lfs[0]
        minvs = []
        for i , v in enumerate(vs):
            # print i
            lf = lfs[i]
            index = np.argmin(lf)
            minv = v[index]
            minvs.append(minv)
        # print minvs
        print(np.median(minvs))
# Lapalma_vrad('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma',ll_lapalma)
# x = np.arange(-300,300,1)

# y= line(x,0,1,25,0.3)

# plt.scatter(x,y)
# plt.show()


# file = filelist_lapalma[0]
# datafile = pf.open(file)
# header = datafile[0].header
# for item in header:
#     print item
#
# print header['DATE-AVG']
# print header['BJD']

# file = filelist_apo[5]
# # datafile = pf.open(file)
# wl,flux = airmass.extractdata(35,file)
#
# flux1 = flux[wl<4600]
# wl1 = wl[wl<4600]
# plt.xlabel('Wavelength (A)')
# plt.ylabel('Normalized flux')
# plt.plot(wl1,flux1)
# # plt.show()
# plt.savefig(r'C:\peter\School\Master Scriptie\figures\4471_omgeving.pdf')
i=1

def aphase(filetime):
    # datafile = pf.open(file)
    period = 6.829
    jdstart = 2454391.
    erphase = int((filetime-jdstart)/period)*0.001
    unc1 = erphase/period
    unc2 = erphase+0.5/period
    # header =  datafile[0].header
    # filetime = Time(header['MJD-OBS'], format='jd').jd
    phase = ((filetime- jdstart)% period )/ period
    return phase
# print fl_lowsnr
# datfile = pf.open(fl_lowsnr[0])
# header = datfile[0].header
# # print header['DATE-OBS']
# # for item in header:
# #     print item
# print 'nr   phase   date    HJD'
# for file in filelist_apo:
#     b_cor, HJD = airmass.barcor(file)
#     datfile = pf.open(file)
#     header = datfile[0].header
#     dd= header['DATE-OBS']
#     phase = aphase(HJD)
#     wl_rebin, flux_rebin = airmass.reduce_spectrum(file, radial_velocity=18.5)
#     linesused = []
#     # print line[2]
#     snr = airmass.snr(file, 4838.0, 4839.0)
#     print i,'  ',phase,'  ',dd,'  ',HJD
#     i+=1
# print 1/3




def plot_TVS_Lapalma(datafile_folder, plot_save_folder, linelist):
    filelist = glob.glob(datafile_folder+'\*.fits')
    for line in linelist:
        print(line[6])
        swl = line[2]-40
        ewl = line[5]+40
        lw,TVS,v,n =airmass.TVS_LaPalma(filelist,line,swl,ewl, v_rad=18.5)
        p = chi2.ppf(0.99, n-1)/(n-1)
        vs,lws = airmass.overplot_LaPalma(filelist,line,v_rad=18.5,startwl=swl,endwl=ewl)
        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec )
        ax1.set_title(line[6])
        ax1.legend()
        # ax1.set_xlim([-600,600])
        # ax1.set_ylim([0.6,1.1])
        ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
        if line[2]==5875.621:
            TVS2 = np.array(TVS)*1.4
            ax2.plot(v, TVS2)
        else:
            ax2.plot(v,TVS)
        ax2.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
        ax2.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
        ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # ax2.set_ylim([0,5])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel('TVS')
        ax2.set_xlim([-600,600])
        plt.savefig(plot_save_folder + r'\\' + line[0] + str(int(np.round(line[1])))+'_TVS.pdf')
        plt.show()
        plt.close()
# testlist = open_masterfiles.apo_demetra_orders()
# import pathlib
# # testlist = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\\*.txt')
# masterfile_ll = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
#              'line4861', 'line4921', 'line6678', 'line4471']
# test_data_file = testlist[0]
# for baseline in [masterfile_ll[0]]:
#     line = baseline+'_order'
#     linedata = getattr(test_data_file,line)
#     wl =linedata.wl
#     v=linedata.v
#     v_cor=linedata.v_cor
#     flux=linedata.flux
#
#     wl,tvs,v,n=airmass.TVS_masterfiles_order(testlist,line)
#     print(len(wl))
#     print(len(v))
#     # print(len(v_cor))
#     # print(len(flux))
#     print(len(tvs))

# class Datafile_mercator_with_orders:
#     observatory = 'MERC'
#     # linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
#     linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
#                      ['Hb', 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'],
#                      ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'],
#                      ['He_I', 4026.1914, 4016, 4020, 4032, 4036, 'He I 4026'],
#                      ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'],
#                      ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'],
#                      ['He_I', 5875.621, 5863.0, 5864.5, 5884.6, 5885.5, 'He I 5875'],
#                      ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4541'],
#                      ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4685'],
#                      ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5411'],
#                      ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4921'],
#                      ['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],
#                      ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'],
#                      ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
#
#     def __init__(self, file,v_rad = 18.5,i='n/a'):
#         fn = os.path.basename(file)
#         data = pf.open(file)
#         self.i =i
#         self.filename = fn[:fn.rfind(".")]
#         self.header = data[0].header
#         self.time_and_date = airmass.timeanddate2(self.header['DATE-OBS'])
#         self.original_filepath = file
#         self.HJD = float(self.header['BJD'])
#         self.phase =  airmass.aphase(self.header['BJD'])
#         self.exptime = self.header['EXPTIME']
#         self.altitude = float(self.header['TELALT'])
#         self.airmass = 1/np.sin(2*np.pi*self.altitude/360)
#         self.baricentric_correction = float(self.header['BVCOR'])
#         self.fwl = airmass.fitfraunlp(file)
#         self.velshift = 299792.458*(self.fwl-5895.92)/5895.92
#         naxis1 = self.header['NAXIS1']
#         crval1 = self.header['CRVAL1']
#         cdelt1 = self.header['CDELT1']
#         self.wl_original = np.exp(np.arange(naxis1) * cdelt1 + crval1 - v_rad / 299792.458)
#         self.flux_original = data[0].data
#         self.wl_rebin, self.flux_rebin = airmass.rebin2(self.wl_original,self.flux_original)
#         self.available_lines = []
#         self.snr_original =airmass.snr(self.wl_original,self.flux_original)
#         self.snr = airmass.snr(self.wl_rebin,self.flux_rebin)
#         for line in self.linelist:
#             linedata,linekey = line_data(line,self.wl_rebin,self.flux_rebin,self.observatory,self.snr,0,0)
#             linedata_original, lk = line_data(line,self.wl_original,self.flux_original,self.observatory,self.snr_original,0,0)
#             setattr(self,linekey , linedata)
#             setattr(self,linekey+'_original',linedata_original)
#             self.available_lines.append(linekey)
#         data.close()

testfile = r"D:\peter\Master_Thesis\Datareduction\Data\LaPalmaData\zet Ori1120151010.fits"
fn = os.path.basename(testfile)
data = pf.open(testfile)
header = data[0].header
print(header)
data.close()
# plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma',ll_lapalma)

# folder = r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\final_spectra\spectra'
# #
# filelist = glob.glob(folder+'\*.fit')
# print(filelist)

# # fp = r'D:\peter\Master_Thesis\Master_Thesis\Other\demetra_test_2\20160304\archive\20160304-225457-Zeta_Ori-600s-1.fit'
# a=pf.open(filelist[0])
# a.info()
# b= a[0].data
# header=a[0].header
# print(header)
# print(airmass.timeanddate2(header['DATE-OBS']))
# def openlinelistfile(listpath):
#     myfile=open(listpath, 'r')
#     b=myfile.read()
#     myfile.close()
# testlist = open_masterfiles.open_linelist(r'D:\peter\Master_Thesis\Master_Thesis\Other\testlist3.txt')
# print(testlist)

# openlinelistfile(r'D:\peter\Master_Thesis\Master_Thesis\Other\testlist.txt')
# v_rad = 18.5
# naxis1 = header['NAXIS1']
# crval1 = header['CRVAL1']
# cdelt1 = header['CDELT1']
# print(naxis1,crval1,cdelt1)
# # print(dir(a[0]))
# from astropy.wcs import WCS
# w = WCS(header, naxis=1, relax=False, fix=False)
# lam = w.wcs_pix2world(np.arange(len(b)), 0)[0]
#
# wl_original = np.arange(naxis1) * cdelt1 + crval1
#               # - (v_rad / 299792.458)
# print(lam)
# print('asd',(lam==wl_original).all())
# # #
# plt.plot(lam,b)
# plt.show()
# plt.close()

# a =pf.open(fp)
# a.close()



# path_fl = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\zo_good'
# fl = open_masterfiles.apo_demetra()
# testfile = fl[2]
# line ='line4713'

# linedata = getattr(testfile, line)
# lineinfo = linedata.lineinfo
# flux = linedata.flux
# wl = linedata.wl
# v = linedata.v_cor
# lw,TVS,v,n =airmass.TVS_masterfiles(fl,line)
# print(lineinfo)
# plt.plot(v,flux)
# plt.show()
# plt.close()

# class Line:
#     def __init__(self, li,lc,wave,fl,velo,nf,vsini,snr,barcor,vrad,norm_boundaries):
#         self.lineinfo = li
#         self.normalization_boundaries_wl=norm_boundaries[0]
#         self.normalization_boundaries_v = norm_boundaries[1]
#         self.wl = wave
#         self.v = velo
#         self.v_cor = np.array(velo)+(barcor+vrad)
#         self.flux = fl
#         self.normalizationflux = nf
#         self.vsini = vsini
#         self.ew_error, self.ew = airmass.equivalent_width(self.v_cor,wave,fl,lc,snr)
#
#
# def line_data(line,wl,flux,observatory,snr,bccor,vrad):
#     if observatory == 'MERC':
#         k = 1
#         barcor = 0
#     elif observatory == 'APO':
#         k = 2
#         barcor = bccor
#     elif observatory == 'APO_DEMETRA':
#         barcor = bccor
#         k = 1
#     else:
#         print('observatory needs to be APO or MERC')
#     center_wl = int(line[k])
#     lw, lf, nf, _,_ = airmass.normalize(wl,flux,line[k+1],line[k+2],line[k+3],line[k+4],line[k+1]-20,line[k+4]+20)
#     v, vsini = airmass.wl_to_velocity(lw, line[k])
#     normalization_wl= [line[k+1],line[k+2],line[k+3],line[k+4]]
#     normalization_v = airmass.wl_to_velocity(normalization_wl, line[k])
#     return Line(line,line[k],lw,lf,v,nf,vsini,snr, barcor,vrad,[normalization_wl,normalization_v]),'line'+str(center_wl)
#
# class single_order:
#     def __init__(self, filepath,order_number,order_number_demetra,zip=True,fullzipfile=None):
#         self.order_number_demetra =order_number_demetra
#         self.order_number = order_number
#         if zip==True:
#             filefolder = zipfile.ZipFile(fullzipfile, "r")
#             a=filefolder.read(filepath)
#
#         else:
#             a = pf.open(filepath)
#         self.header = a[0].header
#         naxis1 = self.header['NAXIS1']
#         crval1 = self.header['CRVAL1']
#         cdelt1 = self.header['CDELT1']
#         self.wl_original = np.arange(naxis1) * cdelt1 + crval1
#         self.flux_original = a[0].data
#         self.wl_start = self.wl_original[0]
#         self.wl_end = self.wl_original[-1]
#         self.wl_avg = np.average([self.wl_start, self.wl_end])
#         a.close()
# def zip_to_list(filepath):
#     zipfilepath = Path(filepath)
#     tempfolder = zipfilepath.parents[0].joinpath('temp').joinpath(zipfilepath.stem)
#     print(tempfolder)
#     preexists = tempfolder.exists()
#     Path(tempfolder).mkdir(parents=True, exist_ok=True)
#     with zipfile.ZipFile(filepath, 'r') as zip_ref:
#         zip_ref.extractall(tempfolder)
#     filelist=glob.glob(str(tempfolder)+r'\*.fit')
#     return filelist,preexists,tempfolder
#
# def remove_temp_folder(tempfolder,preexists=False):
#     if preexists==False:
#         shutil.rmtree(tempfolder)
#
# def open_linelist(path):
#     a = open(path, 'rb')
#     b = pickle.load(a)
#     a.close()
#     return b
#
# class Datafile_apo_demetra_with_orders:
#     observatory = 'APO_DEMETRA'
#     # linelist = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
#     linelist_standard = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
#      ['Hb', 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'],
#      ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'],
#      ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5875'],
#      ['He_II', 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4541'],
#      ['He_II', 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4685'],
#      ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5411'],
#      ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'],
#      ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4921'],
#      ['He_I', 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'],
#      ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'],
#      ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
#
#     def __init__(self, orderfiles,fullspecfile,ll_file=None,v_rad = 18.5,i='n/a',mark = 0,zip=True):
#         if ll_file == None:
#             self.linelist = self.linelist_standard
#         else:
#             self.linelist = open_linelist(ll_file)
#
#         fn = os.path.basename(fullspecfile)
#         data = pf.open(fullspecfile)
#         self.original_filepath = fullspecfile
#         self.i =i
#         self.k=1
#         self.mark = mark
#         # self.mark_explanation = '0 = no weird stuff,      1 = not usable due to very poor SNR,    2 = Not usable for EW, TVS, Quotient due to insufficient SNR,   3 =   Shows weird feature in Halpha'
#         self.filename = fn[:fn.rfind(".")]
#         self.header = data[0].header
#         self.time_and_date = airmass.timeanddate2(self.header['DATE-OBS'])
#         self.baricentric_correction, self.HJD = airmass.barcor(fullspecfile,JDOBS=self.header['JD-MID'])
#         self.phase =  airmass.aphase(self.HJD)
#         self.exptime = airmass.exposuretime(fullspecfile)
#         self.airmass, self.alt, JD = airmass.airmass(fullspecfile,JDOBS=self.header['JD-MID'])
#         try:
#             frwl = airmass.fitfraun_demetra(fullspecfile)
#         except RuntimeError:
#             frwl = 5895.92
#         self.fwl = frwl
#         self.velshift = 299792.458*(self.fwl-5895.92)/5895.92
#         naxis1 = self.header['NAXIS1']
#         crval1 = self.header['CRVAL1']
#         cdelt1 = self.header['CDELT1']
#         self.wl_original = np.arange(naxis1) * cdelt1 + crval1 - v_rad / 299792.458
#         self.flux_original = data[0].data
#         self.wl_rebin, self.flux_rebin = airmass.rebin2(self.wl_original,self.flux_original)
#         self.available_lines = []
#         self.snr_original =airmass.snr(self.wl_original,self.flux_original)
#         self.snr = airmass.snr(self.wl_rebin,self.flux_rebin)
#         onr = 1
#         ords = []
#         if zip==True:
#             filefolder = zipfile.ZipFile(orderfiles, "r")
#             fzf = orderfiles
#             orderfilelist = filefolder.namelist()
#         else:
#             orderfilelist=orderfiles
#             fzf = None
#
#         for filename in orderfilelist:
#             file_name = os.path.basename(filename)
#             order_number_demetra= os.path.splitext(file_name)[0][-2:]
#             od = single_order(filename, order_number=onr, order_number_demetra=order_number_demetra,zip=zip,fullzipfile=fzf)
#             ords.append(od)
#             onr += 1
#         self.orders = ords
#         for line in self.linelist:
#             linedata,linekey = line_data(line,self.wl_rebin,self.flux_rebin,self.observatory,self.snr, self.baricentric_correction,-18.5)
#             linedata_original, lk = line_data(line,self.wl_original,self.flux_original,self.observatory,self.snr_original,0,0)
#             ol = sorted(self.orders, key=lambda x: np.abs(x.wl_avg-line[1]))
#             line_order = ol[0]
#             normalization_wl = [line[self.k + 1], line[self.k + 2], line[self.k + 3], line[self.k + 4]]
#             if line_order.wl_start<normalization_wl[0] and line_order.wl_end>normalization_wl[-1]:
#                 linedata_order,lk = line_data(line,line_order.wl_original,line_order.flux_original,self.observatory,self.snr, self.baricentric_correction,-18.5)
#                 setattr(self,linekey+'_order',linedata_order)
#             else:
#                 print('line out of order bounds, no order line was made for',line[0],self.filename)
#             setattr(self,linekey , linedata)
#             setattr(self,linekey+'_original',linedata_original)
#             self.available_lines.append(linekey)
#         data.close()

# test_zip_archive = r"D:\peter\Master_Thesis\Master_Thesis\Data\demetra\demetra_test\ziptest\ZetOri20160317-2_20160317T194557.zip"
# test_full_spec = r"D:\peter\Master_Thesis\Master_Thesis\Data\demetra\demetra_test\ziptest\ZetOri_03_17_2_20160317T194557.fit"
# llfile = "D:\peter\Master_Thesis\Datareduction\Converted_Data\linelists\linelist_apo.txt"
#
# testfolder = r'D:\peter\Master_Thesis\Master_Thesis\Data\demetra\demetra_test\ziptest\full_series_test\\'
# # zipfiles = glob.glob(testfolder+r'*.zip')
# # print(zipfiles)
# # full_spec_files = glob.glob(testfolder+r'*.fit')
#
# def time_from_filename(filename):
#     a = Path(filename)
#     b=str(a.stem)[-15:]
#     c=b[:8]+b[9:]
#     return(int(c))
# # time_from_filename(zipfiles[0])
# zipfiles = sorted(glob.glob(testfolder+r'*.zip'), key=lambda x: time_from_filename(x))
# full_spec_files = sorted(glob.glob(testfolder+r'*.fit'), key=lambda x: time_from_filename(x))
# savefolder = r'D:\peter\Master_Thesis\Master_Thesis\Data\demetra\demetra_test\ziptest\full_class_obj\\'
# for i in range(len(zipfiles)):
#     if str(Path(zipfiles[i]).stem)[-15:] ==str(Path(full_spec_files[i]).stem)[-15:]:
#         orders_class_object=Datafile_class.Datafile_apo_demetra_with_orders(zipfiles[i],full_spec_files[i],ll_file=llfile)
#         dl, dl2 = airmass.split_date(orders_class_object.header['DATE-OBS'])
#         savename = savefolder + orders_class_object.observatory + '_' + dl[0] + dl[1] + dl[2] + dl[3] + '.txt'
#         workfileresource = open(savename, 'wb')
#         pickle.dump(orders_class_object, workfileresource)
#         workfileresource.close()
#     else:
#         print('YOU SUCK')

# test_demetra_order_object = Datafile_class.Datafile_apo_demetra_with_orders(test_zip_archive,test_full_spec,ll_file=llfile)
# apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
#              'line4861', 'line4921', 'line6678', 'line4471']
#
# order_line_ha =test_demetra_order_object.line6562_order
# full_line_ha = test_demetra_order_object.line6562
# wlo=order_line_ha.wl
# fluxo=order_line_ha.flux
# wlf=full_line_ha.wl
# fluxf=full_line_ha.flux
# plt.plot(wlo,fluxo)
# plt.plot(wlf,fluxf)
# plt.show()
# plt.close

# print(zipfile.is_zipfile(bla))
# filelist,preexist,tempfoldername = zip_to_list(bla)
# print(filelist)
# print(preexist)
# print(tempfoldername)
# bla2=bla.parents[0]
# bla3 =bla.stem
# bla4 = bla.parents[0].joinpath('temp').joinpath(bla.stem)
# print(bla4)

# print(bla.stem)
# parents[]
# class Datafile_apo_demetra_orders:
#     observatory = 'APO_DEMETRA'
#
#     def __init__(self, filelist, ll_file):
#         self.linelist = open_linelist(ll_file)
#         for i,order_file in enumerate(filelist):
#             a = pf.open(order_file)
#             header = a[0].header
#             naxis1 = header['NAXIS1']
#             crval1 = header['CRVAL1']
#             cdelt1 = header['CDELT1']
#             wl_original = np.arange(naxis1) * cdelt1 + crval1
#             flux_original = a[0].data
#             wl_start = wl_original[0]
#             wl_end = wl_original[-1]
#             wl_avg = np.avg([wl_start,wl_end])

# linelist=open_linelist(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\linelists\linelist_apo.txt')
# file = zipfile.ZipFile("zipfile.zip", "r")
# for name in file.namelist():
#     data = file.read(name)
#
# with zipfile.ZipFile("sample.zip", mode="r") as archive:
#     archive.extractall(path="output_dir/")
#
#
# # new_list = sorted(orig_list, key=lambda x: x.count, reverse=True)
#
#
# filefolder = r'D:\peter\Master_Thesis\Master_Thesis\Data\demetra\demetra_test\single_order_test\\'
# filefolder_main= r'D:\peter\Master_Thesis\Master_Thesis\Data\demetra\demetra_test\single_order_test\full\\'
# filelist=glob.glob(filefolder+r'*.fit')
# fullspec = glob.glob(filefolder_main+r'*.fit')[0]
# filefolder_zip = r"D:\peter\Master_Thesis\Datareduction\Data\Demetra\final_spectra\back_to_back_stacked\ZetOri20160325-2_20160325T210354.zip"
# a =Datafile_apo_demetra_with_orders(filelist,fullspec)
#
#
# filename = Path('/some/path/somefile.txt')
# from pathlib import Path
# filename_wo_ext = filename.with_suffix('')
#
# lwl = a.line6562_order.wl
# lfl =a.line6562_order.flux
# plt.plot(lwl,lfl)
# plt.show()
# plt.close()

# filelist.sort(key=lambda x: os.path.splitext(os.path.basename(x))[-2:])
# i=1
# orders=[]
# for file in filelist:
#     file_name = os.path.basename(file)
#     fn2 = os.path.splitext(file_name)[0]
#     order_number_demetra = fn2[-2:]
#     print(fn2, fn2[-2:])
#     od = single_order(file,order_number=i,order_number_demetra=order_number_demetra)
#     orders.append(od)
#     i+=1
#
#
#
# testwl=4861.333
# orders.sort(key=lambda x: np.abs(x.wl_avg-testwl))
# orders.sort(key=lambda x: x.wl_start)
# for order in orders:
#     print(order.order_number,order.order_number_demetra,'@@@@', order.wl_avg,'@@@@',order.header)
# file1path=r'D:\peter\Master_Thesis\Master_Thesis\Data\demetra\demetra_test\single_order_test\ZetOri20160325-2_20160325210354_40.fit'
# a = single_order(file1path)
# header = a.header
# naxis1 = header['NAXIS1']
# crval1 = header['CRVAL1']
# cdelt1 = header['CDELT1']
# wl_original = a.wl_original
# flux_original = a.flux_original
#
# print(a.wl_start,a.wl_avg,a.wl_end)
#
# plt.plot(wl_original,flux_original)
# plt.show()
# plt.close()
