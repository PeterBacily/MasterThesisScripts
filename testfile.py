from __future__ import division
import os
# os.environ['PYSYN_CDBS'] = 'C:\Users\Peter\Anaconda2\envs\p27\Lib\site-packages\pysynphot\data\cbds\grp\hst\cdbs'
os.environ.get('PYSYN_CDBS')
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
from scipy.optimize import *
from scipy.stats import chi2
from PyAstronomy import pyasl
import matplotlib.style
import datareduc
import pickle
import warnings
from SavitzkyGolay import savitzky_golay
import open_masterfiles
# fl_eshel_clean_folder = r'D:\Peter\Master Thesis\Data\eShelData\data\clean' #Zonder twee spectra met rare Halpha spike
# fl_clean = glob.glob(fl_eshel_clean_folder+r'\*.fit')
# ll_TVS_eshel = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
# # dict = {}
# # dict[1]=['a','b']
# # dict[2]=['c','d']
# #
# # print dict[2]
# # dict = {'Python' : [1,2], 'C++' : [3,4], 'Java' : [5,6]}
# # f = open("dict.txt","w")
# # filepath = r'D:\Peter\Master Thesis\Data\masterfiles\dict_apo_files.txt'
# fl_eshel_goodSNR_folder = r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\eShelData\data'
# filelist = glob.glob(fl_eshel_goodSNR_folder + r'\*.fit')
# # filelist = glob.glob(rMaster Thesis\Data\LaPalmaData\*.fits')
# # f = open(filepath,'r')
# # dict = ast.literal_eval(f.read())
# # f.close()
# # print dict
# # BCCor,HJD = airmass.barcor(fl_clean[9])
# # HJD_rounded = round(HJD-2457000,3)
# # print dict[HJD_rounded]
# # print dict['Python']
# # f.write( str(dict) )
# # f.close()
# # line = ll_TVS_eshel[1]
# # file1 = fl_clean[0]
# # file2=fl_clean[1]
# # swl = line[3] - 40
# # ewl = line[6] + 40
# # vs,lfs = airmass.overplot([file1,file2],'-',line, '-',v_rad=18.5,startwl=swl,endwl=ewl)
# # f, (ax1) = plt.subplots(1, sharex=True)
# # for i, spec in enumerate(lfs):
# #     ax1.plot(vs[i], spec, linewidth=1.0)
# # ax1.set_title(line[7])
# # plt.show()
# # ax1.legend()
# # ax1.set_xlim([-600,600])
# # spec2 = spec[(v > -300) & (v < 300)]
# # mini = np.floor(10 * 0.9 * np.amin(spec2)) / 10
# # maxi = np.ceil(10 * 1.01 * np.amax(spec2)) / 10
# # ax1.set_ylim([mini, maxi])
# # ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
# # ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
#
# # datafile = pf.open(file)
# # header = datafile[0].header
# # for item in header:
# #     print item
# # snr = header['SNR50']
# # print type(snr), snr
#     # phase =  aphase(header['BJD'])
#     # naxis1 = header['NAXIS1']
#     # crval1 = header['CRVAL1']
#     # cdelt1 = header['CDELT1']
#
# # a = 0.9
# # b= 0.1
# #
# # print ((b-a)+1)%1
#
# # a = np.array([1,2,3])
# #
# # print a.shape
#
# # file = fl_clean[3]
# # wl1,flux1 = airmass.reduce_spectrum(file)
# # wl2,flux2 = airmass.reduce_spectrum2(file)
# #
# # plt.plot(wl1,flux1,color='red',linewidth = 1.0)
# # plt.plot(wl2,flux2,color = 'blue',linewidth = 1.0)
# # plt.ylim([0,1.4])
# # plt.show()
# # fl = glob.glob(r'D:\Peter\Master Thesis\Data\masterfiles\test\*.txt')
# #
# # testfile = open(fl[0],'r')
# # b=pickle.load(testfile)
# #
# # plt.plot(b.line4026.v,b.line4026.flux)
# # plt.show()
#
# # linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
# # ll_TVS_eshel = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
# #
# # masterfilelist = open_masterfiles.apo()
# # pdffile = matplotlib.backends.backend_pdf.PdfPages(r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\figures\testfile.pdf')
# # for file in masterfilelist:
# #     f, (ax1, ax2) = plt.subplots(2, sharex=True)
# #     flux = file.line4861.flux
# #     v = file.line4861.v
# #     vrad = -18.5
# #     bcor = file.baricentric_correction
# #     v2 =  v+(vrad+bcor)
# #     sg_f = savitzky_golay(flux, 25, 4)
# #     ax1.plot(v2,sg_f)
# #     ax2.plot(v2,flux)
# #     pdffile.savefig(f)
# #     print 'da'
# #
# # pdffile.close()
#
# # line = ll_TVS_eshel[1]
# # k=2
# # for file in filelist:
# #     wl, flux = airmass.extractdata(35,file)
# #     lw, lf, nf = airmass.normalize(wl, flux, line[k + 1], line[k + 2], line[k + 3], line[k + 4], line[k + 1] - 20,
# #                                    line[k + 4] + 20)
# #     v, vsini = airmass.wl_to_velocity(lw, line[k])
# #     sg_f = savitzky_golay(lf, 25, 4)
# #     plt.plot(v,sg_f)
# # plt.show()
# # plt.close()
# # file = fl_clean[3]
# #
# # wl, flux = airmass.extractdata(35, file)
# # wlmin = wl[0]
# #
# # linelist2 = [line for line in linelist if (line[2]>wl[0] and line[5]<wl[-1])]
# #
# # # print linelist
# # #
# # # print wl[0]
# # #
# # # print wl[-1]
# # #
# # # print linelist2
# # #
# # # for line in linelist:
# # #     if not (line[2]>wl[0] and line[5]<wl[-1]):
# # #         print line
# # ll_TVS_eshel = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
# # #
# # # for line in linelist2:
# # #     dxlist = [np.abs(i[2]-line[1]) for i in ll_TVS_eshel]
# # #     mdx = min(dxlist)
# # #     if mdx >2.:
# # #         print line,',',
# # # print len(linelist2), len(ll_TVS_eshel)
# #
# # # lleft = [i for i, j in zip(linelist2, ll_TVS_eshel) if np.abs(i[2]-j[1])>2.]
# # #
# # # print lleft
# #
# #
# #
# # apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']
# #
# # mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471', 'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
# #
# # # for line in apo_lines:
# # #     if line in mercator_lines:
# # #         print line,' yes'
# # #     else:
# # #         print line, ' no'
# # apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
# #              'line4861', 'line4921', 'line6678', 'line4471']
# #
# # lines2 = ['line6562','line4861','line4471', 'line4713' ,'line4921', 'line5875', 'line6678','line4541','line4685','line5411','line5801','line5592']
# #
# # for line in lines2:
# #     if line not in apo_lines:
# #         print line
#
# # a =12
# #
# # if a>100:
# #     print 'yes'
# # elif a>20:
# #     print 'doubleyes'
# # else:
# #     print 'no'
# L_1 = 110
# L_2 =100
#
# AIC_1 = 2 * 2 - 2 * np.log(L_1)
# AIC_2 = 2 * 2 - 2 * np.log(L_2)
# probfactor = np.exp((AIC_1 - AIC_2) / 2)
# print 1/probfactor

# list = ['a','b','c','d','e']
#
# print list[:2]
# print list[2:]

apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
                      'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']


# apo_master_files = open_masterfiles.apo()
# master_files = apo_master_files
# line =apo_lines[0]
# wl, TVS, v, n = airmass.TVS_masterfiles(master_files, line)
# vsini = 127
# plt.plot(wl,TVS)
# plt.show()
# plt.close()
# plt.plot(v,TVS)
# plt.show()
# plt.close()

datareduc.plot_TVS_Lapalma_masterfile(r'D:\peter\Master_Thesis\Datareduction\Plots\test',show='on',save='off',sg='off',oneline='off')