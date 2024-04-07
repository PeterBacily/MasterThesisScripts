from __future__ import division
import glob
import astropy.io.fits as pf
import airmass
import matplotlib.style
import pickle
import Datafile_class
import os
import Path_check
import pathlib
# import matplotlib.pyplot as plt
# from astropy.time import Time
# import math
# import calendar
# import numpy as np
# from scipy.optimize import *
# from scipy.stats import chi2
# from PyAstronomy import pyasl
# import datareduc

folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)

[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)


# Hier was je gebleven

matplotlib.style.use('classic')
vsini = 127
c_light = 299792.458
# ['He_I', 2 ,6678],['CIII', 4647 +4651]
# ll2 = [['H_alpha', 2, 6562.819, 6554, 6556, 6578, 6579], ['H_beta',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I',30, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I',10, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He_II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He_II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O_III',14, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C_IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
# ll_laplma = [[r'H$\alpha$', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'H$\beta$', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'H$\gamma$', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4546.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4690.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
# ll2 = [[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6]]
ll2 = [[r'H$\beta$', 35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0], ['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll3 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$', 35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll4 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579],[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll_lapalma = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
ll_TVS_eshel = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
ll_nieuwe_lijnen = [['He_I', 35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, r'He I 6678']]

# ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, r'He I 6678']
fl_eshel_clean_folder = str(Data_folder)+r'\eShelData\data\clean' #Zonder twee spectra met rare Halpha spike
filelist_lapalma_folder = str(Data_folder)+r'\LaPalmaData'
fl_eshel_all_folder = str(Data_folder)+r'\eShelData\data\AlleSpectra'
fl_eshel_goodSNR_folder = str(Data_folder)+r'\eShelData\data'
fl_all = glob.glob(fl_eshel_all_folder+r'\*.fit')
fl_clean = glob.glob(fl_eshel_clean_folder+r'\*.fit')
fl_goodSNR = glob.glob(fl_eshel_goodSNR_folder+r'\*.fit')
filelist_lapalma = glob.glob(filelist_lapalma_folder+r'\*.fits')

datafile_folder_merc = str(converted_Data_folder)+r'\mercator\test\\'
datafile_folder_apo = str(converted_Data_folder)+r'\apo\test\\'
test_datafile_folder = str(converted_Data_folder)+r'\test\\'
def bjd(file):
    fits = pf.open(file)
    header = fits[0].header
    HJD =  float(header['BJD'])
    fits.close()
    return HJD


mark1 = [1, 37]
mark2 = [7, 10, 11, 32, 33, 34, 35, 36]
mark3 = [8, 9]
sortedfl_lapalma = sorted(filelist_lapalma, key=lambda x: bjd(x), reverse=False)
k = 1

#hier was je geb
for file in fl_all:
    mark = 0
    if k in mark1:
        mark = 1
    elif k in mark2:
        mark = 2
    elif k in mark3:
        mark = 3
    a = Datafile_class.Datafile_apo(file, i=k, mark=mark)
    dl, dl2 = airmass.split_date(a.header['DATE-OBS'])
    savename = datafile_folder_apo+a.observatory+'{num:02d}'.format(num=a.i)+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
    workfileresource = open(savename, 'wb')
    pickle.dump(a, workfileresource)
    workfileresource.close()
    print(k, mark)
    k+=1

quit()

# # # print sortedfl_lapalma
k=1
for file in sortedfl_lapalma:
    print(k)
    startdate = airmass.timeanddatelp(file)
    a = Datafile_class.Datafile_mercator(file, i=k)
    dl, dl2 = airmass.split_date(a.header['DATE-OBS'])
    savename = datafile_folder_merc+a.observatory+'{num:02d}'.format(num=a.i)+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
    pickle.dump(a, open(savename, 'w'))
    k+=1
# testfile_apo = fl_clean[12]
# a = Datafile_class.Datafile_apo(testfile_apo)
# testfile_merc = filelist_lapalma[3]
# a = Datafile_class.Datafile_mercator(testfile_merc, i=1)
# dl,dl2 = airmass.split_date(a.header['DATE-OBS'])
# testname = a.observatory+'{num:02d}'.format(num=a.i)+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
# print testname
#
# print a.i
# print a.phase
# print a.snr
# print a.line6562.ew
# for item in a.header:
#     print item
# print a.header['DATE-OBS']
# dl,dl2 = airmass.split_date(a.header['DATE-OBS'])
# print dl,dl2
# plt.plot(a.wl_original,a.flux_original,label = 'original')
# plt.plot(a.wl_rebin,a.flux_rebin, label = 'rebinned')
# plt.legend()
# plt.show()

# a = Datafile_class.Datafile_mercator()
# # make_line_data_dict(a.linelist,a.wl_rebin,a.flux_rebin,a.observatory)
# print a.filename
# savename = datafile_folder+a.filename+'.txt'
# pickle.dump(a, open(savename,'w'))

# line = ll_lapalma[3]
# wl = a.wl_rebin
# flux = a.flux_rebin
# obs = a.observatory
# line2,lw,lf,v,nf,vsini = line_data(line,wl,flux,obs)
# c = Line(line2,lw,lf,v,nf,vsini)


#
# b = pickle.load(datafile_folder+'test.txt','r')
# testfile = open(testsave,'r')
# b=pickle.load(testfile)
#
# print b.filename

# fl = glob.glob(r'D:\Peter\School\Master Thesis\Data\masterfiles\test\*.txt')
#
# testfile = open(fl[0],'r')
# b=pickle.load(testfile)
#
# print b.snr_original,b.snr

# plt.plot(b.line4026.v,b.line4026.flux)
# plt.show()
