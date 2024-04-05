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

c_light = 299792.458
# ['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],
linelist_apo = [ ['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2],  ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
linelist = [4861.01,4471.47,6678.15,4921.93,4647.42,5015.68,5411.52,4650.85,4541.52,5592.25]
ll2= [['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0],['H_gamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0],['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6],['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0],['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914],['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471],['He_II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],['C_IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
filelist_lapalma = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
filelist_apo =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
nall = [['Na doublet',5892.95,5888.5,5889,5894.2,5894.75]]
ll3 = [[r'H$\beta$', 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I', 5875.621, 5865, 5869.1, 5879.9, 5882]]
# ll_laplma = [[r'Halpha', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'Hbeta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hgamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I', 5875.621, 5860.6241285010492, 5861.5995122694139, 5860.6241285010492, 5870.5989810949914], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4681.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5806.33334373, 5808.314546]]
ll_lapalma = [[r'Ha', 6562.819, 6551, 6552, 6578, 6579], [r'Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hy', 4340.472, 4322, 4324, 4357, 4360], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8], ['He_II', 4541.6, 4498, 4499, 4580, 4581], ['He_II', 4685.804, 4679, 4680, 4690, 4691], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1]]

def my_sin2(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*(x  + phase))  + offset

# datafolder = r'C:\peter\School\Master Scriptie\Data\EW\apo/'
datafolder = r'C:\peter\School\Master Scriptie\Data\EW\lp/'
datafile_list = glob.glob(datafolder + '*_nn.npy')
dat0 = np.load(datafile_list[0])
print dat0[3]
phases0 = np.array(dat0[1])
# print phases0
print 'nr', 'phase    ',
for file in datafile_list:
    dat = np.load(file)
    # print 'asdf', dat[3]
    linename = dat[3][0][0]+' ' + str(int(dat[3][0][1]))+'            '
    print linename,
print
# print 'asdf'
for i in range(len(phases0)):
    print i+1,' ', round(phases0[i],3),' ',
    for file in datafile_list:
        dat = np.load(file)
        ew = round(dat[0][i],3)
        er = round(dat[2][i],3)
        print ew,'$\pm$',er,'   ',
    print
# nrcolumns = 8
# print r'\begin{tabular}{c'+'|c'*(nrcolumns)+r'}\\'
# print r'Line & \multicolumn{'+str(nrcolumns)+r'}{|c|}{Equivalent width}\\'
# print r'\hline'
# print r'\hline'
# print 'Phase',
# inds = phases0.argsort()
# print inds
# for phase in phases0[inds][:8]:
#     print '&',round(phase,3),
# print r'\\'
# print r'\hline'
#
#
# for line in ll_lapalma:
#     linename = linename = line[0]+str(int(line[1]))
#
#     linename2 = linename[:-4]+' '+linename[-4:]
#     linename3 =linename2.replace('_', ' ')
#     print linename3,
#     dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lp\\'+linename+'.npy')
#     # dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lp\\'+linename+'.npy')
#     # lapalmahe5875.npy
#     ew1 = np.array(dat1[0])
#     eb1 = np.array(dat1[2])
#     phases1 = np.array(dat1[1])
#     ewsorted = ew1[inds]
#     ebsorted = eb1[inds]
#     phasessored = phases1[inds]
#     for i, ew in enumerate(ewsorted[:8]):
#         print '&',round(ew,2),'$\pm$',round(ebsorted[i],2),
#
#     p1=[0.005,0.5, 1]
#
#     fit1 = curve_fit(my_sin2, phases1, ew1, p0=p1,  sigma=eb1, absolute_sigma=True)
#     print r' \\'
#     # print '&',fit1[0][1], '$\pm$',fit1[1][1][1],r'\\'
# print '\end{tabular}'
#
#
#
#
#
# print r'\begin{tabular}{c'+'|c'*(nrcolumns+1)+r'}\\'
# print r'Line & \multicolumn{'+str(nrcolumns)+r'}{|c|}{Equivalent width}&\phi_{0, best fit}\\'
# print r'\hline'
# print r'\hline'
# inds = phases0.argsort()
# print 'Phase',
# for phase in phases0[inds][8:]:
#     print '&',round(phase,3),
# print r'& \\'
# print r'\hline'
#
# for line in ll_lapalma:
#     linename = line[0]+str(int(line[1]))
#
#     linename2 = linename[:-4]+' '+linename[-4:]
#     linename3 =linename2.replace('_', ' ')
#     print linename3,
#     dat1 = np.load('C:\peter\School\Master Scriptie\Data\EW\lp\\'+linename+'.npy')
#     ew1 = np.array(dat1[0])
#     eb1 = np.array(dat1[2])
#     phases1 = np.array(dat1[1])
#     # print '********', len(phases1), '********'
#     ewsorted = ew1[inds]
#     ebsorted = eb1[inds]
#     phasessored = phases1[inds]
#     for i, ew in enumerate(ewsorted[8:]):
#         print '&',round(ew,2),'$\pm$',round(ebsorted[i],2),
#
#     p1=[0.005,0.5, 1]
#
#     # fit1 = curve_fit(my_sin2, phases1, ew1, p0=p1,  sigma=eb1, absolute_sigma=True)
#     print '&',round(fit1[0][1],3), '$\pm$',round(fit1[1][1][1],3),r'\\'
#     # print r'\\'
# print '\end{tabular}'
#
#
#
#
