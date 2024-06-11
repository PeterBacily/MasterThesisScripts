from __future__ import division
import matplotlib.pyplot as plt
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
import os
import Datafile_class
matplotlib.style.use('classic')
import open_masterfiles


# merc_linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
# apo_linelist = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]

merc_linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4016, 4020, 4032, 4036, 'He I 4026'], ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'], ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5884.6, 5885.5, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
apo_linelist = [['Ha', 35, 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]

def normalize(wave,flux,a,b,c,d,startwl,endwl):
    normwave = np.hstack((wave[(wave>a)&(wave<b)],wave[(wave>c)&(wave<d)]))
    normflux = np.hstack((flux[(wave>a)&(wave<b)],flux[(wave>c)&(wave<d)]))
    #fit line trough slice
    slope,height = np.polyfit(normwave,normflux,1)
    # print 'slope and height are', slope, height
    fit = np.poly1d([slope,height])
    linewave = wave[(wave>startwl)&(wave<endwl)]
    # print wave[0],startwl,wave[-1],endwl
    # print wave
    lineflux = flux[(wave>(startwl))&(wave<(endwl))]
    normlineflux = []
    # print linewave
    for i,j in enumerate(linewave):
        normlineflux.append(lineflux[i]/fit(j))
    fluxarray = np.array(normlineflux)
    nnf = []
    for k, nwl in enumerate(normwave):
        nnf.append(normflux[k]/fit(nwl))
    return linewave, fluxarray, nnf, lineflux, fit
print('1')
apo_files = open_masterfiles.apo()
# print apo_files
print('2')
merc_files = open_masterfiles.mercator()
print('3')
merc_file =apo_files[1]
# line = ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861']
for i,line in enumerate(apo_linelist):
    line_oud = apo_linelist[i]
    k = 2
    start = line[k+1]-20
    stop = line[k+4]+20
    wl = merc_file.wl_rebin
    flux = merc_file.flux_rebin

    lw,nlf,nnf,lf,fit = normalize(wl,flux,line[k+1],line[k+2],line[k+3],line[k+4],line[k+1]-20,line[k+4]+20)
    cond1 = (lw>line[k+1])&(lw<line[k+2])
    cond2 = (lw>line[k+3])&(lw<line[k+4])
    normwave1 = lw[cond1]
    normflux1= lf[cond1]
    normwave2 = lw[cond2]
    normflux2= lf[cond2]


    # lw_oud,nlf_oud,nnf_oud,lf_oud,fit_oud = normalize(wl,flux,line_oud[k+1],line_oud[k+2],line_oud[k+3],line_oud[k+4],line_oud[k+1]-20,line_oud[k+4]+20)
    # cond1_oud = (lw>line_oud[k+1])&(lw<line_oud[k+2])
    # cond2_oud = (lw>line_oud[k+3])&(lw<line_oud[k+4])
    # normwave1_oud = lw[cond1_oud]
    # normflux1_oud= lf[cond1_oud]
    # normwave2_oud = lw[cond2_oud]
    # normflux2_oud= lf[cond2_oud]



    # ((lw>line[k+3])&(lw<line[k+4]))]

    normflux3 = lf[lw<line[k+1]]
    normwave3 = lw[lw<line[k+1]]

    normflux4 = lf[(lw>line[k+2])&(lw<line[k+3])]
    normwave4 = lw[(lw>line[k+2])&(lw<line[k+3])]

    normflux5 = lf[lw>line[k+4]]
    normwave5 = lw[lw>line[k+4]]
    # normflux = lf[((lw>line[k+1])&(lw<line[k+2])) or ((lw>line[k+3])&(lw<line[k+4]))]
    f, (ax1, ax2) = plt.subplots(2, sharex=True)
    # ax1.plot(lw,lf)
    ax1.plot(normwave1,normflux1, c='r', label = 'continuum data slice')
    ax1.plot(normwave2,normflux2, c='r')
    ax1.plot(normwave3,normflux3, c='b', label = 'data')
    ax1.plot(normwave4,normflux4, c='b')
    ax1.plot(normwave5,normflux5, c='b')
    # ax1.plot(normwave1_oud,normflux1_oud,c='g')
    # ax1.plot(normwave2_oud, normflux2_oud, c='g')
    # ax1.plot(lw, fit_oud(lw), c='g', label='fitted continuum')
    ax1.plot(lw,fit(lw),c='red', linestyle='--', label = 'fitted continuum')
    ax1.legend(prop={'size': 10}, loc='best')
    # ax2.plot
    ax2.plot(lw,nlf,c='r')
    ax1.set_title('Raw Flux')
    ax2.set_title('Normalized Flux')
    # ax2.plot(lw_oud,nlf_oud,c='g')
    plt.suptitle(line[-1])
    plt.show()
    plt.close()




