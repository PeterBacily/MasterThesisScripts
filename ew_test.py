from __future__ import division
import matplotlib.pyplot as plt
import glob
import SavitzkyGolay
import pyfits as pf
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
import open_masterfiles
import warnings
def chunks(l, n):
    # For item i in a range that is a length of l,
    chunklist = []
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        chunk = l[i:i+n]
        chunklist.append(chunk)
    return chunklist

def my_sin(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*x  + phase)  + offset

# apo_lines  = ['line6562', 'line4861', 'line4713', 'line5875', 'line4541', 'line4685', 'line5411', 'line4471', 'line4921', 'line6678', 'line5592', 'line5801']
# mercator_lines = ['line6562', 'line4861', 'line4340', 'line4026', 'line4471', 'line4713', 'line5875', 'line4541', 'line4685', 'line5411', 'line4921', 'line6678', 'line5592', 'line5801']
def plot_EW(obs='MERCATOR',apowantedmarks = None):
    apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']
    mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471', 'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
    custom_lines = []
    custom_files = []

    figsavefolder = r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\figures\EWs\verslag\\'
    # split_up_lines = list(chunks(mercator_lines,3))
    # print chunks(apo_lines,3)
    # obs = 'APO'
    # obs = 'MERCATOR'
    if obs =='APO':
        if apowantedmarks == None:
            master_files = open_masterfiles.apo()
        else:
            master_files = open_masterfiles.apo(wantedmarks=apowantedmarks)
        lines = apo_lines
    elif obs == 'MERCATOR':
        master_files = open_masterfiles.mercator()
        lines =mercator_lines
    elif obs =='BOTH':
        apo_files = open_masterfiles.apo()
        mercator_files = open_masterfiles.mercator()
        lines = apo_lines
    else:
        master_files = custom_files
        lines = custom_lines
        warnings.warn('observatory is not APO, MERCATOR or BOTH, if this is not intentional look at obs variable')

    if obs is not 'BOTH':
        for line in lines:
            # fig, axs = plt.subplots(nrows=len(chunk), ncols=1,sharex=True,figsize=(5, 8))
            # fig.subplots_adjust(top = 0.83)
            # size = fig.get_size_inches()
            # savename = ''
            phases = []
            ews = []
            ew_error = []
            for file in master_files:
                linedata = getattr(file,line)
                lineinfo = linedata.lineinfo
                phases.append(file.phase)
                ews.append(linedata.ew)
                ew_error.append(linedata.ew_error)
                er, ew, wl_linepart, v_linepart, flux_linepart = airmass.testew(linedata.v_cor,linedata.wl,linedata.flux,lineinfo[2],file.snr)
                print line, linedata.ew-ew
                plt.plot(v_linepart,flux_linepart)
            plt.show()
            #
            #     chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats2(np.array(ews),np.array(phases),np.array(ew_error))
            #     print lineinfo
            #     p1 = [0.005, 0.5, 1]
            #     fit1 = curve_fit(my_sin, phases, ews, p0=p1, sigma=ew_error, absolute_sigma=True)
            #     fit2 = curve_fit(airmass.flat_line, phases,ews,p0=[1],sigma=ew_error,absolute_sigma=True)
            #     t = np.linspace(0, 1, 50)
            #     data_fit = my_sin(t, *fit1[0])
            #     ax.set_title(lineinfo[-1])
            #     l1 = ax.errorbar(phases, ews, yerr=ew_error, fmt='o', c='r', label='Data')
            #     l2, = ax.plot(t, data_fit, c='r', label='Sinusodial fit')
            #     l3 = ax.axhline(fit2[0][0],linestyle = '--',c='gray', label = 'Weighted Average')
            #     aldat = np.append(ews, data_fit)
            #     ax.locator_params(axis='y',nbins=7)
            #     lowlim,uplim =np.round(np.min(aldat),decimals=2)-0.01,np.round(np.max(aldat),decimals=2)+0.01
            #     l4 = ax.text(0.75, 0.95*(uplim-lowlim)+lowlim,str('%.2E' % (1./probfactor)),label = 'likelihood ratio')
            #     ax.set_ylim(uplim,lowlim)
            #     savename+= lineinfo[0]+line[-4:]+'_'
            # savename += obs+'_EW.pdf'
            # print savename
            # plt.figlegend((l1,l2,l3),('Data','Sinusodial fit','Weighted Average'),bbox_to_anchor = [1., 0.92], prop={'size': 10},fancybox=True, framealpha=1)
            # plt.suptitle('Normalized Equivalent width \n'+obs, size=16)
            # # plt.tight_layout()
            # plt.savefig(figsavefolder+savename)
            # # plt.show()
            # plt.close()
plot_EW(obs='APO')