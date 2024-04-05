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
from scipy.stats import chi2
from PyAstronomy import pyasl
import matplotlib.style
import datareduc
import pickle
import os
import Datafile_class
import open_masterfiles
matplotlib.style.use('classic')
vsini =127
c_light = 299792.458



fl_eshel_clean_folder = r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\eShelData\data\clean' #Zonder twee spectra met rare Halpha spike
filelist_lapalma_folder = r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\LaPalmaData'
fl_eshel_all_folder = r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\eShelData\data\AlleSpectra'
fl_eshel_goodSNR_folder = r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\eShelData\data'
fl_all = glob.glob(fl_eshel_all_folder+r'\*.fit')
fl_clean = glob.glob(fl_eshel_clean_folder+r'\*.fit')
fl_goodSNR = glob.glob(fl_eshel_goodSNR_folder+r'\*.fit')
filelist_lapalma = glob.glob(filelist_lapalma_folder+r'\*.fits')

goodfile = fl_all[25]

wl, data, header = airmass.extractdata(20,goodfile,return_header='on')
wl_full,data_full,header_full = airmass.extractdata(35,goodfile,return_header='on')


lp_file = open_masterfiles.mercator()[4]
lp_wl_long,lp_flux_long = lp_file.wl_rebin, lp_file.flux_rebin

start = wl[0]
stop = wl[-1]
wl2 = wl_full[(wl_full>start) &(wl_full< stop)]
data2 = data_full[(wl_full>start) &(wl_full< stop)]
data2_norm = data2/np.mean(data2)

lp_wl= lp_wl_long[(lp_wl_long>start) &(lp_wl_long< stop)]
lp_flux= lp_flux_long[(lp_wl_long>start) &(lp_wl_long< stop)]
lp_flux_norm= lp_flux/np.mean(lp_flux)

f, (ax1, ax2) = plt.subplots(2, sharex=True)
#
ax1.plot(wl,data, label = 'Single Order')
ax1.plot(wl2,data2_norm,label='Order Merged')
ax2.plot(lp_wl,lp_flux_norm)
ax1.set_title('eShel Spectrum')
ax2.set_title('HERMES Spectrum')
ax2.set_xlabel('Wavelength ($\AA$)')
ax1.set_ylabel('Relative Flux')
ax2.set_ylabel('Relative Flux')
ax1.legend(loc=0)
plt.show()