from __future__ import division
import matplotlib.pyplot as plt
import glob
import warnings
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
from scipy import interpolate
from scipy.optimize import *
from PyAstronomy import pyasl
from SavitzkyGolay import savitzky_golay
import scipy.stats as ss
from pysynphot import observation
from pysynphot import spectrum
import Datafile_class
import airmass
c_light = 299792.458

fl_clean = glob.glob(r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\eShelData\data\clean\*.fit')
filelist = glob.glob(r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\LaPalmaData\*.fits')
filelist_lapalma = glob.glob(r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\LaPalmaData\*.fits')
filelist2 =  glob.glob(r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data/eShelData/data/*.fit')
filepath_eshel_spectra_info = r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\Data\masterfiles\dict_apo_files.txt'

def extractdata(j,file,header = 'off'):
    datafile = pf.open(file)
    data = datafile[j].data
    header = datafile[j].header
    currentwl = header['CRVAL1']
    step =header['CDELT1']
    wl=[]
    for i in range(len(data)):
        wl.append(currentwl)
        currentwl+=step
    wl2 = np.array(wl)
    if header == 'on':
        return wl2, data, header
    else:
        return wl2, data
def raw_apo_image():
    file =fl_clean[8]
    j=0
    datafile = pf.open(file)
    data = datafile[j].data
    header = datafile[j].header
    # print data
    # for file in fl_clean:
    # wl,flux,header = airmass.extractdata(j,file,header='on')
    for item in header:
        print(item)
    print(header['DATE-OBS'])
    # plt.plot(wl,flux)
    # plt.show()
    # plt.close()

    plt.imshow(data, cmap='gray',origin='lower')
    plt.title('Raw image of APO15')
    plt.axis('off')
    plt.savefig(r'E:\Peter\School\MasterThesis\2018-05-20\Master Thesis\figures\CompareApoMercator\raw_apo_image.pdf',format='pdf', dpi=1200)
    plt.show()

# file = filelist_lapalma[0]
# datafile = pf.open(file)
