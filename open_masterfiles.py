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
import pickle
import os
# import Datafile_class
import sys
from Datafile_class import *
sys.modules['Line'] = Line
sys.modules['Datafile_mercator'] = Datafile_mercator
sys.modules['Datafile_apo'] = Datafile_apo
def apo(path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\apo\test\\', wantedmarks = None,manual_filelist=None):
    if manual_filelist == None:
        fl = glob.glob(path+r'*.txt')
    else:
        fl=manual_filelist
    datafiles = []
    # for line in testfile:
    #     print line
    # print b.line6562.ew
    for file in fl:
        a = open(file, 'rb')
        # for line in testfile:
        #     print line
        b = pickle.load(a)
        # print b.line6562.ew
        datafiles.append(b)
        a.close()
    if wantedmarks is None:
        newlist=datafiles
    else:
        newlist = [x for x in datafiles if x.mark in wantedmarks]
    sortednewlist = sorted(newlist,key=lambda x: x.i)
    return sortednewlist
def mercator(path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\mercator\test\\',manual_filelist=None):
    if manual_filelist == None:
        fl = glob.glob(path+r'*.txt')
    else:
        fl=manual_filelist
    datafiles = []

    for file in fl:
        a = open(file, 'rb')
        # for line in testfile:
        #     print line
        b = pickle.load(a)
        # print b.line6562.ew
        datafiles.append(b)
        a.close()
    sortednewlist = sorted(datafiles,key=lambda x: float(x.HJD))
    return sortednewlist

def apo_demetra(path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\test\demetra\\',manual_filelist=None,sorted='off'):
    if manual_filelist == None:
        fl = glob.glob(path+r'*.txt')
    else:
        fl=manual_filelist
    datafiles = []

    for file in fl:
        a = open(file, 'rb')
        # for line in testfile:
        #     print line
        b = pickle.load(a)
        # print b.line6562.ew
        datafiles.append(b)
        a.close()
    if sorted == 'on':
        sortednewlist = sorted(datafiles,key=lambda x: x.i)
        return sortednewlist
    else:
        return datafiles

def apo_demetra_orders(path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\\',manual_filelist=None,sort_data_files='off'):
    if manual_filelist == None:
        fl = glob.glob(path+r'*.txt')
    else:
        fl=manual_filelist
    datafiles = []

    for file in fl:
        a = open(file, 'rb')
        # for line in testfile:
        #     print line
        b = pickle.load(a)
        # print b.line6562.ew
        datafiles.append(b)
        a.close()
    # print(datafiles[0].header['DATE-OBS'])
    if sort_data_files == 'on':
        sortednewlist = sorted(datafiles,key=lambda x: x.header['JD-MID'])
        return sortednewlist
    else:
        return datafiles





def open_linelist(path):
    a = open(path, 'rb')
    b = pickle.load(a)
    a.close()
    return b