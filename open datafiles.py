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
def open_masterfiles_apo(wantedmarks = [0,3],path = r'D:\Peter\Master Thesis\Data\masterfiles\apo\\',manual_filelist=None):
    if manual_filelist == None:
        fl = glob.glob(path+r'*.txt')
    else:
        fl=manual_filelist
    datafiles = []

    for file in fl:
        a = open(file, 'r')
        # for line in testfile:
        #     print line
        b = pickle.load(a)
        # print b.line6562.ew
        datafiles.append(b)
    newlist = [x for x in datafiles if x.mark in wantedmarks]
    sortednewlist = sorted(newlist,key=lambda x: x.i)
    return sortednewlist

def open_masterfiles_mercator(path = r'D:\Peter\Master Thesis\Data\masterfiles\mercator\\',manual_filelist=None):
    if manual_filelist == None:
        fl = glob.glob(path+r'*.txt')
    else:
        fl=manual_filelist
    datafiles = []

    for file in fl:
        a = open(file, 'r')
        # for line in testfile:
        #     print line
        b = pickle.load(a)
        # print b.line6562.ew
        datafiles.append(b)
    sortednewlist = sorted(datafiles,key=lambda x: x.i)
    return sortednewlist
