from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
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

apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
             'line4861', 'line4921', 'line6678', 'line4471']
apo_master_files = open_masterfiles.apo()

merc_master_files = open_masterfiles.mercator()

obs = 'MERCATOR'



if obs == 'APO':
    master_files = apo_master_files

elif obs == 'MERCATOR':
    master_files = merc_master_files
number_of_lines_1 = 5
print r'\begin{table*}'
print r'\caption{',obs,' 1 }'
print r'\label{}'
print r'\begin{tabular}{l r c',
for i in range(number_of_lines_1):
    print ' c',
print r'}'
print r'\hline'
print r'\hline'
print r'BJD-2457300 & \# & Phase',
for line in apo_lines[:number_of_lines_1]:
    testfile = master_files[-5]
    testlinedata = getattr(testfile, line)
    testlineinfo = testlinedata.lineinfo
    print ' & ',testlineinfo[-1],
print r'\\'

print '\hline'

for file in master_files:
    print '%.3f' %(file.HJD-2457300),'&',file.i,'&', '%.2f' % file.phase,
    for line in apo_lines[:number_of_lines_1]:
        linedata = getattr(file, line)
        lineinfo = linedata.lineinfo
        print ' & ', '%.2f'%linedata.ew + r'$\pm$'+ '%.2f'%linedata.ew_error,
    print r'\\'
print r'\hline'
print r'\end{tabular}'
print r'\end{table*}'

number_of_lines_2 = len(apo_lines)-number_of_lines_1

print r'\begin{table*}'
print r'\caption{',obs,' 2 }'
print r'\label{}'
print r'\begin{tabular}{r c',
for i in range(number_of_lines_2):
    print ' c',
print r'}'
print r'\hline'
print r'\hline'
print r'\# & Phase',
for line in apo_lines[number_of_lines_1:]:
    testfile = master_files[-5]
    testlinedata = getattr(testfile, line)
    testlineinfo = testlinedata.lineinfo
    print ' & ',testlineinfo[-1],
print r'\\'

print '\hline'

for file in master_files:
    print file.i,'&', '%.2f' % file.phase,
    for line in apo_lines[number_of_lines_1:]:
        linedata = getattr(file, line)
        lineinfo = linedata.lineinfo
        print ' & ', '%.2f'%linedata.ew + r'$\pm$'+ '%.2f'%linedata.ew_error,
    print r'\\'
print r'\hline'
print r'\end{tabular}'
print r'\end{table*}'


#     ews.append(linedata.ew)
#     ew_error.append(linedata.ew_error)