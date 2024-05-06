from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import glob
import SavitzkyGolay
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
import open_masterfiles
import warnings

files = open_masterfiles.apo()

clusters = []
cluster = []

i=1
imax = len(files)
print(imax)
while i<imax:
    file0 = files[i-1]
    file1 = files[i]
    t0 = file0.HJD
    t1 = file1.HJD
    dt = (t1-t0)*24*60
    # print i
    cluster.append(file0)
    if dt>20:
        clusters.append(cluster)
        cluster = []
    i+=1
cluster.append(file1)
clusters.append(cluster)

# for file in files:
for file in files:
    print(file.i, file.time_and_date)
#
# # k=1
# a = open(file, 'r')
# clusters = pickle.load(open(r'D:\Peter\School\Master Thesis\Data\apo_file_stack_list.txt','r'))
# for cluster in clusters:
#     print '----'
#     wls = []
#     for file in cluster:
#         wl =file.wl_rebin
#         flux =file.flux_rebin
#         wls.append(wl)
#         print type(wl), type(flux)
#         # k+=1

cluster = clusters[2]
HJDs = []
fluxs = []
for file in cluster:
    HJDs.append(file.HJD)
    fluxs.append(file.flux_rebin)
    plt.plot(file.wl_original,file.flux_original, label = 'APO' + str(file.i))

# wl = cluster[0].wl_rebin
# flux = np.sum(np.array(fluxs),axis=0)/len(cluster)

# plt.plot(wl,flux,label='avg')
# plt.legend()
plt.legend()

plt.ylabel('Relative Flux')
plt.xlabel('Wavelength ($\AA$)')
plt.show()
plt.close()
# a = [1,2,3,4,5,6]
# b= [1,2,3,4,5,6]
# c=[1,2,3,4,5,6]
# d = [1,2,3,4,5,6]
#
# array = np.array([a,b,c,d])
# sumarray = np.sum(array, axis=0)
# print sumarray


# hjdavg = np.mean(HJDs)
# print HJDs
# print hjdavg
#
# testinstance = Datafile_class.datafile_stack_apo(cluster)
# print testinstance.HJD
# pickle.dump(clusters, open(r'D:\Peter\School\Master Thesis\Data\apo_file_stack_list.txt','w'))

# for i,file in enumerate(files[:-1]):
#     cluster.append(file)
#     t1 = file.HJD
#     t2 = files[i+1].HJD
#     dt = (t2-t1)*24*60
#     print dt
#     if dt >20:
#         print '----'
#         clusters.append(cluster)
#         cluster = []
# clusters[-1].append(files[-1])
# for a in clusters:
#     print '----'
#     HJD0=a[0].HJD
#     for file in a:
#         print (file.HJD-HJD0)*24*60