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
print 'asdf'
nrow = 3
ncol = 2

def my_sin(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*x  + phase)  + offset

def my_line(x,a,b):
    return a*x+b

def flat_line(x,a):
    return a


datafolder = r'C:\peter\School\Master Scriptie\Data\EW\apo/'

datafile_list = glob.glob(datafolder + '*[!nn].npy')
sorted_datafilelist = sorted(datafile_list, key=lambda x: np.load(x)[9])
print(len(sorted_datafilelist))
print datafile_list
print sorted_datafilelist
# print len(datafile_list)
fig, axs = plt.subplots(nrows=nrow, ncols=ncol,sharex=True,figsize=(10,9))
# fig.suptitle('Normalized Equivalent Width', size=16)
k=0
for i, row in enumerate(axs):
    for j, ax in enumerate(row):
        file1 = np.load(sorted_datafilelist[k%len(sorted_datafilelist)])
        k+=1
        ew = file1[0]
        phase = file1[1]
        error = file1[2]
        lineused = file1[3][0]
        rel_llh = file1[9]
        print lineused
        linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
        # linename.replace('_',' ')
        print linename
        p1=[0.005,0.5, 1]
        fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
        t= np.linspace(0,1,50)
        data_fit = my_sin(t, *fit1[0])
        ax.set_title(linename)
        ax.errorbar(phase, ew, yerr=error, fmt='o', label = 'HERMES', c='r')
        ax.plot(t,data_fit, label='Fit HERMES',c='r')
        aldat = np.append(ew, data_fit)
        ax.locator_params(axis='y',nbins=7)
        lowlim,uplim =np.round(np.min(aldat),decimals=2)-0.01,np.round(np.max(aldat),decimals=2)+0.01
        ax.text(0.8, 0.8*(uplim-lowlim)+lowlim,str(int(round(1/rel_llh,0)) ))
        ax.set_ylim(uplim,lowlim)
plt.tight_layout()
fig.get_axes()[0].annotate(r'Normalized Equivalent Width', (0.5, 0.972),xycoords='figure fraction', ha='center',fontsize=16)
fig.text(0.25, 0.01, 'Phase (6.83d)', ha='center', va='center',size=14)
fig.text(0.74, 0.01, 'Phase (6.83d)', ha='center', va='center',size=14)
plt.show()
plt.close()



#
# f = open(datafolder + 'chisquared_manual.txt', 'w')
# f.write('line name'+' \t '+'chi2 sin' +'\n' )
# for file in datafile_list:
#     file1 =np.load(file)
#     ew = file1[0]
#     print ew
#     phase = file1[1]
#     error = file1[2]
#     lineused = file1[3][0]
#     # print lineused
#     linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
#     p1=[0.005,0.5, 1]
#     fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
#     p2 = [0,1]
#     fit2 = curve_fit(my_line, phase, ew, p0=p2,  sigma=error, absolute_sigma=True)
#     p3 =[1]
#     fit3 = curve_fit(flat_line , phase, ew, p0=p3,  sigma=error, absolute_sigma=True)
#     # t= np.linspace(0,1,50)
#     data_fit_sin = my_sin(np.array(phase), *fit1[0])
#     data_fit_line = my_line(np.array(phase), *fit2[0])
#     data_fit_flat = flat_line(np.array(phase), *fit3[0])
#
#     residuals = (ew - data_fit_sin)/error
#     chisq = np.sum(residuals**2) / (len(ew)-1-3)
#     chisq2 = airmass.redchisqg(ew,data_fit_sin,deg=3)
#
#     print chisq, chisq2, chisq-chisq2
#     chi2_sin,p_sin = ss.chisquare(ew,data_fit_sin,ddof=3)
#     chi2_line,p_line = ss.chisquare(ew,data_fit_line,ddof=2)
#     chi2_flat,p_flat = ss.chisquare(ew,data_fit_flat,ddof=1)
#     dif1 = p_sin-p_line
#     dif2 = p_sin-p_flat
#     line_to_write = linename+' \t '+ str(chisq)+ '\n'
#     # print line_to_write
#     f.write(line_to_write)
#     print linename
#     print ew
#     print data_fit_sin
#     print ew-data_fit_sin
#     print error
# #






    # print ew
    # print len(ew), len(phase),len(error)



