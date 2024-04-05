from __future__ import division
import glob
import pyfits as pf
import numpy as np
import os
import csv
import matplotlib.pyplot as plt

#set working directory
workingdirectory = r'C:/peter/School/Master Scriptie/Data'
filelist = glob.glob(workingdirectory + "/*.fit")
print filelist

# select and open a file from filelist
datafile = pf.open(filelist[0])

#call the header of the raw image
header = datafile[0].header


# j=0 for raw image j=1.....j=34 for individual orders, j=35 for all orders merged and j=36 for datareduction parameters
# if j=36 comment out the wl array making part
j=2
hd2 = datafile[j].header
dt2 = datafile[j].data
print dt2
#print the header
# for item in hd2:
#     print item
# print hd2
# # make arrway of Wavelengths
currentwl = hd2['CRVAL1']
step =hd2['CDELT1']
wl=[]
for i in range(len(dt2)):
    wl.append(currentwl)
    currentwl+=step
#
# # plot
# plt.plot(wl,dt2)
# plt.xlabel('Wavelength [A]')
# plt.ylabel('Counts')
# plt.show()
print type(wl)
wl2 = np.array(wl)
print type(dt2), type(wl2)
print len(dt2)
test = dt2[(wl2>6600) & (wl2<6650)]
print len(test)
x= wl2[(wl2>6600) & (wl2<6650)]
print len(x)
# print wl2>6600, wl2<6650, (wl2>6600) & (wl2<6650)
# x = wl[(wl>6600) & (wl<6650)]
# print len(x), len(test)
plt.plot(x,test)
plt.show()