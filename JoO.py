from __future__ import division
import pyfits as pf
import glob
import airmass
# import datareduc
import matplotlib.pyplot as plt
import numpy as np
import Datafile_class
import open_masterfiles
import pickle
c_light = 299792.458

# filelist = glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
# print filelist
# asdf = pf.open(filelist[12])
# newlist = sorted(filelist, key=lambda x: airmass.airmass(x)[2])
# print airmass.snr(filelist[12])
obs = 'MERC'

if obs == 'APO':
    files = open_masterfiles.apo(wantedmarks=[0, 1, 2, 3])
    print r'\begin{tabular}{ l|| r| r| r| r|r|r|r|r|r|l }'
    print r'\# & Date & HJD & T$_{\textrm{exp}}$  & SNR & Airmass & Alt & Phase & BC& v$_{\textrm{ISM}}$ & Notes\\'
    print r' APO & 2016 & $-$2457000 & (s)& & & (deg)& p = 6.83 d& (km/s) & (km/s) &          \\'
    print r'\hline'
    files.sort(key=lambda x: x.i)
    for file in files:
        if file.mark == 0:
            note = ''
        else:
            note = str(file.mark)
        print obs, file.i, '&', file.time_and_date, '&', "{:.3f}".format(file.HJD - 2457000), '&', "{:.0f}".format(
            file.exptime), '&', "{:.0f}".format(file.snr), '&', "{:.1f}".format(file.airmass), '&', "{:.0f}".format(
            file.alt), '&', "{:.3f}".format(file.phase), '&', "{:.0f}".format(
            file.baricentric_correction), '&', "{:.0f}".format(file.velshift), '&', note, r'\\'


elif obs == 'MERC':
    files = open_masterfiles.mercator()
    print r'\begin{tabular}{ l|| r| r| r| r|r|r|r|r|r }'
    print r'\# & Date & HJD & T$_{\textrm{exp}}$  & SNR & Airmass & Alt & Phase & BC& v$_{\textrm{ISM}}$ \\'
    print r' & 2016 & $-$2457000 & (s)& & & (deg)& p = 6.83 d& (km/s) & (km/s) \\'
    print r'\hline'
    files.sort(key=lambda x: x.i)
    for file in files:
        print obs, file.i, '&', file.time_and_date, '&', "{:.3f}".format(file.HJD - 2457000), '&', "{:.0f}".format(
            file.exptime), '&', "{:.0f}".format(file.snr), '&', "{:.1f}".format(file.airmass), '&', "{:.0f}".format(
            file.altitude), '&', "{:.3f}".format(file.phase), '&', "{:.0f}".format(
            file.baricentric_correction), '&', "{:.0f}".format(file.velshift), r'\\'


# snrs = []
# airmasses = []
# for file in newlist[1:-1]:
#     i+=1
#     fwl = airmass.fitfraun(file)
#     velshift = c_light*(fwl-5895.92)/5895.92
#     SNR = round(airmass.snr(file),1)
#     date = airmass.timeanddate(file)[0]
#     expt = int(round(airmass.exposuretime(file)))
#     am = round(airmass.airmass(file)[0],1)
#
#     alt = round(round(airmass.airmass(file)[1],1))
#     # JD = airmass.airmass(file)[2]
#     BCCor,HJD = airmass.barcor(file)
#     phase = round(airmass.aphase(HJD),3)
#     airmasses.append(am)
#     snrs.append(SNR)
#     print str(i), '&' ,str(date),'&' ,str(round(HJD-2457000,3)), '&' ,str(expt), '&' ,str(int(np.round(SNR))), '&' ,str(am),'&' ,str(alt), '&' ,str(phase),'&' ,str(int(np.round(BCCor))),'&', str(int(np.round(velshift))),r'\\'
# print r'\end{tabular}'
#


# plt.title('Slice of data used to calculate SNR')
# plt.xlabel('Wavelength ($\AA$)')
# plt.ylabel('Normalized Flux')
# plt.show()
# plt.close()
#
# plt.scatter(airmasses,snrs)
# plt.title('SNR as a function of airmass')
# plt.xlabel('Airmass')
# plt.ylabel('snr')
# plt.show()
# plt.close()