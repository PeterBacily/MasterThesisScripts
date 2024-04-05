from __future__ import division
import pyfits as pf
import glob
import airmass
# import datareduc
import matplotlib.pyplot as plt
import numpy as np
c_light = 299792.458

filelist = glob.glob(r'D:\Peter\Master Thesis\Data\eShelData\data\AlleSpectra/*.fit')
print filelist
asdf = pf.open(filelist[12])
newlist = sorted(filelist, key=lambda x: airmass.airmass(x)[2])
# print airmass.snr(filelist[12])

print r'\begin{tabular}{ l|| r| r| r| r|r|r|r|r|r }'
print r'\# & Date & HJD & T$_{\textrm{exp}}$  & SNR & Airmass & Alt & Phase & BC&  v$_{\textrm{ISM}}$\\'
print r' - & 2016 & $-$2457000 & (s)& & & (deg)& p = 6.83 d& (km/s) & (km/s)'
i = 0
# for file in filelist:
snrs = []
airmasses = []
phases = []
dict = {}
dict['header'] = ['i','Telescope','date','HJD-2457000','T_exp','SNR','Airmass','Altitude','phase','BCcor','vshift']
for file in newlist:
    i+=1
    try:
        fwl = airmass.fitfraun(file)
    except RuntimeError:
        fwl = 5895.92
    velshift = c_light*(fwl-5895.92)/5895.92
    SNR = round(airmass.snr(file, 4991,5001),1)
    date = airmass.timeanddate(file)[0]
    expt = int(round(airmass.exposuretime(file)))
    am = round(airmass.airmass(file)[0],1)

    alt = round(round(airmass.airmass(file)[1],1))
    # JD = airmass.airmass(file)[2]
    BCCor,HJD = airmass.barcor(file)
    # phase = round(airmass.aphase(HJD),3)
    ap = airmass.aphase(HJD)
    phase = format(airmass.aphase(HJD), '.3f')
    phases.append(phase)
    airmasses.append(am)
    snrs.append(SNR)
    HJD_rounded = round(HJD-2457000,3)
    dictentry = [i,'APO', str(date), str(round(HJD-2457000,3)), str(expt), str(SNR), str(am),str(int(alt)), str(ap), str(BCCor), str(velshift)]
    dict[HJD_rounded] = dictentry
    print str(i), '&' ,str(date),'&' ,str(round(HJD-2457000,3)), '&' ,str(expt), '&' ,str(int(np.round(SNR))), '&' ,str(am),'&' ,str(int(alt)), '&' ,str(phase),'&' ,str(int(np.round(BCCor))),'&', str(int(np.round(velshift))),r' & \\'
print r'\end{tabular}'
# plt.scatter(phases,np.zeros(len(phases)))
# plt.show()
filepath = r'D:\Peter\Master Thesis\Data\masterfiles\dict_apo_files.txt'
f = open(filepath,"w")
f.write( str(dict) )
f.close()

# print len(dict['header']), dict['header']
# print len(dict[HJD_rounded]), dict[HJD_rounded]
# y = [1]*len(phases)
# plt.scatter(phases,y)
# plt.xlim(0,1)
# plt.show()
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