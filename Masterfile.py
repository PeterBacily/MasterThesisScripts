from __future__ import division
import pyfits as pf
import glob
import airmass
# import datareduc
import matplotlib.pyplot as plt
import numpy as np
c_light = 299792.458
ll = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579],[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]

filelist = glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/*.fit')
legend = [['filename','date','HJD','SNR','T_exp','WL_ism','v_ism','airmass','altitude','phase','BCCor'],\
          ['wl (rebinned, baricentric corrected, radial velocity shifted)','v','flux' ]\
          ]
data = []
linelegend = []
for line in ll:
    linelegend.append([line[0]+' lw',line[0]+' v',line[0]+' lf',line])
legend.append(linelegend)

for file in filelist:
    filedata = []
    fwl = airmass.fitfraun(file)
    velshift = c_light*(fwl-5895.92)/5895.92
    SNR = airmass.snr(file)
    date = airmass.timeanddate(file)[0]
    expt = airmass.exposuretime(file)
    am = airmass.airmass(file)[0]
    BCCor = airmass.barcor(file)[0]
    alt = airmass.airmass(file)[1]
    # JD = airmass.airmass(file)[2]
    BCCor,HJD = airmass.barcor(file)
    phase = airmass.aphase(HJD)
    filedata.append([file,date,HJD,SNR,expt,fwl,velshift,am,alt,phase,BCCor])
    wl, flux = airmass.reduce_spectrum(file, radial_velocity=18.5)
    filedata.append([wl,flux])
    ld = []
    for line in ll:
        lw, lf, nf = airmass.normalize(wl,flux,line[3],line[4],line[5],line[6],line[3]-5,line[6]+5)
        v,vsini = airmass.wl_to_velocity(lw, line[2])
        # print v
        ld.append([lw,v,lf])
    filedata.append(ld)
    # print ld
    data.append(filedata)

# print legend[2][1][1]
# print data[2][1][1]
MF = [legend,data]
np.save(r'C:\peter\School\Master Scriptie\Data\masterfiles\\eShel.npy',MF)