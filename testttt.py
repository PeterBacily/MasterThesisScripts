from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.style
# import matplotlib as mpl
import glob
import pyfits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
from scipy.optimize import *
from scipy.stats import chi2
from PyAstronomy import pyasl
matplotlib.style.use('classic')
vsini =127
c_light = 299792.458
fl_lowsnr = glob.glob('C:\peter\School\Master Scriptie\Data\eShelData\data\lowsnr/*.fit')
# ['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],
linelist_apo = [ ['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2],  ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
linelist = [4861.01,4471.47,6678.15,4921.93,4647.42,5015.68,5411.52,4650.85,4541.52,5592.25]
ll2= [['H beta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0],['H_gamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0],['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6],['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0],['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914],['He_II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471],['He_II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0],['C_IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
filelist_lapalma = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')
filelist_apo =  glob.glob(r'C:/peter/School/Master Scriptie/Data/eShelData/data/lowsnr*.fit')
nall = [['Na doublet',5892.95,5888.5,5889,5894.2,5894.75]]
ll3 = [[r'H$\beta$', 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I', 5875.621, 5865, 5869.1, 5879.9, 5882]]
# ll_laplma = [[r'Halpha', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'Hbeta', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'Hgamma', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I', 5875.621, 5860.6241285010492, 5861.5995122694139, 5860.6241285010492, 5870.5989810949914], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4537.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4681.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5806.33334373, 5808.314546]]
ll_lapalma = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
ll = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll_lapalma2 = [[r'H$\alpha$', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'H$\beta$', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'H$\gamma$', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4546.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4690.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
# print len(ll_lapalma2)
# print len(ll_lapalma)


for line in ll_lapalma:
    if line[0] == 'He_I':
        line.append('He I '+ str(int(np.round(line[1]))))
    if line[0] == 'He_II':
        line.append('He II '+ str(int(np.round(line[1]))))


def line(x, x0, a1, b1, tau1):
    return a1 * np.exp(-tau1 * np.exp(-((x - x0) / (2 * b1)) ** 2))

def Lapalma_vrad(datafile_folder, plot_save_folder, linelist):
    filelist = glob.glob(datafile_folder+'\*.fits')
    for line in linelist:
        print line[6]
        swl = line[2]-10
        ewl = line[5]+10
        vs,lfs =airmass.overplot_LaPalma(filelist,line,swl,ewl, v_rad=18.5)
        # print vs[0], lfs[0]
        minvs = []
        for i , v in enumerate(vs):
            # print i
            lf = lfs[i]
            index = np.argmin(lf)
            minv = v[index]
            minvs.append(minv)
        # print minvs
        print np.median(minvs)
# Lapalma_vrad('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma',ll_lapalma)
# x = np.arange(-300,300,1)

# y= line(x,0,1,25,0.3)

# plt.scatter(x,y)
# plt.show()


# file = filelist_lapalma[0]
# datafile = pf.open(file)
# header = datafile[0].header
# for item in header:
#     print item
#
# print header['DATE-AVG']
# print header['BJD']

# file = filelist_apo[5]
# # datafile = pf.open(file)
# wl,flux = airmass.extractdata(35,file)
#
# flux1 = flux[wl<4600]
# wl1 = wl[wl<4600]
# plt.xlabel('Wavelength (A)')
# plt.ylabel('Normalized flux')
# plt.plot(wl1,flux1)
# # plt.show()
# plt.savefig(r'C:\peter\School\Master Scriptie\figures\4471_omgeving.pdf')
i=1

def aphase(filetime):
    # datafile = pf.open(file)
    period = 6.829
    jdstart = 2454391.
    erphase = int((filetime-jdstart)/period)*0.001
    unc1 = erphase/period
    unc2 = erphase+0.5/period
    # header =  datafile[0].header
    # filetime = Time(header['MJD-OBS'], format='jd').jd
    phase = ((filetime- jdstart)% period )/ period
    return phase
# print fl_lowsnr
# datfile = pf.open(fl_lowsnr[0])
# header = datfile[0].header
# # print header['DATE-OBS']
# # for item in header:
# #     print item
# print 'nr   phase   date    HJD'
# for file in filelist_apo:
#     b_cor, HJD = airmass.barcor(file)
#     datfile = pf.open(file)
#     header = datfile[0].header
#     dd= header['DATE-OBS']
#     phase = aphase(HJD)
#     wl_rebin, flux_rebin = airmass.reduce_spectrum(file, radial_velocity=18.5)
#     linesused = []
#     # print line[2]
#     snr = airmass.snr(file, 4838.0, 4839.0)
#     print i,'  ',phase,'  ',dd,'  ',HJD
#     i+=1

# print 1/3




def plot_TVS_Lapalma(datafile_folder, plot_save_folder, linelist):
    filelist = glob.glob(datafile_folder+'\*.fits')
    for line in linelist:
        print line[6]
        swl = line[2]-40
        ewl = line[5]+40
        lw,TVS,v,n =airmass.TVS_LaPalma(filelist,line,swl,ewl, v_rad=18.5)
        p = chi2.ppf(0.99, n-1)/(n-1)
        vs,lws = airmass.overplot_LaPalma(filelist,line,v_rad=18.5,startwl=swl,endwl=ewl)
        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec )
        ax1.set_title(line[6])
        ax1.legend()
        # ax1.set_xlim([-600,600])
        # ax1.set_ylim([0.6,1.1])
        ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
        if line[2]==5875.621:
            TVS2 = np.array(TVS)*1.4
            ax2.plot(v, TVS2)
        else:
            ax2.plot(v,TVS)
        ax2.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
        ax2.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
        ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # ax2.set_ylim([0,5])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel('TVS')
        ax2.set_xlim([-600,600])
        plt.savefig(plot_save_folder + r'\\' + line[0] + str(int(np.round(line[1])))+'_TVS.pdf')
        plt.show()
        plt.close()

# plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma',ll_lapalma)

folder = r'D:\Peter\Master Thesis\Data\LaPalmaData'

filelist = glob.glob(folder+'\*.fits')
print filelist
