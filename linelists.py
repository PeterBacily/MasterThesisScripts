from __future__ import division
import matplotlib.pyplot as plt
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
import matplotlib.style
matplotlib.style.use('classic')
vsini =127
c_light = 299792.458

ll = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
# ll_laplma = [[r'H$\alpha$', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'H$\beta$', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'H$\gamma$', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4546.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4690.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
# ll2 = [[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6]]
ll2 = [[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll3 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579],[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll4 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579],[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll_lapalma = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]

# for line in ll:
#     for l2 in ll_lapalma:
#         if line[2]==l2[1]:
#             line[0] = l2[0]
#             line.append(l2[6])
# print ll

fl_clean = glob.glob(r'D:\Peter\Master Thesis\Data\eShelData\data\clean\*.fit')
filelist2 =  glob.glob(r'D:\Peter\Master Thesis\Data\eShelData\data\*.fit')
print len(filelist2)
filelist = glob.glob(r'D:\Peter\Master Thesis\Data\processed\data_vrad18\*')
print len(filelist)
print '--------------'
print filelist
# for file in filelist2:
#     print airmass.snr(file,4991,5001)
# print  filelist2
# file = filelist2[-6]
#
# print airmass.snr(file, 4991, 5001)


def plot_TVS_eShel(filelist, plot_save_folder, linelist):
    # print datafile_folder
    # filelist = glob.glob(datafile_folder+'\*.fit')
    print filelist
    for line in linelist:
        # print line[6]
        swl = line[3]-40
        ewl = line[6]+40
        lw,TVS,v,n =airmass.TVS(filelist,line,swl,ewl, v_rad=18.5)
        p = chi2.ppf(0.99, n-1)/(n-1)
        vs,lws = airmass.overplot(filelist,'-',line, '-',v_rad=18.5,startwl=swl,endwl=ewl)
        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        for i,spec in enumerate(lws):

            ax1.plot(vs[i],spec )

        ax1.set_title(line[6])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        print v.shape, spec.shape
        print v[12]
        print spec[12]
        spec2 = spec[(v > -300)&(v <300)]
        mini = np.floor(10*0.9*np.amin(spec2))/10
        maxi = np.ceil(10*np.amax(spec2))/10
        ax1.set_ylim([mini,maxi])
        ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4
        ax2.plot(v, TVS)
        # else:
        #     ax2.plot(v,TVS)
        ax2.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
        ax2.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
        # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        TVS2 = np.array(TVS)[(v > -300)&(v <300)]
        maxi2 = np.ceil(np.amax(TVS2))
        # print np.amax(spec)
        # print maxi2
        ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel('TVS')
        ax2.set_xlim([-600,600])
        # plt.savefig(plot_save_folder + r'\\LaPalma' + line[0] + str(int(np.round(line[1])))+'_TVS.pdf')
        plt.show()
        plt.close()
# print fl_clean
# print range(len(fl_clean))

line = ['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563']
swl = line[3] - 40
ewl = line[6] + 40
# print fl_clean
vs, lws = airmass.overplot(fl_clean, '-', ['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], '-', v_rad=18.5, startwl=swl, endwl=ewl)
for i, spec in enumerate(lws):
    plt.plot(vs[i], spec)
plt.show()
plt.close()
# plot_TVS_eShel(fl_clean,'savefolder', [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563']])