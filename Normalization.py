from __future__ import division
import time
import matplotlib.pyplot as plt
import glob
import warnings
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
from scipy.stats import *
from scipy.optimize import *
import open_masterfiles
import Datafile_class

def overplot_masterfiles_with_wl(filelist,line,separate_lines=False):
    vsini = 127
    lfs = []
    v_cors = []
    vs=[]
    wls=[]
    a = 0.0
    if separate_lines:
        ap = 0.1
    elif not separate_lines:
        ap = 0.0
    else:
        raise SyntaxError('separate_lines needs to be Boolean')
    for file in filelist:
        linedata = getattr(file, line)
        flux = linedata.flux
        wl = linedata.wl
        v_cor = linedata.v_cor
        v= linedata.v
        wls.append(wl)
        lfs.append(np.array(flux) + a)
        v_cors.append(v_cor)
        vs.append(v)
        a+=ap
    return wls,vs, v_cors,lfs


def wl_to_velocity(wavelengths, linecenter):
    v = []
    for l in wavelengths:
        velo = (299792.5*(l/linecenter - 1))
        v.append(velo)
    return np.array(v)

def velocity_to_wl(velocities,linecenter,voff):
    wls = []
    for v in velocities:
        wl=linecenter*(1+((v+voff)/299792.5))
        wls.append(wl)
    return np.array(wls)


def plot_TVS_eShel_masterfile(linelist, plot_save_folder, custom_class_object_list=None,show='on',save='off',datafilefolder=None,datareductionprogram='Demetra',norm_boundaries='on'):
    # print datafile_folder
    if datareductionprogram == 'AudeLA':
        k=1
        if datafilefolder==None:
            filelist = open_masterfiles.apo()
        else:
            filelist = open_masterfiles.apo(path=datafilefolder)
    elif datareductionprogram =='Demetra':
        k=0
        if datafilefolder == None:
            filelist = open_masterfiles.apo_demetra()
        else:
            filelist = open_masterfiles.apo_demetra(path=datafilefolder)
    if custom_class_object_list == None:
        pass
    elif isinstance(custom_class_object_list, list):
        filelist=custom_class_object_list
    vrad=18.5
    v_offset = 0
    print(v_offset)
    for line in linelist:
        lineinfo = getattr(filelist[0], line).lineinfo
        # print filelist
        # swl = line[3]-40
        # ewl = line[6]+40
        # lw,TVS,v,n =airmass.TVS_masterfiles(filelist,line)
        sgn = 101  # window size for SavitzkyGolay (must be odd integer)
        # TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        # print(v[sgn] - v[0])
        wls,vs, v_cors,lfs = overplot_masterfiles_with_wl(filelist,line)
        f, (ax1, ax2) = plt.subplots(2)
        for i,spec in enumerate(lfs):
            ax1.plot(vs[i],spec,linewidth=1.0 )
            ax2.plot(wls[i],spec,linewidth=1.0)
        ax1.set_title(lineinfo[6+k])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec2 = spec[(vs[0]>-300)& (vs[0]<300)]
        mini = np.floor(10*0.9*np.amin(spec2))/10
        maxi = np.ceil(10*1.1*np.amax(spec2))/10
        ax1.set_ylim([mini,maxi])
        ax2.set_ylim([mini, maxi])
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4] = airmass.wl_to_velocity_2([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k],v_offset)
            ax1.axvline(normv_1, color='k', linestyle='dashed', linewidth=1)
            ax1.axvline(normv_2, color='k', linestyle='dashed', linewidth=1)
            ax1.axvline(normv_3, color='k', linestyle='dashed', linewidth=1)
            ax1.axvline(normv_4, color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(lineinfo[2+k], color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(lineinfo[3+k], color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(lineinfo[4+k], color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(lineinfo[5+k], color='k', linestyle='dashed', linewidth=1)

        # ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        # ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4
        # else:
        #     ax2.plot(v,TVS)
        # if norm_boundaries == 'on':
        #     [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
        #     ax2.axvline(normv_1, color='k', linestyle='dashed', linewidth=1)
        #     ax2.axvline(normv_2, color='k', linestyle='dashed', linewidth=1)
        #     ax2.axvline(normv_3, color='k', linestyle='dashed', linewidth=1)
        #     ax2.axvline(normv_4, color='k', linestyle='dashed', linewidth=1)
        # ax2.axvline(vsini, color='k', linestyle=':', linewidth=1)
        # ax2.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # print len(v)
        # print len(TVS)
        # print v
        # TVS2 = np.array(TVS)[(v>-200)& (v<200)]
        # print TVS2
        # print np.amax(TVS2)
        # maxi2 = np.ceil(np.amax(TVS2))
        # ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel('Normlized Flux',size=16)
        vlim1,vlim2 =normv_1-200,normv_4+200
        [wl_lim1,wl_lim2] = velocity_to_wl([vlim1,vlim2],lineinfo[1],v_offset)
        ax2.set_xlim([wl_lim1,wl_lim2])
        ax1.set_xlim([vlim1, vlim2])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\APO_'+datareductionprogram+'_' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def generate_custom_class_object_list(datafilefolder,linelist):
    filelist = glob.glob(datafilefolder+r'*.fit')
    class_object_list=[]
    for file in filelist:
        a = Datafile_class.Datafile_apo_demetra(file,ll_file=linelist,mark='Normalization_object_custom')
        class_object_list.append(a)
    return class_object_list
def line(x,x0,a1,b1,tau1):
    return a1*np.exp(-tau1*np.exp(-((x-x0)/2*b1)**2))
def line2(x,x0,a1,b1,tau1):
    return a1*np.exp(-tau1*np.exp(-((x-x0)/2*b1)**2))

def line3(x, x0,y0, a, sigma):
    return y0 + a * np.exp(-(x - x0)**2 / (2 * sigma**2))


def wl_shift(objectlist,line):
    vshiftlist =[]
    for object in objectlist:
        linedata = getattr(object, line)
        lineinfo = linedata.lineinfo
        lc =lineinfo[1]
        flux = linedata.flux
        wl = linedata.wl
        v = linedata.v_cor

        parms, pcov = curve_fit(line3,v,flux,p0 = (0, 0.9, 0.1, 10))
        plt.plot(v,flux)
        plt.plot(v,line3(v,*parms))
        vshift = parms[0]
        vshiftlist.append(vshift)
    avg_vshift=np.average(vshiftlist)
    std_vshift = np.std(vshiftlist)
    print(line,avg_vshift,std_vshift)
    plt.show()
    plt.close()
        # return parms,pcov




demetra_file_dir = r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\Zet_Ori_Data_Zet_Ori_Response\final_spectra\good\\'
linelist_apo = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
     ['Hb', 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'],
     ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'],
     ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5875'],
     ['He_II', 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4541'],
     ['He_II', 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4685'],
     ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5411'],
     ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'],
     ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4921'],
     ['He_I', 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'],
     ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'],
     ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]

linelist_apo_new = [['Ha', 6562.819, 6550.2, 6551.2, 6578, 6579, r'H$\alpha$ 6563'],
     ['Hb', 4861.333, 4848.0, 4849.0, 4877.0, 4878.0, r'H$\beta$ 4861'],
     ['He_I', 4713.1457, 4708.3, 4709.3, 4718, 4719, 'He I 4713'],
     ['He_I', 5875.621, 5867, 5868, 5881.7, 5882.7, 'He I 5875'],
     ['He_II', 4541.6, 4537.5, 4538.5, 4546.6, 4547.6, 'He II 4541'],
     ['He_II', 4685.804, 4682.3, 4683.3, 4690.5, 4691.5, 'He II 4685'],
     ['He_II', 5411.521, 5405.2, 5406.2, 5417.0, 5418.2, 'He II 5411'],
     ['He_I', 4471.4802, 4463.0, 4464, 4476.2, 4477.2, 'He I 4471'],
     ['He_I', 4921.93, 4914.2, 4915.2, 4928.2, 4929.2, 'He I 4921'],
     ['He_I', 6678.15, 6670, 6671, 6685, 6686, 'He I 6678'],
     ['O_III', 5592.37, 5588.3, 5589.3, 5597.5, 5598.5, 'O III 5592'],
     ['C_IV', 5801.33, 5793.8, 5794.8, 5817.1, 5818.1, 'C IV 5801']]


linelist_apo2 = [['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861']]
apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
             'line4861', 'line4921', 'line6678', 'line4471']
apo_lines2 = ['line6562']
objectlist = generate_custom_class_object_list(demetra_file_dir,r'D:\peter\Master_Thesis\Datareduction\Converted_Data\linelists\linelist_apo.txt')
# for line in apo_lines:
#     wl_shift(objectlist,line)

plot_TVS_eShel_masterfile(apo_lines,r'D:\peter\Master_Thesis\Datareduction\Plots\normalization\test\apo\demetra',custom_class_object_list=objectlist)