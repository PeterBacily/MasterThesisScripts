from __future__ import division
import matplotlib.pyplot as plt
import glob
import warnings
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
from scipy import interpolate
from scipy.optimize import *
from PyAstronomy import pyasl
import scipy.stats as ss
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import SavitzkyGolay
warnings.simplefilter('ignore')
from SavitzkyGolay import savitzky_golay
from pysynphot import observation
from pysynphot import spectrum
warnings.resetwarnings()
import airmass
import open_masterfiles



def overplot_masterfiles_with_wl(filelist,line,separate_lines=False):
    vsini = 127
    lfs = []
    vs = []
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
        v = linedata.v_cor
        wls.append(wl)
        lfs.append(np.array(flux) + a)
        vs.append(v)
        a+=ap
    return wls,vs, lfs


def wl_to_velocity(wavelengths, linecenter):
    v = []
    for l in wavelengths:
        velo = (299792.5*(l/linecenter - 1))
        v.append(velo)
    return np.array(v)

def velocity_to_wl(velocities,linecenter):
    wls = []
    for v in velocities:
        wl=linecenter*(1+(v/299792.5))
        wls.append(wl)
    return np.array(wls)


def plot_TVS_eShel_masterfile(linelist, plot_save_folder,show='off',save='on',datafilefolder=None,datareductionprogram='Demetra',norm_boundaries='on'):
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
    for line in linelist:
        lineinfo = getattr(filelist[0], line).lineinfo
        # print filelist
        # swl = line[3]-40
        # ewl = line[6]+40
        lw,TVS,v,n =airmass.TVS_masterfiles(filelist,line)
        sgn = 101  # window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        # print(v[sgn] - v[0])
        wls,vs,lws = overplot_masterfiles_with_wl(filelist,line)
        f, (ax1, ax2) = plt.subplots(2)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec,linewidth=1.0 )
            ax2.plot(wls[i],spec,linewidth=1.0)
        ax1.set_title(lineinfo[6+k])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec2 = spec[(v>-300)& (v<300)]
        mini = np.floor(10*0.9*np.amin(spec2))/10
        maxi = np.ceil(10*1.01*np.amax(spec2))/10
        ax1.set_ylim([mini,maxi])
        ax2.set_ylim([mini, maxi])
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
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
        ax2.set_xlim([normv_1-200,normv_4+200])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\APO_'+datareductionprogram+'_' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()