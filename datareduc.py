from __future__ import division
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
from collections import defaultdict
import numpy as np
import airmass
from scipy.optimize import *
from scipy.stats import chi2
from PyAstronomy import pyasl
import matplotlib.style
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import SavitzkyGolay
warnings.filterwarnings("ignore", category=DeprecationWarning)
import ast
import open_masterfiles
import os
import Path_check


folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)
[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
# print(converted_Data_folder, Data_folder, Plots_folder, Scripts_folder)


matplotlib.style.use('classic')
vsini =127
c_light = 299792.458
# ['He_I', 2 ,6678],['CIII', 4647 +4651]
# ll2 = [['H_alpha', 2, 6562.819, 6554, 6556, 6578, 6579], ['H_beta',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I',30, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I',10, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He_II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He_II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O_III',14, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C_IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
# ll_laplma = [[r'H$\alpha$', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'H$\beta$', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'H$\gamma$', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4546.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4690.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
# ll2 = [[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6]]
ll2 = [[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll3 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579],[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll4 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579],[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll_lapalma = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
ll_TVS_eshel = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
fl_clean = glob.glob(str(Data_folder) + r'\eShelData\data\clean\*.fit')
filelist = glob.glob(str(Data_folder) + r'\LaPalmaData\*.fits')
filelist_lapalma = glob.glob(str(Data_folder)+r'\LaPalmaData\*.fits')
filelist2 = glob.glob(str(Data_folder)+r'/eShelData/data/*.fit')
filepath_eshel_spectra_info = str(Data_folder)+r'\masterfiles\dict_apo_files.txt'
try:
    f = open(filepath_eshel_spectra_info,'r')
    dict_eshel = ast.literal_eval(f.read())
    f.close()
except:
    pass
# del filelist2[0]
# del filelist2[-6:]
# del filelist2[]
# s_fl = filelist
# s_fl.sort(key=lambda file: airmass.snr(file), reverse=True)

# for file in s_fl[:10]:
#     print airmass.snr(file)
# for file in filelist:
#     wl,flux = airmass.reduce_spectrum(file,radial_velocity=-18.5)
#     JD = airmass.airmass(file)[2]
#     a = np.array([wl,flux])
#     np.savetxt(r'C:/peter/School/Master Scriptie/Data/processed/data_vradn18/'+str(JD)+'zeta_ori.txt',a)

# print filelist[12]

# files = [filelist[12]]

# for line in ll:
#     wl, flux= airmass.extractdata(line[1],file)
#     nff1 = flux[(wl>line[3])&(wl<line[4])]
#     nff2 = flux[(wl>line[5])&(wl<line[6])]
#     nff = np.hstack((nff1,nff2))
#     print slope, height
#     normfactor = np.average(nff)
#     x= wl[(wl>line[3])&(wl<line[6])]
#     y= flux[(wl>(line[3]))&(wl<(line[6]))]/normfactor
#     plt.plot(x,y)
#     plt.show()
# for file in filelist2:
#     for line in [ll[1]]:
#         wl, flux = airmass.extractdata(35,file)
#         wl_cor = airmass.wlcorrection(wl)
#         wl_rebin,flux_rebin = airmass.rebin(wl_cor,flux)
#         startwl = line[3]
#         endwl = line[6]
#         plt.plot(wl,flux)
#         plt.plot(wl_rebin,flux_rebin)
#         plt.show()
#         normwl, normflux,slope,height = airmass.normalize(wl,flux,line[3],line[4],line[5],line[6],startwl,endwl)
#         v,vsini = airmass.wl_to_velocity(normwl,line[2],file)
#         print slope, height
#         plt.title(str(line[2]))
#         plt.plot(v,normflux)
#         plt.show()


def plot_TVS(filelist,linelist,vrad=18.5):
    for line in linelist:
        lws,TVS = airmass.TVS(filelist,line,line[3]-20,line[6]+20,v_rad = vrad)
        v,vsini = airmass.wl_to_velocity(lws, line[2])
        plt.title(line[0]+' '+str(round(line[2],1))+' TVS')
        plt.axvline(x=-vsini)
        plt.axvline(x=vsini)
        plt.xlim(-500,500)
        plt.xlabel('Velocity(km/s)',size=15)
        plt.ylabel('$\sigma$'+ '/'+'$\sigma_{exp}}$', size=15)
        # print 'testinguuuuuuuu', len(v),len(TVS)
        plt.plot(v,TVS)
        plt.show()
        plt.close()

#
# def plot_quotient_eShel(file1,file2,line,i):
#     swl = line[3]-40
#     ewl = line[6]+40
#     lw,quo = airmass.quotient(file1,file2,line,swl,ewl)
#     # print lw
#     # print quo
#     v,vsini = airmass.wl_to_velocity(lw, line[2])
#     dt = np.round(airmass.barcor(file2)[1]-airmass.barcor(file1)[1],1)
#     # print airmass.barcor(file2)[1]
#     # print airmass.barcor(file1)[1]
#     plt.title( 'Quotient ' +str(i+1)+'/'+str(i)+' dt = '+ str(dt) + ' d')
#     plt.axvline(x=-vsini)
#     plt.axvline(x=vsini)
#     plt.xlim(-500,500)
#     plt.ylim([0.9,1.15])
#     plt.xlabel('Velocity(km/s)',size=15)
#     plt.ylabel('Quotient Flux', size=15)
#     plt.plot(v,quo, label =line[0]+' '+str(int(np.round(line[2]))) )
#     # plt.show()
#     # plt.close()
# for i in range(len(filelist2)-1):
#     for line in ll3:
#         plot_quotient(filelist2[i], filelist2[i+1], line,i)
#     # plt.plot(lws,qf)
#     plt.legend()
#     plt.savefig(r'C:\peter\School\Master Scriptie\figures\Quotients\eshel\\'+str(i+1)+'_'+str(i)+'_quotient.png')
#     # plt.show()
#     plt.close()
#     print i, '/', len(filelist2)-1
#
# for i in range(len(filelist2)-1):
#     for line in ll4:
#         plot_quotient(filelist2[i], filelist2[i+1], line,i)
#     # plt.plot(lws,qf)
#     plt.legend()
#     plt.savefig(r'C:\peter\School\Master Scriptie\figures\Quotients\eshel\\'+str(i+1)+'_'+str(i)+'3best_lines_quotient.pdf')
#     # plt.show()
#     plt.close()
#     print i, '/', len(filelist2)-1

# pyasl.baryCorr

def plot_faunhover_lines(filelist, startwl = 5886, endwl = 5899):
    for file in filelist:
        wl,flux = airmass.reduce_spectrum(file,radial_velocity=18.5)

        wln, fluxn, _ = airmass.normalize(np.array(wl),np.array(flux),5887.3,5888.3,5891.2,5894.7,startwl,endwl)
        plt.plot(wln,fluxn)
        plt.xlim(startwl,endwl)
        plt.axvline(x=5895.92)
        plt.axvline(x=5889.95)
    plt.show()
    plt.close
    line = ['Fraunhover_Na_doublet',35,5892.95,5887.3,5888.3,5891.2,5894.7]

    lws,TVS = airmass.TVS(filelist,line,line[3]-20,line[6]+20)
    # v,vsini = airmass.wl_to_velocity(lws, line[2], filelist[1])
    plt.title(line[0]+' '+str(round(line[2],1))+' TVS')
    plt.axvline(x=5895.92)
    plt.axvline(x=5889.95)
    plt.xlim(startwl,endwl)
    plt.xlabel('Velocity(km/s)')
    plt.ylabel('$\sigma$'+ '/'+'$\sigma_{exp}}$', size=15)
    # print 'testinguuuuuuuu', len(v),len(TVS)
    plt.plot(lws,TVS)
    plt.show()
    plt.close()

def plot_faunhover_lapalma(filelist, startwl = 5886, endwl = 5899):
    line = ['Na doublet',35,5892.95,5888.5,5889,5894.2,5894.75]
    for file in filelist:
        datafile = pf.open(file)
        header = datafile[0].header
        naxis1 = header['NAXIS1']
        crval1 = header['CRVAL1']
        cdelt1 = header['CDELT1']
        flux = datafile[0].data
        wl = np.exp(np.arange(naxis1)*cdelt1 + crval1 )
        wln, fluxn, _ = airmass.normalize(np.array(wl),np.array(flux),5888.5,5889,5894.2,5894.75,5000,6500)
        plt.plot(wln,fluxn)
        plt.xlim(startwl,endwl)
        # plt.axvline(x=5895.92)
        # plt.axvline(x=5889.95)
    plt.show()
    plt.close
    lws,TVS = airmass.TVS_LaPalma(filelist,line,line[3]-20,line[6]+20)
    v,vsini = airmass.wl_to_velocity(lws, line[2])
    plt.title(line[0]+' '+str(round(line[2],1))+' TVS')
    # plt.axvline(x=5895.92)
    # plt.axvline(x=5889.95)
    plt.xlim(startwl,endwl)
    plt.xlabel('Wavelength $\AA$')
    plt.ylabel('$\sigma$'+ '/'+'$\sigma_{exp}}$', size=15)
    # print 'testinguuuuuuuu', len(v),len(TVS)
    plt.plot(lws,TVS)
    plt.show()
    plt.close()


# plot_faunhover_lapalma(filelist_lapalma)
# plot_quotient(filelist[2],filelist[1],ll[0])

# plot_faunhover_lines(filelist[1:])

def plot_full_spectrum(lapalma_file,apo_file):
    lp_file = pf.open(lapalma_file)
    lp_flux = lp_file[0].data
    header = lp_file[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    lp_wl = np.exp(np.arange(naxis1)*cdelt1 + crval1)
    lp_wl_2, lp_flux_2 = airmass.remove_nan(lp_wl,lp_flux)
    lp_wl_rb, lp_flux_rb = airmass.rebin(lp_wl_2,lp_flux_2)

    # apo_wl, apo_flux = airmass.reduce_spectrum(apo_file)
    apo_wl, apo_flux = airmass.extractdata(35,apo_file)
    f, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(100,50))
    f.suptitle('$\zeta$ Orionis',size=200)
    ax1.plot(apo_wl, apo_flux)
    ax1.set_title('eShel Spectrum',size=100)
    ax2.plot(lp_wl_rb, lp_flux_rb/np.average(lp_flux_rb))
    ax2.set_title('MERCATOR Spectrum',size=100)
    ax1.xaxis.set_ticks(np.arange(min(apo_wl),max(apo_wl+50),100))
    ax1.tick_params(labelsize=60)
    ax2.xaxis.set_ticks(np.arange(min(apo_wl),max(apo_wl+50),100))
    ax2.tick_params(labelsize=60)
    ax2.set_xlabel('Wavelength $\AA$', fontsize=60)
    ax1.set_xlim(min(apo_wl),max(apo_wl))
    ax2.set_xlim(min(apo_wl),max(apo_wl))

    ax1.set_ylim([0,1.5])
    ax2.set_ylim([0,1.5])
    plt.tight_layout()
    f.subplots_adjust(top=0.80)
    # plt.savefig(r'C:\peter\School\Master Scriptie\figures\zet_ori_spec_apo_lp.pdf')
    plt.show()
    plt.close()
# plot_full_spectrum(filelist_lapalma[0],fl_clean[12])
# x = np.linspace(0,1,50)
# y=airmass.sl(x,0.5,np.log(9),0.1,7)
# plt.plot(x,y)
# plt.show()

def double_line(x,x0,a1,b1,tau1,a2,b2,tau2):
    # return a1*np.exp(-tau1*np.exp(-((x-x0)/2*b1)**2))+a2*np.exp(-tau2*np.exp(-((x-x0-14.332)/2*b2)**2))+c
    return a1*np.exp(-tau1*np.exp(-((x-x0)/2*b1)**2))+a2*np.exp(-tau2*np.exp(-((x-x0+5.97)/2*b2)**2))

# x = np.linspace(5870,5905,500)

def line(x,x0,a1,b1,tau1):
    return a1*np.exp(-tau1*np.exp(-((x-x0)/2*b1)**2))


def fitfraun(file):
    apo_wl,apo_flux = airmass.extractdata(35,file)
    dat_x = apo_wl[(apo_wl>5880)&(apo_wl<5905)]
    dat_y = apo_flux[(apo_wl>5880)&(apo_wl<5905)]
    parms, pcov = curve_fit(double_line,dat_x,dat_y,p0=(5896.92,0.55,2.5,0.4,0.55,2.5,0.4))
    # wl_shift.append(parms[0])
    # x = np.linspace(5880,5905,200)
    # y = double_line(x,parms[0],parms[1],parms[2],parms[3],parms[4],parms[5],parms[6])
    # plt.plot(dat_x,dat_y)
    # plt.plot(x,y)
    # plt.show()
    # plt.close()
    return parms[0]

# plt.plot(x,y)
# plt.xlim(5875,5900)
# plt.show()
def rebin_snr(lapalmafile):
    # a=
    lp_file = pf.open(lapalmafile)
    lp_flux = lp_file[0].data
    header = lp_file[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    lp_wl = np.exp(np.arange(naxis1)*cdelt1 + crval1)
    lp_wl_2, lp_flux_2 = airmass.remove_nan(lp_wl,lp_flux)
    lp_wl2_ar = np.array(lp_wl_2)
    lp_flux2_ar = np.array(lp_flux_2)
    # lp_wl_3, lp_flux_3, _ = airmass.normalize(lp_wl2_ar[(lp_wl2_ar>6070)&(lp_wl2_ar<6130)],lp_flux2_ar[(lp_wl2_ar>6070)&(lp_wl2_ar<6130)],6080,6100,6100,6120,6080,6120)
    # print lp_flux_3
    # print np.average(lp_flux_3)/np.std(lp_flux_3)
    # lp_wl_rb, lp_flux_rb = airmass.rebin2(lp_wl_2,lp_flux_2)



def plot_full_spectrum2(lapalma_file,apo_file):
    # lp_file = pf.open(lapalma_file)
    # lp_flux = lp_file[0].data
    # header = lp_file[0].header
    # naxis1 = header['NAXIS1']
    # crval1 = header['CRVAL1']
    # cdelt1 = header['CDELT1']
    # lp_wl = np.exp(np.arange(naxis1)*cdelt1 + crval1)
    # lp_wl_2, lp_flux_2 = airmass.remove_nan(lp_wl,lp_flux)
    # lp_wl_rb, lp_flux_rb = airmass.rebin(lp_wl_2,lp_flux_2)
    # apo_wl, apo_flux = airmass.reduce_spectrum(apo_file)
    lp_wl_rb, lp_flux_rb = airmass.extractdata(35,lapalma_file)
    apo_wl, apo_flux = airmass.extractdata(35,apo_file)
    f, (ax1, ax2) = plt.subplots(2, sharex=True,figsize=(100,50))
    f.suptitle('$\zeta$ Orionis',size=30)
    ax1.plot(apo_wl, apo_flux)
    ax1.set_title('eShel Spectrum',size=20)
    ax2.plot(lp_wl_rb, lp_flux_rb/np.average(lp_flux_rb))
    ax2.set_title('HERMES Spectrum',size=20)
    ax1.xaxis.set_ticks(np.arange(min(apo_wl),max(apo_wl+50),100))
    ax1.yaxis.set_ticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4])
    ax2.xaxis.set_ticks(np.arange(min(apo_wl),max(apo_wl+50),100))
    ax2.yaxis.set_ticks([0,0.2,0.4,0.6,0.8,1,1.2,1.4])
    ax2.set_xlabel('Wavelength $\AA$', fontsize=10)
    ax1.set_xlim([min(apo_wl),max(apo_wl)])
    ax2.set_xlim([min(apo_wl),max(apo_wl)])
    ax1.set_ylim([0,1.5])
    ax2.set_ylim([0,1.5])
    ax2.tick_params(labelsize=10)
    ax1.tick_params(labelsize=10)
    plt.tight_layout()
    f.subplots_adjust(top=0.85,left=0.02,right=0.98)

    # plt.savefig(r'C:\peter\School\Master Scriptie\figures\zet_ori_spec_red2.pdf')
    plt.show()
    plt.close()

# -----------------------------------------------------------------------
# for line in ll:
#     # del filelist2[-6]
#     # del
#     swl = line[3]-40
#     ewl = line[6]+40
#     lw,TVS,v,n =airmass.TVS(filelist2,line,swl,ewl, v_rad=18.5)
#     p = chi2.ppf(0.99, n-1)/(n-1)
#     vs,lws = airmass.overplot(filelist2,filelist_lapalma,line,line,v_rad=18.5,startwl=swl,endwl=ewl)
#     f, (ax1, ax2) = plt.subplots(2, sharex=True)
#
#     for i,spec in enumerate(lws):
#         ax1.plot(vs[i],spec,label=str(i) )
#
#     ax1.set_title(line[0] +' '+ str(int(np.round(line[2]))))
#     ax1.set_xlim([-600,600])
#     ax1.set_ylim([0.6,1.1])
#     # ax1.legend()
#     ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
#     ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
#     if line[2]==5875.621:
#         TVS2 = np.array(TVS)*1.4
#         ax2.plot(v, TVS2)
#     else:
#         ax2.plot(v,TVS)
#     ax2.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
#     ax2.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
#     ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
#     ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
#     ax2.set_ylim([0,5])
#     ax2.set_xlabel('V (km/s)',size='14')
#     ax1.set_ylabel('Normlized Flux',size='14')
#     ax2.set_ylabel('TVS',size='14')
#     ax2.set_xlim([-600,600])
#     plt.savefig(r'C:\peter\School\Master Scriptie\figures\TVS\eshel\\hign_snr_spectra' + str(int(np.round(line[2])))+'.pdf')
#     # plt.show()
#     print line
#     plt.close()

# --------------------------------------------------------
#
# for line in [ll2[0]]:
#     swl = line[3]-40
#     ewl = line[6]+40
#     lw,TVS,v,n =airmass.TVS(filelist2[1:-1],line,swl,ewl, v_rad=18.5)
#     p = chi2.ppf(0.99, n-1)/(n-1)
#     vs,lws = airmass.overplot2(filelist2[1:-2],filelist_lapalma,line,line,v_rad=18.5,startwl=swl,endwl=ewl)
#     f, (ax1, ax2) = plt.subplots(2, sharex=True)
#     for i,spec in enumerate(lws):
#         ax1.plot(vs[i],spec,label=str(i) )
#     ax1.set_title(line[0] +' '+ str(int(np.round(line[2]))))
#     ax1.legend()
#     # ax1.set_xlim([-600,600])
#     # ax1.set_ylim([0.6,1.1])
#     ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
#     ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
#     if line[2]==5875.621:
#         TVS2 = np.array(TVS)*1.4
#         ax2.plot(v, TVS2)
#     else:
#         ax2.plot(v,TVS)
#     ax2.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
#     ax2.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
#     ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
#     ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
#     ax2.set_ylim([0,5])
#     ax2.set_xlabel('V (km/s)')
#     ax1.set_ylabel('Normlized Flux')
#     ax2.set_ylabel('TVS')
#     ax2.set_xlim([-600,600])
#     plt.savefig(r'C:\peter\School\Master Scriptie\figures\TVS\eshel\\norm2' + str(int(np.round(line[2])))+'.png')
#     plt.show()
#     plt.close()
def plot_TVS_eShel(datafile_folder, plot_save_folder, linelist,show='off',save='on',sg='on',oneline='on'):
    # print datafile_folder
    filelist = glob.glob(datafile_folder+'\*.fit')
    # print filelist
    for line in linelist:
        print(line[7])
        swl = line[3]-40
        ewl = line[6]+40
        lw,TVS,v,n =airmass.TVS(filelist,line,swl,ewl, v_rad=18.5)
        sgn = 101  # window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        print(v[sgn] - v[0])
        p = chi2.ppf(0.99, n-1)/(n-1)
        vs,lws = airmass.overplot(filelist,'-',line, '-',v_rad=18.5,startwl=swl,endwl=ewl)
        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec,linewidth=1.0 )
        ax1.set_title(line[7])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec2 = spec[(v>-300)& (v<300)]
        mini = np.floor(10*0.9*np.amin(spec2))/10
        maxi = np.ceil(10*1.01*np.amax(spec2))/10
        ax1.set_ylim([mini,maxi])
        ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4
        ax2.plot(v, TVS, color='b')
        if sg == 'on':
            ax2.plot(v,TVS_smoothed,color='r',linestyle='dashed')
        if oneline == 'on':
            ax2.axhline(y=1, color='gray', linestyle='--')
        # else:
        #     ax2.plot(v,TVS)
        ax2.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax2.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # print len(v)
        # print len(TVS)
        # print v
        TVS2 = np.array(TVS)[(v>-200)& (v<200)]
        # print TVS2
        # print np.amax(TVS2)
        maxi2 = np.ceil(np.amax(TVS2))
        ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
        ax2.set_xlim([-600,600])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\APO' + line[0] + str(int(np.round(line[2])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def plot_TVS_Lapalma(datafile_folder, plot_save_folder, linelist,show='off',save='on',sg='on',oneline='on'):
    filelist = glob.glob(datafile_folder+'\*.fits')
    for line in linelist:
        print(line[6])
        swl = line[2]-40
        ewl = line[5]+40
        lw,TVS,v,n =airmass.TVS_LaPalma(filelist,line,swl,ewl, v_rad=18.5)
        sgn= 151 #window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS,sgn,4)
        # print v[sgn]-v[0]
        p = chi2.ppf(0.99, n-1)/(n-1)
        vs,lws = airmass.overplot_LaPalma(filelist,line,v_rad=18.5,startwl=swl,endwl=ewl)
        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec,linewidth=1.0 )
        ax1.set_title(line[6])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec2 = spec[(v>-300)& (v<300)]
        mini = np.floor(10*0.9*np.amin(spec2))/10
        # maxi = np.ceil(10*1.01*np.amax(spec2))/10
        maxi=np.ceil(np.amax(spec2) / 0.05) * 0.05
        ax1.set_ylim([mini,maxi])
        ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4
        ax2.plot(v, TVS, color='b')
        if sg == 'on':
            ax2.plot(v,TVS_smoothed,color='r',linestyle='dashed')
        if oneline == 'on':
            ax2.axhline(y=1, color='gray', linestyle='--')
        # else:
        #     ax2.plot(v,TVS)
        ax2.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax2.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # print len(v)
        # print len(TVS)
        # print v
        TVS2 = np.array(TVS)[(v>-200)& (v<200)]
        # print TVS2
        # print np.amax(TVS2)
        maxi2 = np.ceil(np.amax(TVS2))
        ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
        ax2.set_xlim([-600,600])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\LaPalma' + line[0] + str(int(np.round(line[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def plot_TVS_eShel_masterfile(linelist, plot_save_folder,show='off',save='on',sg='on',oneline='off', siglvlline=0.01,datafilefolder=None,datareductionprogram='Demetra',norm_boundaries='on',vrange=None):
    # print datafile_folder
    if datareductionprogram == 'AudeLA':
        k=0
        if datafilefolder==None:
            filelist = open_masterfiles.apo()
        else:
            filelist = open_masterfiles.apo(path=datafilefolder)[1:]
    elif datareductionprogram =='Demetra':
        k=0
        if datafilefolder == None:
            filelist = open_masterfiles.apo_demetra()
        else:
            filelist = open_masterfiles.apo_demetra(path=datafilefolder)

    bccor = filelist[0].baricentric_correction
    vrad= -18.5
    # velo_shift = bccor+vrad
    velo_shift=0

    for line in linelist:
        lineinfo = getattr(filelist[0], line).lineinfo
        # print filelist
        # swl = line[3]-40
        # ewl = line[6]+40
        lw,TVS,v,n =airmass.TVS_masterfiles(filelist,line)
        sgn = 101  # window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        # print(v[sgn] - v[0])
        vs,lws = airmass.overplot_masterfiles(filelist,line)

        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec,linewidth=1.0 )
        ax1.set_title(lineinfo[6+k])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec2 = spec[(v>-300)& (v<300)]
        mini = np.floor(100*0.98*np.amin(spec2))/100
        maxi = np.ceil(100*1.02*np.amax(spec2))/100
        ax1.set_ylim([mini,maxi])
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
            ax1.axvline(normv_1+velo_shift, color='k', linestyle='dashed', linewidth=1)
            ax1.axvline(normv_2+velo_shift, color='k', linestyle='dashed', linewidth=1)
            ax1.axvline(normv_3+velo_shift, color='k', linestyle='dashed', linewidth=1)
            ax1.axvline(normv_4+velo_shift, color='k', linestyle='dashed', linewidth=1)

        ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4
        ax2.plot(v, TVS, color='b')
        if sg == 'on':
            ax2.plot(v,TVS_smoothed,color='r',linestyle='dashed')
        if oneline == 'on':
            ax2.axhline(y=1, color='gray', linestyle='--')
        if isinstance(siglvlline, float):
            Nfiles = len(filelist)
            p = siglvlline
            siglvl = airmass.TVS_significance_level(Nfiles, p)
            ax2.axhline(y=siglvl, color='red', linestyle='--')
        # else:
        #     ax2.plot(v,TVS)
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
            ax2.axvline(normv_1+velo_shift, color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(normv_2+velo_shift, color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(normv_3+velo_shift, color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(normv_4+velo_shift, color='k', linestyle='dashed', linewidth=1)
        ax2.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax2.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # print len(v)
        # print len(TVS)
        # print v
        TVS2 = np.array(TVS)[(v>-200)& (v<200)]
        # print TVS2
        # print np.amax(TVS2)
        maxi2 = np.ceil(np.amax(TVS2))
        ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
        if vrange == None:
            ax2.set_xlim([normv_1 - 200, normv_4 + 200])
        else:
            ax2.set_xlim([-vrange, vrange])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\APO_'+datareductionprogram+'_' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def plot_TVS_orders(linelist, plot_save_folder,show='off',save='on',sg='off',oneline='on', siglvlline=0.01,datafilefolder=None,norm_boundaries='on',vrange=None,style=None,from_order=True,es_top=0,es_bottom=0):
    k=0
    if datafilefolder == None:
        print('Datafolder automatically set (this is not the right one, select a datafolder)')
        filelist = open_masterfiles.apo_demetra_orders()
    else:
        filelist = open_masterfiles.apo_demetra_orders(path=datafilefolder)
        print('b')
    print(filelist)
    bccor = filelist[0].baricentric_correction
    vrad= -18.5
    velo_shift = 0
    if style is not None:
        plt.style.use(style)

    for baseline in linelist:
        if from_order is True:
            line=baseline+'_order'
        else:
            line=baseline
        # line=baseline
        lineinfo = getattr(filelist[0], line).lineinfo
        # print filelist
        # swl = line[3]-40
        # ewl = line[6]+40
        lw,TVS,v,n =airmass.TVS_masterfiles_order(filelist,line)
        print(len(v))
        print(len(TVS))
        sgn = 101  # window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        # print(v[sgn] - v[0])
        vs,lws = airmass.overplot_masterfiles_order(filelist,line)
        f,(ax1,ax2) = plt.subplots(2,sharex=True)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec,linewidth=1 )
        ax1.set_title(lineinfo[6+k])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec2 = lws[0][(vs[0]>-1000)& (vs[0]<1000)]
        mini = np.floor(20*np.amin(spec2))/20
        maxi = np.ceil(20*np.amax(spec2))/20
        # extra_space=0.05
        ax1.set_ylim([mini-es_bottom,maxi+es_top])
        # ax1.set_ylim([0.65,1.05])
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
            ax1.axvspan(normv_1+velo_shift, normv_2+velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            ax1.axvspan(normv_3 + velo_shift, normv_4 + velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            # ax1.axvline(normv_1+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax1.axvline(normv_2+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax1.axvline(normv_3+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax1.axvline(normv_4+velo_shift, color='k', linestyle='dashed', linewidth=1)
        ax1.axvline(vsini, color='0.5', linestyle=':', linewidth=1)
        ax1.axvline(-vsini, color='0.5', linestyle=':', linewidth=1)
        # ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        # ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4

        ax2.plot(v, TVS,linewidth=1)
        if sg == 'on':
            ax2.plot(v,TVS_smoothed,color='r',linestyle='dashed')
        if oneline == 'on':
            ax2.axhline(y=1, color='gray', linestyle='--')
        if isinstance(siglvlline, float):
            Nfiles = len(filelist)
            p = siglvlline
            siglvl = airmass.TVS_significance_level(Nfiles, p)
            ax2.axhline(y=siglvl, color='salmon', linestyle='--')
        # else:
        #     ax2.plot(v,TVS)
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
            ax2.axvspan(normv_1+velo_shift, normv_2+velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            ax2.axvspan(normv_3 + velo_shift, normv_4 + velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            # ax2.axvline(normv_1+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax2.axvline(normv_2+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax2.axvline(normv_3+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax2.axvline(normv_4+velo_shift, color='k', linestyle='dashed', linewidth=1)
        ax2.axvline(vsini, color='0.5', linestyle=':', linewidth=1)
        ax2.axvline(-vsini, color='0.5', linestyle=':', linewidth=1)
        # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # print len(v)
        # print len(TVS)
        # print v
        TVS2 = np.array(TVS)[(v>-200)& (v<200)]
        # print TVS2
        # print np.amax(TVS2)
        maxi2 = np.ceil(np.amax(TVS2))
        ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
        if vrange== None:
            ax2.set_xlim([normv_1-200,normv_4+200])
        else:
            ax2.set_xlim([-vrange,vrange])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\APO_orders_' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def plot_TVS_orders_lines_together(linelist, plot_save_folder,show='off',save='on',sg='off',oneline='on', siglvlline=0.01,datafilefolder=None,norm_boundaries='on',vrange=None,style=None,from_order=True,es_top=0,es_bottom=0):
    k=0
    if datafilefolder == None:
        print('Datafolder automatically set (this is not the right one, select a datafolder)')
        filelist = open_masterfiles.apo_demetra_orders()
    else:
        filelist = open_masterfiles.apo_demetra_orders(path=datafilefolder)
        print('b')
    print(filelist)
    bccor = filelist[0].baricentric_correction
    vrad= -18.5
    velo_shift = 0
    if style is not None:
        plt.style.use(style)
    axes = {}
    fig = plt.figure()
    axes['ax0'] = fig.add_subplot(212)
    for i,baseline in enumerate(linelist):
        if from_order is True:
            line=baseline+'_order_rebin'
        else:
            line=baseline
        # line=baseline
        lineinfo = getattr(filelist[0], line).lineinfo
        nlines = len(linelist)
        # subplotnum = (100*nlines)+20+((i+1)*2)
        subplotnum = 200+nlines*10+i+1
        axes[f'ax{i + 1}'] = fig.add_subplot(subplotnum)
        # print filelist
        # swl = line[3]-40
        # ewl = line[6]+40
        lw,TVS,v,n =airmass.TVS_masterfiles_order(filelist,line)
        print(len(v))
        print(len(TVS))
        sgn = 101  # window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        # print(v[sgn] - v[0])
        vs,lws = airmass.overplot_masterfiles_order(filelist,line)
        # f,(ax1,ax2) = plt.subplots(2,sharex=True)
        for j,spec in enumerate(lws):
            axes[f'ax{i + 1}'].plot(vs[j],spec,linewidth=1 )
        axes[f'ax{i + 1}'].set_title(lineinfo[6+k])
        # ax1.legend()
        axes[f'ax{i + 1}'].set_xlim([-600,600])
        spec2 = lws[0][(vs[0]>-1000)& (vs[0]<1000)]
        mini = np.floor(20*np.amin(spec2))/20
        maxi = np.ceil(20*np.amax(spec2))/20
        # extra_space=0.05
        # axes[f'ax{i + 1}'].set_ylim([mini-es_bottom,maxi+es_top])
        axes[f'ax{i + 1}'].set_ylim([0.65,1.1])
        axes[f'ax{i + 1}'].xaxis.get_ticklabels()[1].set_visible(False)
        axes[f'ax{i + 1}'].xaxis.get_ticklabels()[-2].set_visible(False)
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
            axes[f'ax{i + 1}'].axvspan(normv_1+velo_shift, normv_2+velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            axes[f'ax{i + 1}'].axvspan(normv_3 + velo_shift, normv_4 + velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            # ax1.axvline(normv_1+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax1.axvline(normv_2+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax1.axvline(normv_3+velo_shift, color='k', linestyle='dashed', linewidth=1)
            # ax1.axvline(normv_4+velo_shift, color='k', linestyle='dashed', linewidth=1)
        axes[f'ax{i + 1}'].axvline(vsini, color='0.5', linestyle=':', linewidth=1)
        axes[f'ax{i + 1}'].axvline(-vsini, color='0.5', linestyle=':', linewidth=1)
        # ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        # ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4

        axes['ax0'].plot(v, TVS,linewidth=1,label=lineinfo[k+6])

        if sg == 'on':
            axes['ax0'].plot(v,TVS_smoothed,color='r',linestyle='dashed')
    axes['ax1'].set_ylabel('Normlized Flux')
    if oneline == 'on':
        axes['ax0'].axhline(y=1, color='gray', linestyle='--')
    if isinstance(siglvlline, float):
        Nfiles = len(filelist)
        p = siglvlline
        siglvl = airmass.TVS_significance_level(Nfiles, p)
        axes['ax0'].axhline(y=siglvl, color='salmon', linestyle='--')
    # else:
    #     ax2.plot(v,TVS)
    if norm_boundaries == 'on':
        [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
        axes['ax0'].axvspan(normv_1+velo_shift, normv_2+velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
        axes['ax0'].axvspan(normv_3 + velo_shift, normv_4 + velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
        # ax2.axvline(normv_1+velo_shift, color='k', linestyle='dashed', linewidth=1)
        # ax2.axvline(normv_2+velo_shift, color='k', linestyle='dashed', linewidth=1)
        # ax2.axvline(normv_3+velo_shift, color='k', linestyle='dashed', linewidth=1)
        # ax2.axvline(normv_4+velo_shift, color='k', linestyle='dashed', linewidth=1)
    axes['ax0'].axvline(vsini, color='0.5', linestyle=':', linewidth=1)
    axes['ax0'].axvline(-vsini, color='0.5', linestyle=':', linewidth=1)
    # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
    # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
    # print len(v)
    # print len(TVS)
    # print v
    TVS2 = np.array(TVS)[(v>-200)& (v<200)]
    # print TVS2
    # print np.amax(TVS2)
    maxi2 = np.ceil(np.amax(TVS2))
    axes['ax0'].set_ylim([0,6])
    axes['ax0'].set_xlabel('V (km/s)')
    axes['ax0'].set_title("TVS")
    axes['ax0'].set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
    axes['ax0'].legend(prop={'size': 10})
    if vrange== None:
        axes['ax0'].set_xlim([normv_1-200,normv_4+200])
    else:
        axes['ax0'].set_xlim([-vrange,vrange])
    if save =='on':
        plt.savefig(plot_save_folder + r'\\APO_orders_' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_TVS.pdf',format='pdf', dpi=1200)
    if show =='on':
        plt.show()
    plt.close()






def plot_snr_test(filelist,boundaries,rebin=True,rebin_size=0.5):
    for file in filelist:
        wl=file.wl_original
        flux = file.flux_original
        if len(boundaries) == 2:
            [a, d] = boundaries
        elif len(boundaries) == 4:
            [a, b, c, d] = boundaries
        else:
            raise Exception('boundaries  needs to be either list of either 2 or 4 values')
        start = a - 2
        stop = d + 2
        slice = flux[(wl > start) & (wl < stop)]
        wlslice = wl[(wl > start) & (wl < stop)]
        snr=airmass.snr_2(wl, flux, boundaries=boundaries, rebin=True, rebin_size=rebin_size, separate=False)
        if len(boundaries) == 2:
            l = int((len(wlslice) - 1) / 2)
            b = l
            c = l + 1
        if rebin is True:
            wl_rebin = np.arange(wlslice[10], wlslice[-10], rebin_size)
            flux_rebin = airmass.rebin_spec(wlslice, slice, wl_rebin)
            lw, lf, _, _, _ = airmass.normalize(wl_rebin, flux_rebin, a, b, c, d, wl_rebin[0], wl_rebin[-1])
            snr = airmass.snr_2(wl, flux, boundaries=boundaries, rebin=True, rebin_size=rebin_size, separate=False)
        else:
            lw, lf, _, _, _ = airmass.normalize(wlslice, slice, a, b, c, d, a, d)
            snr = airmass.snr_2(wl, flux, boundaries=boundaries, rebin=False, rebin_size=rebin_size, separate=False)
        note=file.mark
        td = file.time_and_date
        print(note, 'snr=',snr)
        plt.plot(lw,lf,label = td+note[-12:])
    plt.legend()
    plt.show()
    plt.close()



def plot_SNR_orders(linelist,filelist,plot_save_folder, file_full_night = None,plot_avg=True,show='on',save='off',norm_boundaries='on',rebin=False,vrange=None,style=None,from_order=True,es_top=0,es_bottom=0,subplotylim = [None,None]):
    k=0
    if file_full_night == None:
        plot_avg=False
    bccor = filelist[0].baricentric_correction
    vrad= -18.5
    date = filelist[0].time_and_date[0:6]
    # velo_shift=bccor+vrad
    velo_shift = 0
    rebin_bin_size = 0.001
    if isinstance(rebin, float):
        rebin_bin_size = rebin
        man_rebin = True
    if isinstance(rebin,bool):
        man_rebin = False
    if style is not None:
        plt.style.use(style)

    for baseline in linelist:
        if from_order is True:
            if rebin is True:
                line = baseline + '_order_rebin'
            else:
                line=baseline+'_order'
        else:
            line=baseline
        lineinfo = getattr(filelist[0], line).lineinfo
        wls,vs,lfs = airmass.overplot_masterfiles_order(filelist,line,manual_rebin=man_rebin,rebin_size=rebin_bin_size,return_wl=True)
        boundary_1,boundary_2,boundary_3,boundary_4 = lineinfo[2 + k], lineinfo[3 + k], lineinfo[4 + k], lineinfo[5 + k]
        boundaries = [boundary_1,boundary_2,boundary_3,boundary_4]

        fig = plt.figure()
        ax1 =fig.add_subplot(211)
        ax2 = fig.add_subplot(223)
        ax3 = fig.add_subplot(224)
        [normv_1, normv_2, normv_3, normv_4], uselessvar = airmass.wl_to_velocity(
            boundaries, lineinfo[1 + k])
        for i,spec in enumerate(lfs):
            snr = airmass.snr_2(wls[i],spec,boundaries=boundaries,rebin=man_rebin,rebin_size=rebin_bin_size,separate=True)
            ax1.plot(vs[i],spec,linewidth=1 ,label = filelist[i].time_and_date)
            ax2.plot(vs[i][(vs[i] > normv_1)&(vs[i]<normv_2)],spec[(vs[i]>normv_1)&(vs[i]<normv_2)])
            ax3.plot(vs[i][(vs[i]>normv_3)&(vs[i]<normv_4)],spec[(vs[i]>normv_3)&(vs[i]<normv_4)])
        if plot_avg is True:
            linedata = getattr(file_full_night, line)
            wl=linedata.wl
            wavenew = np.arange(wl[10], wl[-10], rebin_bin_size)
            linecenter = lineinfo[1]
            v_new = airmass.wl_to_velocity(wavenew, linecenter)[0]
            flux = linedata.flux
            flux_rebin = airmass.rebin_spec(wl, flux, wavenew)
            ax1.plot(v_new, flux_rebin, linewidth=1.5,color='black',label=date+' Average')
            ax2.plot(v_new[(v_new > normv_1) & (v_new < normv_2)], flux_rebin[(v_new > normv_1) & (v_new < normv_2)],linewidth=2,color='black')
            ax3.plot(v_new[(v_new > normv_3) & (v_new < normv_4)], flux_rebin[(v_new > normv_3) & (v_new < normv_4)],linewidth=2,color='black')
        ax1.set_title(lineinfo[6+k]+'  '+ date, fontsize='x-large')
        ax1.set_xlim([-vrange,vrange])
        spec2 = lfs[0][(vs[0]>-1000)& (vs[0]<1000)]
        mini = np.floor(20*np.amin(spec2))/20
        maxi = np.ceil(20*np.amax(spec2))/20
        ax1.set_ylim([mini-es_bottom,maxi+es_top])

        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([boundary_1,boundary_2,boundary_3,boundary_4],lineinfo[1+k])
            ax1.axvspan(normv_1+velo_shift, normv_2+velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            ax1.axvspan(normv_3 + velo_shift, normv_4 + velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
        ax1.axvline(vsini, color='0.5', linestyle=':', linewidth=1)
        ax1.axvline(-vsini, color='0.5', linestyle=':', linewidth=1)
        ax2.set_ylim(subplotylim[0],subplotylim[1])
        ax3.set_ylim(subplotylim[0], subplotylim[1])
        fig.supylabel('Normalized flux',fontsize='large')
        fig.supxlabel('v (km/s)',fontsize = 'large')
        ax2.xaxis.set_major_locator(plt.MaxNLocator(6)) # sets number of ticks on x axis otherwise its too many overlapping ticks
        ax3.xaxis.set_major_locator(plt.MaxNLocator(6)) # sets number of ticks on x axis to match left plot
        box = ax1.get_position() # get values of ax1 box
        ax1.set_position([box.x0, box.y0, box.width * 0.9, box.height]) # use values of ax1 box to shrink it horizontally in order to fit legend box
        ax1.legend(prop={'size':8}, loc='center left', bbox_to_anchor=(1, 0.5),labelspacing = 2)
                  # fancybox=True, shadow=True
        if save =='on':
            plt.savefig(plot_save_folder + r'\\APO_orders_' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_'+date.replace(" ", "_")+'_SNR.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()

def plot_TVS_Lapalma_masterfile(linelist, plot_save_folder,show='off',save='on',sg='on',oneline='on',norm_boundaries = 'on', siglvlline=0.01,datafilefolder=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\mercator\test\\',vrange=None,style=None):
    filelist = open_masterfiles.mercator(path=datafilefolder)
    k=0
    barcor=filelist[0].baricentric_correction
    v_rad=-18.5
    velo_shift=barcor+v_rad
    for line in linelist:
        lineinfo = getattr(filelist[0], line).lineinfo
        # swl = line[2]-40
        # ewl = line[5]+40
        lw,TVS,v,n =airmass.TVS_masterfiles(filelist,line)
        sgn= 151 #window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS,sgn,4)
        # print v[sgn]-v[0]
        vs,lws = airmass.overplot_masterfiles(filelist,line)
        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        for i,spec in enumerate(lws):
            ax1.plot(vs[i],spec,linewidth=1.0 )
        ax1.set_title(lineinfo[6])
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec2 = spec[(v>-300)& (v<300)]

        # maxi = np.ceil(np.amax(spec2) / 0.05) * 0.05
        mini = np.floor(20*np.amin(spec2))/20
        # maxi = np.ceil(20*np.amax(spec2)+1)/20
        # mini = 0.75
        maxi = 1.05
        ax1.set_ylim([mini,maxi])
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
            ax1.axvspan(normv_1+velo_shift, normv_2+velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            ax1.axvspan(normv_3 + velo_shift, normv_4 + velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
        ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4
        ax2.plot(v, TVS, color='b')
        if norm_boundaries == 'on':
            [normv_1,normv_2,normv_3,normv_4],uselessvar = airmass.wl_to_velocity([lineinfo[2+k],lineinfo[3+k],lineinfo[4+k],lineinfo[5+k]],lineinfo[1+k])
            ax2.axvspan(normv_1+velo_shift, normv_2+velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
            ax2.axvspan(normv_3 + velo_shift, normv_4 + velo_shift, facecolor='0.95', edgecolor='0', linestyle='--',alpha=1)
        if sg == 'on':
            ax2.plot(v,TVS_smoothed,color='r',linestyle='dashed')
        if oneline == 'on':
            ax2.axhline(y=1, color='gray', linestyle='--')
        if isinstance(siglvlline, float):
            Nfiles = len(filelist)
            p = siglvlline
            siglvl= airmass.TVS_significance_level(Nfiles,p)
            ax2.axhline(y=siglvl, color='red', linestyle='--')
        # else:
        #     ax2.plot(v,TVS)
        ax2.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax2.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
        # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
        # print len(v)
        # print len(TVS)
        # print v
        TVS2 = np.array(TVS)[(v>-200)& (v<200)]
        # print TVS2
        # print np.amax(TVS2)
        maxi2 = np.ceil(np.amax(TVS2))
        ax2.set_ylim([0,3])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
        if vrange == None:
            ax2.set_xlim([normv_1 - 200, normv_4 + 200])
        else:
            ax2.set_xlim([-vrange, vrange])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\LaPalma' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def plot_TVS_together(linelist=None, filefolder_apo = None,filefolder_merc=None,orders=True,save='off',show='on',plot_save_folder='',oneline='off',sg='off', siglvlline=0.01,):
    linelist_standard = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    # mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
    #                   'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
    apo_master_files = open_masterfiles.apo_demetra_orders(filefolder_apo)
    merc_master_files = open_masterfiles.mercator(filefolder_merc)

    if linelist is None:
        lines =linelist_standard
    else:
        lines=linelist
    vsini = 127

    # for line in [apo_lines[1]]:
    for line in lines:
        f, axarr = plt.subplots(2, 2, sharex=True, figsize = (8,6))
        k=1
        for i, obs in enumerate(['APO', 'MERCATOR']):
            if obs == 'APO':
                master_files = apo_master_files
                if orders is True:
                    line_extension='_order'
                    lw, TVS, v, n = airmass.TVS_masterfiles_order(master_files, line+line_extension)
                else:
                    line_extension=''
                    wl, TVS, v, n = airmass.TVS_masterfiles(master_files, line+line_extension)
            elif obs == 'MERCATOR':
                master_files = merc_master_files
                line_extension = ''
                wl, TVS, v, n = airmass.TVS_masterfiles(master_files, line+line_extension)
            # wl, TVS, v, n = airmass.TVS_masterfiles(master_files, line)
            lineinfo = getattr(master_files[0], line+line_extension).lineinfo

            normwl_edges = [lineinfo[k+1],lineinfo[k+2],lineinfo[k+3],lineinfo[k+4]]
            norm_v_edges =airmass.wl_to_velocity(normwl_edges,lineinfo[k])[0]
            sgn = 91  # window size for SavitzkyGolay (must be odd integer)
            TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
            # print v[sgn]-v[0]
            # p = chi2.ppf(0.99, n - 1) / (n - 1)
            if orders is True:
                vs, lfs = airmass.overplot_masterfiles_order(master_files,line)
            else:
                vs, lfs = airmass.overplot_masterfiles(master_files, line+line_extension)
            # f, (axarr[0, i], axarr[0, i]) = plt.subplots(2, sharex=True)
            spec2 = []
            for j, spec in enumerate(lfs):
                axarr[0, i].plot(vs[j], spec, linewidth=1.0)
                spec2.append(spec[(v > -300) & (v < 300)])

            plt.suptitle(lineinfo[-1], size=24)
            axarr[0, i].set_title(r'$\bf{'+obs+r'}$' +'\n Spectra ')
            axarr[1, i].set_title('TVS ')
            # axarr[0, i].legend()
            # axarr[0, i].set_xlim([-600,600])
            spec3 = np.array(spec2)
            mini = np.floor(20 * 0.99 * np.amin(spec3)) / 20
            maxi = np.ceil(20 * 1.01 * np.amax(spec3)) / 20
            # print maxi, np.amax(spec3)
            axarr[0, i].set_ylim([mini, maxi])
            axarr[0, i].axvline(vsini, color='k', linestyle=':', linewidth=1)
            axarr[0, i].axvline(-vsini, color='k', linestyle=':', linewidth=1)
            for edge in norm_v_edges:
                axarr[0, i].axvline(edge, color='green', linestyle=':', linewidth=1)
            # if line[2]==5875.621:
            #     TVS2 = np.array(TVS)*1.4
            axarr[1, i].plot(v, TVS, color='b')
            if sg == 'on':
                axarr[1, i].plot(v, TVS_smoothed, color='r', linestyle='dashed')
            if oneline == 'on':
                axarr[1, i].axhline(y=1, color='gray', linestyle='--')
            if isinstance(siglvlline, float):
                Nfiles = len(master_files)
                p = siglvlline
                siglvl = airmass.TVS_significance_level(Nfiles, p)
                axarr[1, i].axhline(y=siglvl, color='red', linestyle='--')
            # else:
            #     axarr[1, i].plot(v,TVS)
            axarr[1, i].axvline(vsini, color='k', linestyle=':', linewidth=1)
            axarr[1, i].axvline(-vsini, color='k', linestyle=':', linewidth=1)
            # axarr[1, i].plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
            # axarr[1, i].plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
            # print len(v)
            # print len(TVS)
            # print v
            TVS2 = np.array(TVS)[(v > -200) & (v < 200)]
            # print TVS2
            # print np.amax(TVS2)
            maxi2 = np.ceil(np.amax(TVS2))
            axarr[1, i].set_ylim([0, maxi2])
            # axarr[1,1].set_ylim([0,3])
            axarr[1, i].set_xlabel('V (km/s)',size = 14)
            axarr[0, 0].set_ylabel('Normlized Flux', size = 14)
            axarr[1, 0].set_ylabel(r'$\sigma_{obs}$' + r' \ ' + r'$\sigma_{exp}$', size=20)
            # axarr[1, i].set_xlim([-600, 600])
            # print line[4:]
            # plt.tight_layout()
            plt.subplots_adjust(left=0.08, bottom=None, right=0.98, top=None,
                            wspace=None, hspace=0.15)
        if save == 'on':
            plt.savefig(plot_save_folder + r'\\TVS_' +  lineinfo[0] + line[4:] + '_TVS.pdf', format='pdf',
                            dpi=1200)
        if show == 'on':
            plt.show()
        plt.close()



# plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma',ll_lapalma)
# D:\Peter\Master Thesis\Data\eShelData\data\clean
# plot_TVS_eShel(r'D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\TVS\eShel\every_good_snr',ll_TVS_eshel)
# plot_TVS_eShel(r'D:\Peter\Master Thesis\Data\eShelData\data\clean',r'D:\Peter\Master Thesis\figures\TVS\eShel\cleanonly',ll_TVS_eshel)

# print filelist2[11],filelist[12]
# plot_full_spectrum2(filelist[3],fl_clean[3])
# rebin_snr(filelist_lapalma[0])
# plot_TVS(s_fl[:10],ll)
# for line in ll:
#     airmass.overplot(filelist2[1:-1], filelist_lapalma, line, ll_laplma[6], v_rad = 18.5,together=False)

# print len(filelist2)

# phases = []
# for file in filelist2[1:]:
#     phs = airmass.phase(file)
#     phases.append(phs)
# ews = airmass.ew(ll,filelist2[1:])
# print ews
# plt.scatter(phases,ews)
# plt.ylim(max(ews),min(ews))
# plt.show()



# wl_shift = []
# fwl = []
# for i,file in enumerate(filelist2[1:-1]):
#     try:
#         a = fitfraun(file)
#         wl_shift.append(a-5895.92)
#         fwl.append(a)
#         print i, 'Na lineshift', a -5895.92
#     except RuntimeError:
#         print i, file
#         continue
# disp =  max(wl_shift)-min(wl_shift)
#
# v_disp = c_light*disp/5895.92
#
# v= c_light*np.average(wl_shift)/5895.92
# varray = c_light*np.array(wl_shift)/5895.92
# for i,item in enumerate(varray):
#     print i,item
# print 'v_disp=',v_disp, '     v_avg=',v, 'standard deviation:',np.std(varray)
# plt.plot(varray)
# plt.show()
# plt.plot()

def plot_quotient_eShel(datafile_folder, plot_save_folder, linelist,overplot = 'off',show='off',save='on', sg='off',oneline='on'):
    # print datafile_folder
    filelist = glob.glob(datafile_folder+'\*.fit')
    fileinfo_header = dict_eshel['header']
    # print fileinfo_header
    # print filelist
    for i in range(len(filelist) - 1):
        file1 = filelist[i]
        file2 = filelist[i + 1]
        t1 = airmass.HJD_rounded(file1)
        t2 = airmass.HJD_rounded(file2)
        dt = t2 - t1
        file1_info = dict_eshel[t1]
        file2_info = dict_eshel[t2]
        phi1 = float(file1_info[8])
        phi2 = float(file2_info[8])
        # print phi1,phi2
        dphi=(phi2 - phi1 + 1) % 1
        file1_id = file1_info[1]+str(file1_info[0])
        file2_id = file2_info[1] + str(file2_info[0])
        fi=[file1_info,file2_info]
        for line in linelist:
            print(line[7])
            swl = line[3]-40
            ewl = line[6]+40
            lw,v,qf =airmass.quotient_eShel(file1,file2,line,swl,ewl, v_rad=18.5)
            sgn = 101  # window size for SavitzkyGolay (must be odd integer)
            qf_smoothed = SavitzkyGolay.savitzky_golay(qf, sgn, 4)
            # print v[sgn] - v[0]
            # p = chi2.ppf(0.99, n-1)/(n-1)
# Maak textbox met info en legenda en shit.
            if overplot == 'on':
                vs,lfs = airmass.overplot([file1,file2],'-',line, '-',v_rad=18.5,startwl=swl,endwl=ewl)
                f, (ax1, ax2) = plt.subplots(2, sharex=True)
                for i,spec in enumerate(lfs):
                    lb = fi[i][1]+str(fi[i][0])
                    ax1.plot(vs[i],spec,linewidth=1.0,label = lb)
                # ax1.set_title(line[7])
                # ax1.legend()
                # ax1.set_xlim([-600,600])
                spec2 = spec[(v>-300)& (v<300)]
                mini = np.floor(10*0.9*np.amin(spec2))/10
                maxi = np.ceil(10*1.01*np.amax(spec2))/10
                ax1.set_ylim([mini,maxi])
                ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
                ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
                ax1.set_ylabel('Normalized Flux', size=14)
                ax1.legend(loc=4,fancybox=True,framealpha=0.5)
                # if line[2]==5875.621:
                #     TVS2 = np.array(TVS)*1.4
            else:
                f, (ax2) = plt.subplots(1, sharex=True)
            plt.suptitle('Quotient spectrum ' + line[7] +'\n' + fi[1][1]+str(fi[1][0])+r' / ' + fi[0][1]+str(fi[0][0]), size=14)
            ax2.plot(v, qf, color='b')
            informationtext = '$\Delta$t = '+str(round(dt,3)) +'d'+ '\n'+'$\Delta\phi$ = ' +str(dphi)
            ax2.text(0.815,0.95,informationtext, transform=ax2.transAxes, verticalalignment='top', bbox = dict(boxstyle='round', facecolor='white', alpha=0.5))
            if sg=='on':
                ax2.plot(v,qf_smoothed,color='r',linestyle='dashed')
            if oneline == 'on':
                ax2.axhline(y=1, color='gray', linestyle='--')
            # else:
            #     ax2.plot(v,TVS)
            ax2.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
            ax2.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
            # ax2.plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
            # ax2.plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
            # print len(v)
            # print len(TVS)
            # print v
            TVS2 = np.array(qf)[(v>-200)& (v<200)]
            # print TVS2
            # print np.amax(TVS2)
            # maxi2 = np.ceil(np.amax(TVS2))
            # ax2.set_ylim([0,maxi2])
            ax2.set_xlabel('V (km/s)')
            ax2.set_ylabel(r'Quotient Flux',size=14)
            ax2.set_xlim([-600,600])
            ax2.ticklabel_format(useOffset=False)
            # y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            # ax2.yaxis.set_major_formatter(y_formatter)

            if save =='on' and overplot == 'on':
                # print fi[1][1]+str(fi[1][0])+r'_' + fi[0][1]+str(fi[0][0])+'_'
                plt.savefig(plot_save_folder +r'\\' + fi[1][1]+str(fi[1][0])+r'_' + fi[0][1]+str(fi[0][0])+'_' + line[0] + str(int(np.round(line[2])))+'_Quotient_Overplot.pdf',format='pdf', dpi=1200)
            elif save =='on' and overplot == 'off':
                plt.savefig(plot_save_folder +r'\\' + fi[1][1]+str(fi[1][0])+r'_' + fi[0][1]+str(fi[0][0])+'_' + line[0] + str(int(np.round(line[2]))) + '_Quotient_bare.pdf',format='pdf', dpi=1200)
            if show =='on':
                plt.show()
            plt.close()

def create_JoO_apo_demetra(filefolder,stacked=True):
    files = open_masterfiles.apo_demetra_orders(path = filefolder,manual_filelist=None,sort_data_files='on')
    print(r'\begin{table*}')
    print(r'\begin{tabular}{ l|| c| c| c| c|c|c|c|c|c|c|c|c }')
    print(r'\# & Date & HJD & T$_{\textrm{exp}}$  & SNR & SNR rate& efficiency&Airmass & Alt & Phase & BC& v$_{\textrm{ISM}}$\\')
    print(r'APO & 2016 & $-$2457000 & (s)&H$\alpha$ & SNR 60s$^{-1}$& SNR$^{2}$ 60s$^{-1}$& &(deg)& (6.83 d)& km\,s$^{-1}$) & km\,s$^{-1}$))\\')
    print(r'\hline')
    # files.sort(key=lambda x: x.i)
    i=1
    for file in files:
        if file.mark == 0:
            note = ''
        else:
            note = str(file.mark)
        # nf_ha=file.line6562.normalizationflux
        # snr_ha = 1/np.std(nf_ha)
        a,b,c,d = 6549, 6550.7, 6576.0, 6578.0
        wave = file.line6562_order_rebin.wl
        flux = file.line6562_order_rebin.flux
        normwave = np.hstack((wave[(wave > a) & (wave < b)], wave[(wave > c) & (wave < d)]))
        normflux = np.hstack((flux[(wave > a) & (wave < b)], flux[(wave > c) & (wave < d)]))
        # fit line trough slice
        # print('ss', startwl, endwl)
        slope, height = np.polyfit(normwave, normflux, 1)
        # print 'slope and height are', slope, height
        fit = np.poly1d([slope, height])
        nnf = []
        for k, nwl in enumerate(normwave):
            nnf.append(normflux[k] / fit(nwl))
        snr_ha = 1 / np.std(nnf)
        snr_per_60 = snr_ha*60/file.exptime
        count_per_60 = (snr_ha**2)*60/file.exptime
        if stacked is True:
            obs_prefix='A'
        elif stacked is False:
            obs_prefix='APO'
        print(obs_prefix,i, '&', file.time_and_date, '&', "{:.3f}".format(file.HJD - 2457000), '&', "{:.0f}".format(
            file.exptime), '&', "{:.0f}".format(snr_ha), '&',"{:.0f}".format(snr_per_60), '&',"{:.0f}".format(count_per_60), '&', "{:.1f}".format(file.airmass), '&', "{:.0f}".format(
            file.alt), '&', "{:.3f}".format(file.phase), '&', "{:.0f}".format(
            file.baricentric_correction), '&', "{:.0f}".format(file.velshift),r'\\')
        i+=1
    print(r'\end{tabular}')
    print(r'\end{table*}')

def create_JoO_apo_audela(filefolder):
    files =  open_masterfiles.apo(wantedmarks = None,path = filefolder,manual_filelist=None)
    print(r'\begin{tabular}{ l|| r| r| r| r|r|r|r|r|r|l }')
    print(r'\# & Date & HJD & T$_{\textrm{exp}}$  & SNR & Airmass & Alt & Phase & BC& v$_{\textrm{ISM}}$ & Notes\\')
    print(r' APO & 2016 & $-$2457000 & (s)& & & (deg)& p = 6.83 d& (km/s) & (km/s) &          \\')
    print(r'\hline')
    # files.sort(key=lambda x: x.i)
    for file in files:
        if file.mark == 0:
            note = ''
        else:
            note = str(file.mark)
        print('APO', file.i, '&', file.time_and_date, '&', "{:.3f}".format(file.HJD - 2457000), '&', "{:.0f}".format(
            file.exptime), '&', "{:.0f}".format(file.snr), '&', "{:.1f}".format(file.airmass), '&', "{:.0f}".format(
            file.alt), '&', "{:.3f}".format(file.phase), '&', "{:.0f}".format(
            file.baricentric_correction), '&', "{:.0f}".format(file.velshift), '&', note, r'\\')


def create_JoO_mercator(filefolder):
    files =  open_masterfiles.mercator(path = filefolder,manual_filelist=None)
    print(r'\begin{table*}')
    print(r'\begin{tabular}{ l|| c| c| c| c|c|c|c|c|c|c|c }')
    print(r'\# & Date & HJD & T$_{\textrm{exp}}$  & SNR & SNR rate& efficiency&Airmass & Alt & Phase & BC& v$_{\textrm{ISM}}$ \\')
    print(r'Mercator & 2015 & $-$2457000 & (s)&H$\alpha$ & SNR 60s$^{-1}$& SNR$^{2}$ 60s$^{-1}$& (deg)& (6.83 d)& km\,s$^{-1}$) & km\,s$^{-1}$)) \\')
    print(r'\hline')
    # files.sort(key=lambda x: x.i)
    i = 1
    note=''
    for file in files:
        a, b, c, d = 6549.7, 6550.7, 6577.0, 6578.0
        wave = file.line6562_rebin.wl
        flux = file.line6562_rebin.flux
        normwave = np.hstack((wave[(wave > a) & (wave < b)], wave[(wave > c) & (wave < d)]))
        normflux = np.hstack((flux[(wave > a) & (wave < b)], flux[(wave > c) & (wave < d)]))
        # fit line trough slice
        # print('ss', startwl, endwl)
        slope, height = np.polyfit(normwave, normflux, 1)
        # print 'slope and height are', slope, height
        fit = np.poly1d([slope, height])
        nnf = []
        for k, nwl in enumerate(normwave):
            nnf.append(normflux[k] / fit(nwl))
        snr_ha = 1 / np.std(nnf)
        snr_per_60 = snr_ha*60/file.exptime
        count_per_60 = (snr_ha**2)*60/file.exptime
        print('MERC', file.i, '&', file.time_and_date, '&', "{:.3f}".format(file.HJD - 2457000), '&', "{:.0f}".format(
            file.exptime), '&', "{:.0f}".format(snr_ha), '&',"{:.0f}".format(snr_per_60), '&',"{:.0f}".format(count_per_60), '&', "{:.1f}".format(file.airmass), '&', "{:.0f}".format(
            file.altitude), '&', "{:.3f}".format(file.phase), '&', "{:.0f}".format(
            file.baricentric_correction), '&', "{:.0f}".format(file.velshift),  r'\\')
        i+=1
    print(r'\end{tabular}')
    print(r'\end{table*}')

def plot_order_stack(data_individual_list,wlpiece= [5335, 5345],rebinstep=0.1,day = 'All',from_order=True):
    groups = defaultdict(list)

    for obj in data_individual_list:
        # print(obj.time_and_date)
        groups[obj.time_and_date[0:5]].append(obj)

    new_list = groups.values()
    list_of_day_data = list(new_list)
    if day == 'All':
        list_of_day_data_slice = list_of_day_data
    else:
        list_of_day_data_slice = [list_of_day_data[day]]
    for testday in list_of_day_data_slice:
        wlarray = airmass.find_order(wlpiece, testday[0]).wl_original[5:-5]
        wl_rebin = np.arange(wlarray[0], wlarray[-1], rebinstep)
        day_data = []
        for observation in testday:
            relevant_order = airmass.find_order(wlpiece, observation)
            if from_order is True:
                wl1 = relevant_order.wl_original
                flux1 = relevant_order.flux_original
                addtext = 'taken from individual orders'
            elif from_order is False:
                wl1 = observation.wl_original
                flux1 = observation.flux_original
                addtext = 'taken from order merged spectra'
            flux_rebin = airmass.rebin_spec(wl1, flux1, wl_rebin)
            snr2 = airmass.snr_2(wl_rebin, flux_rebin, boundaries=wlpiece)
            print(snr2)
            day_data.append([wl1, flux1, wl_rebin, flux_rebin, snr2,observation])

        spec_list = []
        weightlist = []
        i = 0
        for spec in day_data:
            spc_avg = np.average(spec[3][10:-10])
            normspec = spec[3] / spc_avg
            spec_list.append(normspec)
            weightlist.append(spec[4] ** 2)
            plt.plot(spec[2], normspec + 0.1 * i,label=spec[5].time_and_date)
            i += 0
        avg_spec = np.average(spec_list, axis=0, weights=weightlist)
        # print('snr_combined  measured = ', airmass.snr_2(spec[2], avg_spec, boundaries=wlpiece))
        # print('snr_combined theoretical = ', np.sqrt(np.sum(weightlist)))

        for value in weightlist:
            print(np.sqrt(value))
        plt.plot(spec[2], avg_spec, c='black',label = 'Day Average',linewidth=1.5)
        plt.xlim(wlpiece[0],wlpiece[1])
        plt.xlabel('Wavelength (Å)')
        plt.ylabel('Relative Flux')
        plt.title('Flat piece of continuum in APO data '+addtext)
        plt.ylim(0.96,1.05)
        plt.ticklabel_format(useOffset=False)
        plt.legend()
        plt.show()
        plt.close()
