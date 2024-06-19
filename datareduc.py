from __future__ import division
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
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
import os
import Path_check
import open_masterfiles
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
f = open(filepath_eshel_spectra_info,'r')
dict_eshel = ast.literal_eval(f.read())
f.close()
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
            plt.savefig(plot_save_folder + r'\\LaPalma' + line[0] + str(int(np.round(line[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def plot_TVS_eShel_masterfile(linelist, plot_save_folder,show='off',save='on',sg='on',oneline='on', siglvlline=0.01,datafilefolder=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\apo\test\\'):
    # print datafile_folder
    filelist = open_masterfiles.apo(path=datafilefolder)
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
        ax1.set_title(lineinfo[7])
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
        if isinstance(siglvlline, float):
            Nfiles = len(filelist)
            p = siglvlline
            siglvl = airmass.TVS_significance_level(Nfiles, p)
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
        ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
        ax2.set_xlim([-600,600])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\APO' + lineinfo[0] + str(int(np.round(lineinfo[2])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
            plt.show()
        plt.close()


def plot_TVS_Lapalma_masterfile(linelist, plot_save_folder,show='off',save='on',sg='on',oneline='on', siglvlline=0.01,datafilefolder=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\mercator\test\\'):
    filelist = open_masterfiles.mercator(path=datafilefolder)
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
        ax2.set_ylim([0,maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$'+ r' \ ' + r'$\sigma_{exp}$',size=16)
        ax2.set_xlim([-600,600])
        if save =='on':
            plt.savefig(plot_save_folder + r'\\LaPalma' + lineinfo[0] + str(int(np.round(lineinfo[1])))+'_TVS.pdf',format='pdf', dpi=1200)
        if show =='on':
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
