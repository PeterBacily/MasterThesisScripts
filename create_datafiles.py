from __future__ import division
import glob
import astropy.io.fits as pf
import numpy as np
import airmass
import matplotlib.style
import pickle
# import Datafile_class
import os
from pathlib import Path
import re
import Path_check
import pathlib
import sys
from astropy.timeseries import LombScargle
import open_masterfiles
from Datafile_class import *
import tqdm
# sys.modules['Line'] = Line
# sys.modules['Datafile_mercator'] = Datafile_mercator
# sys.modules['Datafile_apo'] = Datafile_apo
# import matplotlib.pyplot as plt
# from astropy.time import Time
# import math
# import calendar
# import numpy as np
# from scipy.optimize import *
# from scipy.stats import chi2
# from PyAstronomy import pyasl
# import datareduc

folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)

[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)




matplotlib.style.use('classic')
vsini = 127
c_light = 299792.458
# ['He_I', 2 ,6678],['CIII', 4647 +4651]
# ll2 = [['H_alpha', 2, 6562.819, 6554, 6556, 6578, 6579], ['H_beta',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He_I',30, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He_I',10, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He_II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He_II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He_II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O_III',14, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C_IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$',35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
# ll_laplma = [[r'H$\alpha$', 6562.819, 6537.819, 6541.819, 6576.819, 6580.819], [r'H$\beta$', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0], [r'H$\gamma$', 4340.472, 4333.0, 4334.0, 4352.0, 4353.0], ['He I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914], ['He I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0], ['He I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6], ['He II', 4541.6, 4537.6154685, 4538.58308527, 4546.6154685, 4547.58471471], ['He II', 4685.804, 4681.82450988, 4682.7985253, 4690.82450988, 4691.79337991], ['He II', 5411.521, 5386.5435547289835, 5387.4960723716758, 5386.5435547289835, 5436.4963825520754], ['O III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV', 5801.33, 5791.33334373, 5792.32731839, 5791.33334373, 5808.314546]]
# ll2 = [[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6]]
ll2 = [[r'H$\beta$', 35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0], ['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll3 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579], [r'H$\beta$', 35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882], ['He I',35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15], ['He II',35, 4541.6, 4523, 4529, 4546, 4548.5], ['He II',35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3], ['He II',35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2], ['O III',35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0], ['C IV',35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5]]
ll4 = [[r'H$\alpha$', 35, 6562.819, 6554, 6556, 6578, 6579],[r'H$\beta$',35, 4861.333, 4838.0-5.0, 4839.0, 4880.0, 4881.0+5.0],['He I',35, 5875.621, 5865, 5869.1, 5879.9, 5882]]
ll_lapalma = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
ll_TVS_eshel = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
ll_nieuwe_lijnen = [['He_I', 35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, r'He I 6678']]

# ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, r'He I 6678']
fl_eshel_clean_folder = str(Data_folder)+r'\eShelData\data\clean' #Zonder twee spectra met rare Halpha spike
filelist_lapalma_folder = str(Data_folder)+r'\LaPalmaData'
fl_eshel_demetra_folder = str(Data_folder)+r'\Demetra\spectra'
fl_eshel_all_folder = str(Data_folder)+r'\eShelData\data\AlleSpectra'
fl_eshel_goodSNR_folder = str(Data_folder)+r'\eShelData\data'
fl_eshel_demetra = glob.glob(fl_eshel_demetra_folder+r'\*.fit')
fl_eshel_all = glob.glob(fl_eshel_all_folder+r'\*.fit')
fl_clean = glob.glob(fl_eshel_clean_folder+r'\*.fit')
fl_goodSNR = glob.glob(fl_eshel_goodSNR_folder+r'\*.fit')
filelist_lapalma = glob.glob(filelist_lapalma_folder+r'\*.fits')

datafile_folder_merc = str(converted_Data_folder)+r'\mercator\\'
datafile_folder_omar = str(converted_Data_folder)+r'\dataset_omar\\'
datafile_folder_apo = str(converted_Data_folder)+r'\apo\\'
datafile_folder_demetra = str(converted_Data_folder)+r'\demetra\\'
datafile_folder_demetra_test= str(converted_Data_folder)+r'\test\demetra\\'
test_datafile_folder = str(converted_Data_folder)+r'\test\\'
datafile_folder_audela_all = str(converted_Data_folder)+r'\AudeLA\all\\'

fl_demetra_all = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\Zet_Ori_Data_Zet_Ori_Response\final_spectra\*.fit')
fl_demetra_good = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\Zet_Ori_Data_Zet_Ori_Response\final_spectra\good\*.fit')

fl_demetra_all_alt = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\Zet_Ori_Data_Altair_Response\final_spectra\*.fit')
fl_demetra_good_alt = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\Demetra\Zet_Ori_Data_Altair_Response\final_spectra\good\*.fit')
fl_apo_audela_all = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\eShelData\data\AlleSpectra\*.fit')
fl_dataset_omar = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\Dataset_Omar\fits\*.fits')



def bjd(file):
    fits = pf.open(file)
    header = fits[0].header
    if 'BJD' in header:
        HJD =  float(header['BJD'])
        fits.close()
        return HJD
    elif 'HJD' in header:
        HJD = float(header['HJD'])
        fits.close()
        return HJD
    else:
        print(file)
        fits.close()


# testfile = pf.open(fl_dataset_omar[0])
# header = testfile[0].header
# for item in header:
#     print(item, header[item])
# sortedfl_omar = sorted(fl_dataset_omar, key=lambda x: bjd(x), reverse=False)
sortedfl_lapalma = sorted(filelist_lapalma, key=lambda x: bjd(x), reverse=False)
k = 1
def create_datafiles_demetra(filelist=fl_eshel_demetra,savefolder=datafile_folder_demetra_test,linelist_file_path=None):
    i=0
    for file in filelist:
        a = Datafile_apo_demetra(file,ll_file=linelist_file_path)
        dl, dl2 = airmass.split_date(a.header['DATE-OBS'])
        savename = savefolder+a.observatory+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(a, workfileresource)
        workfileresource.close()
        i+=1
def create_datafiles_audela(filelist=fl_apo_audela_all,savefolder=datafile_folder_audela_all,linelist_file_path=None,vshift=True):
    i=0
    for file in filelist:
        a = Datafile_apo(file,ll_file=linelist_file_path,velo_shift=vshift)
        dl, dl2 = airmass.split_date(a.header['DATE-OBS'])
        savename = savefolder+a.observatory+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(a, workfileresource)
        workfileresource.close()
        i+=1
def time_from_filename(filename):
    a = Path(filename)
    b=str(a.stem)[-15:]
    c=b[:8]+b[9:]
    return(int(c))
def create_datafiles_demetra_orders(datafolder,savefolder,linelist_file,snrtreshhold=None,vshift=True,filename_prefix='',note=0,rebin_size=0.5):
    zipfiles = sorted(glob.glob(datafolder+r'*.zip'), key=lambda x: time_from_filename(x))
    print(zipfiles)
    full_spec_files = sorted(glob.glob(datafolder+r'*.fit'), key=lambda x: time_from_filename(x))
    print(full_spec_files)
    for i in range(len(zipfiles)):
        if str(Path(zipfiles[i]).stem)[-15:] ==str(Path(full_spec_files[i]).stem)[-15:]:
            print(str(i))
            orders_class_object=Datafile_apo_demetra_with_orders(zipfiles[i],full_spec_files[i],ll_file=linelist_file,velo_shift=vshift,mark=note,rebin_size=rebin_size)
            dl, dl2 = airmass.split_date(orders_class_object.header['DATE-OBS'])
            snr=orders_class_object.snr_original
            savename = savefolder + filename_prefix + orders_class_object.observatory + '_' + dl[0] + dl[1] + dl[2] + dl[3] + '.txt'
            savefolder_snr=savefolder+r'snr_'+str(snrtreshhold)+r'\\'
            savename_snr = savefolder_snr + filename_prefix + orders_class_object.observatory + '_' + dl[0] + dl[1] + dl[2] + dl[3] + '.txt'
            workfileresource = open(savename, 'wb')
            pickle.dump(orders_class_object, workfileresource)
            workfileresource.close()
            if snrtreshhold is not None:
                Path(savefolder_snr).mkdir(parents=True, exist_ok=True)
                if snr>snrtreshhold:
                    workfileresource = open(savename_snr, 'wb')
                    pickle.dump(orders_class_object, workfileresource)
                    workfileresource.close()
        else:
            print('YOU SUCK')

def create_test_version_datafiles_demetra_orders(datafolder,linelist_file,snrtreshhold=None):
    zipfiles = sorted(glob.glob(datafolder+r'*.zip'), key=lambda x: time_from_filename(x))
    print(zipfiles)
    full_spec_files = sorted(glob.glob(datafolder+r'*.fit'), key=lambda x: time_from_filename(x))
    print(full_spec_files)
    testfile_list=[]
    for i in range(len(zipfiles)):
        if str(Path(zipfiles[i]).stem)[-15:] ==str(Path(full_spec_files[i]).stem)[-15:]:
            print(str(i))
            orders_class_object=Datafile_apo_demetra_with_orders(zipfiles[i],full_spec_files[i],ll_file=linelist_file)
            testfile_list.append(orders_class_object)
        else:
            print('YOU SUCK')
    return testfile_list

def create_datafile_eshel(filelist=fl_eshel_all,i=0):
    mark1 = [1, 37]
    mark2 = [7, 10, 11, 32, 33, 34, 35, 36]
    mark3 = [8, 9]
    k=1
    for file in filelist:
        mark = 0
        if k in mark1:
            mark = 1
        elif k in mark2:
            mark = 2
        elif k in mark3:
            mark = 3
        a = Datafile_apo(file, i=k, mark=mark)
        dl, dl2 = airmass.split_date(a.header['DATE-OBS'])
        savename = datafile_folder_apo+a.observatory+'{num:02d}'.format(num=a.i)+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(a, workfileresource)
        workfileresource.close()
        print(k, mark)
        k+=1

def create_datafiles_lapalma(filelist=sortedfl_lapalma,save_folder=datafile_folder_merc,linelist_file=None):
    k=1
    for file in filelist:
        print(k)
        startdate = airmass.timeanddatelp(file)
        a = Datafile_mercator(file, i=k,ll_file=linelist_file)
        dl, dl2 = airmass.split_date(a.header['DATE-OBS'])
        savename = save_folder+a.observatory+'{num:02d}'.format(num=a.i)+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(a, workfileresource)
        workfileresource.close()
        k+=1


def create_datafiles_lapalma_omar(filelist=fl_dataset_omar,save_folder=datafile_folder_omar,linelist_file=None):
    k=1
    for file in tqdm.tqdm(filelist):
        startdate = airmass.timeanddatelp(file)
        a = Datafile_mercator_omar(file,i=k,ll_file=linelist_file)
        dl, dl2 = airmass.split_date(a.header['DATE-OBS'])
        savename = save_folder+a.observatory+'{num:03d}'.format(num=a.i)+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(a, workfileresource)
        workfileresource.close()
        k+=1

def make_data_grid(masterfilelist,line,v_min,v_max,rebin_size=0.5,selectionstring = 'All',selectionmode=None):
    linekey = line+'_original'
    snr_region = [5224, 5239]
    rebinv_lim = 1000
    firstfile = masterfilelist[0]
    lastfile = masterfilelist[-1]
    [yr_f,m_f,d_f,t_f]=airmass.split_date(firstfile.header['DATE-OBS'])[0]
    [yr_l, m_l, d_l, t_l] = airmass.split_date(lastfile.header['DATE-OBS'])[0]
    li = getattr(firstfile,linekey).lineinfo
    pi = [['v_min',v_min], ['v_max',v_max],['Rebin binsize (A)',rebin_size] ]
    si = [['Selection',selectionstring],['First observation date',yr_f+'-'+m_f+'-'+d_f],['Last observation date',yr_l+'-'+m_l+'-'+d_l],['Selection mode',selectionmode]]
    centerwl = li[1]
    rebinwl_lim = np.round(airmass.velocity_to_wl([-rebinv_lim,rebinv_lim],centerwl),decimals=1)
    wavenew = np.arange(rebinwl_lim[0],rebinwl_lim[1],rebin_size)
    snr_wavenew = np.arange(snr_region[0],snr_region[1],rebin_size)
    v_rebin = airmass.wl_to_velocity(wavenew, centerwl)[0]
    speed_index = np.where((v_rebin > v_min) & (v_rebin < v_max))
    wl_bound = wavenew[speed_index]
    speed_bound = v_rebin[speed_index]
    fluxarraylist = []

    bjdlist = []
    headerlist = []
    snrlist =  []
    for file in masterfilelist:
        full_wl = file.wl_original
        full_flux = file.flux_original
        snr_bound = np.where((full_wl > (snr_region[0]-1)) & (full_wl < (snr_region[1]+1)))
        snr_wl = full_wl[snr_bound]
        snr_flux = full_flux[snr_bound]
        header = file.header
        BJD = file.BJD
        wl = getattr(file,linekey).wl
        v= getattr(file,linekey).v
        bvcor_check = file.bc_from_header
        if bvcor_check is False:
            v = v-file.baricentric_correction
            wl = airmass.velocity_to_wl(v,centerwl)
        flux = getattr(file,linekey).flux
        flux_rebin = airmass.rebin_spec(wl,flux,wavenew)
        snr_flux_rebin = airmass.rebin_spec(snr_wl,snr_flux,snr_wavenew)
        snr = airmass.SNR_3(snr_wavenew,snr_flux_rebin,boundaries=snr_region,rebin=False,separate=False)
        snrlist.append(snr)
        flux_bound = flux_rebin[speed_index]
        bjdlist.append(BJD)
        headerlist.append(header)
        fluxarraylist.append(flux_bound)
    pi.append(['snr_average',np.average(snrlist)])
    datadict = dict(flux=fluxarraylist, wl=wl_bound, v=speed_bound, BJD=bjdlist, header=headerlist, snrlist=snrlist,
                    li=li, paraminfo=pi,selectioninfo =si)
    # datadict[flux]=fluxarraylist
    # datadict[wl]=wl_bound
    # datadict[v]=speed_bound
    # # datadict[bjd]=bjdlist
    # datadict[header]=headerlist
    return datadict


def make_data_grid_with_degradation(masterfilelist,line,v_min,v_max,R,snr_desired, selectionstring='All', selectionmode=None,rebin_size=0.5):
    linekey = line+'_original'
    snr_region = [5224, 5239]
    rebinv_lim = 1000
    firstfile = masterfilelist[0]
    lastfile = masterfilelist[-1]
    [yr_f,m_f,d_f,t_f]=airmass.split_date(firstfile.header['DATE-OBS'])[0]
    [yr_l, m_l, d_l, t_l] = airmass.split_date(lastfile.header['DATE-OBS'])[0]
    li = getattr(firstfile,linekey).lineinfo
    pi = [['v_min',v_min], ['v_max',v_max],['Rebin binsize (A)',rebin_size],['Spectral resolution',R],['Desired SNR',snr_desired]]
    si = [['Selection', selectionstring], ['First observation date', yr_f + '-' + m_f + '-' + d_f],
          ['Last observation date', yr_l + '-' + m_l + '-' + d_l], ['Selection mode', selectionmode]]
    centerwl = li[1]
    rebinwl_lim = np.round(airmass.velocity_to_wl([-rebinv_lim,rebinv_lim],centerwl),decimals=1)
    snr_wavenew = np.arange(snr_region[0], snr_region[1], rebin_size)
    wavenew = np.arange(rebinwl_lim[0],rebinwl_lim[1],rebin_size)
    v_rebin = airmass.wl_to_velocity(wavenew, centerwl)[0]
    speed_index = np.where((v_rebin > v_min) & (v_rebin < v_max))
    wl_bound = wavenew[speed_index]
    speed_bound = v_rebin[speed_index]
    fluxarraylist = []
    snrlist = []
    bjdlist = []
    headerlist = []
    for file in masterfilelist:
        full_wl = file.wl_original
        full_flux = file.flux_original
        snr_bound = np.where((full_wl > (snr_region[0]-1)) & (full_wl < (snr_region[1]+1)))
        snr_wl = full_wl[snr_bound]
        snr_flux = full_flux[snr_bound]
        snr_wl_deg, snr_flux_deg = airmass.degrade_spectrum_noise_first(snr_wl, snr_flux, spectral_resolution=R,
                                                                desired_snr=snr_desired, pre_rebin=None)
        header = file.header
        BJD = file.BJD
        wl = getattr(file,linekey).wl
        v= getattr(file,linekey).v
        bvcor_check = file.bc_from_header
        if bvcor_check is False:
            v=v-file.baricentric_correction
            wl = airmass.velocity_to_wl(v,centerwl)
        flux = getattr(file,linekey).flux
        wl_deg,flux_deg = airmass.degrade_spectrum_noise_first(wl,flux,spectral_resolution=R,desired_snr=snr_desired,pre_rebin=None)
        flux_deg_rebin = airmass.rebin_spec(wl_deg,flux_deg,wavenew)
        flux_deg_bound = flux_deg_rebin[speed_index]
        snr_flux_rebin = airmass.rebin_spec(snr_wl_deg, snr_flux_deg, snr_wavenew)
        snr = airmass.SNR_3(snr_wavenew,snr_flux_rebin,boundaries=snr_region,rebin=False,separate=False)
        snrlist.append(snr)
        bjdlist.append(BJD)
        headerlist.append(header)
        fluxarraylist.append(flux_deg_bound)
    datadict =  dict(flux = fluxarraylist, wl = wl_bound, v = speed_bound,BJD= bjdlist, header = headerlist,snrlist=snrlist,
                    li=li, paraminfo=pi,selectioninfo =si)
    # datadict[flux]=fluxarraylist
    # datadict[wl]=wl_bound
    # datadict[v]=speed_bound
    # # datadict[bjd]=bjdlist
    # datadict[header]=headerlist
    return datadict


def make_data_grid_apo(masterfilelist,line,v_min,v_max,rebin_size=0.5,selectionstring = 'All',selectionmode=None):
    linekey = line+'_original'
    snr_region = [5224, 5239]
    rebinv_lim = 1000
    firstfile = masterfilelist[0]
    lastfile = masterfilelist[-1]
    [yr_f,m_f,d_f,t_f]=airmass.split_date(firstfile.header['DATE-OBS'])[0]
    [yr_l, m_l, d_l, t_l] = airmass.split_date(lastfile.header['DATE-OBS'])[0]
    li = getattr(firstfile,linekey).lineinfo
    pi = [['v_min',v_min], ['v_max',v_max],['Rebin binsize (A)',rebin_size]]
    si = [['Selection', selectionstring], ['First observation date', yr_f + '-' + m_f + '-' + d_f],
          ['Last observation date', yr_l + '-' + m_l + '-' + d_l], ['Selection mode', selectionmode]]
    centerwl = li[1]
    rebinwl_lim = np.round(airmass.velocity_to_wl([-rebinv_lim,rebinv_lim],centerwl),decimals=1)
    wavenew = np.arange(rebinwl_lim[0],rebinwl_lim[1],rebin_size)
    snr_wavenew = np.arange(snr_region[0],snr_region[1],rebin_size)
    v_rebin = airmass.wl_to_velocity(wavenew, centerwl)[0]
    speed_index = np.where((v_rebin > v_min) & (v_rebin < v_max))
    wl_bound = wavenew[speed_index]
    speed_bound = v_rebin[speed_index]
    fluxarraylist = []
    bjdlist = []
    headerlist = []
    snrlist =  []
    for file in masterfilelist:
        full_wl = file.wl_original
        full_flux = file.flux_original
        snr_bound = np.where((full_wl > (snr_region[0]-1)) & (full_wl < (snr_region[1]+1)))
        snr_wl = full_wl[snr_bound]
        snr_flux = full_flux[snr_bound]
        header = file.header
        BJD = file.HJD
        wl = getattr(file,linekey).wl
        v= getattr(file,linekey).v
        flux = getattr(file,linekey).flux
        flux_rebin = airmass.rebin_spec(wl,flux,wavenew)
        snr_flux_rebin = airmass.rebin_spec(snr_wl,snr_flux,snr_wavenew)
        snr = airmass.SNR_3(snr_wavenew,snr_flux_rebin,boundaries=snr_region,rebin=False,separate=False)
        snrlist.append(snr)
        flux_bound = flux_rebin[speed_index]
        bjdlist.append(BJD)
        headerlist.append(header)
        fluxarraylist.append(flux_bound)
    pi.append(['snr_average',np.average(snrlist)])
    datadict = dict(flux=fluxarraylist, wl=wl_bound, v=speed_bound, BJD=bjdlist, header=headerlist, snrlist=snrlist,
                    li=li, paraminfo=pi,selectioninfo =si)
    # datadict[flux]=fluxarraylist
    # datadict[wl]=wl_bound
    # datadict[v]=speed_bound
    # # datadict[bjd]=bjdlist
    # datadict[header]=headerlist
    return datadict



def make_ls_brick(fluxbrick_filepath,frequencyarray = None):
    a = open(fluxbrick_filepath, 'rb')
    b = pickle.load(a)
    a.close()
    fluxgrid = np.array(b['flux'])
    fluxgrid_T = np.transpose(fluxgrid)
    BJDlist = np.array(b['BJD'])
    v = np.array(b['v'])
    header = b['header']
    lineinfo = b['li']
    snrlist = b['snrlist']
    li = b['li']
    pi = b['paraminfo']
    si = b['selectioninfo']
    power_ls_list = []
    min_freq = 1 / 10
    max_freq = 1/2
    # looping over single observation per v bin
    # Here use Lomb Scargle
    for single_vbin in fluxgrid_T:
        # print(single_vbin)
        # print(len(BJDlist), len(single_vbin), len(single_vbin_err))

        # autopower code
        if frequencyarray is None:
            frequency_LS, power_LS = LombScargle(BJDlist, single_vbin).autopower(
                minimum_frequency=min_freq,
                maximum_frequency=max_freq)
            if len(frequency_LS)<1000:
                frequency_LS = np.arange(1/10,1/2, 0.0001)
                power_LS = LombScargle(BJDlist, single_vbin).power(np.arange(1/10,1/2, 0.0001))
        else:
            frequency_LS = frequencyarray
            power_LS = LombScargle(BJDlist, single_vbin).power(frequencyarray)
        power_ls_list.append(power_LS)
    power_ls_array = np.asarray(power_ls_list)
    lombscl_dict = dict(powerarray=power_ls_array, frequency=frequency_LS, v=v, BJD=BJDlist, header=header, snrlist=snrlist,
                    li=li, paraminfo=pi,selectioninfo =si)
    # output_filepath=output_filefolder+lineinfo[0]+str(int(lineinfo[1]))+'rebin'+str(pi[2][1])+'_ls_brick.txt'
    # workfileresource = open(output_filepath, 'wb')
    # pickle.dump(lombscl_dict, workfileresource)
    # workfileresource.close()
    output_filename = lineinfo[0]+str(int(lineinfo[1]))+'rebin'+str(pi[2][1])+'_ls_brick.txt'
    return lombscl_dict, output_filename

def makenote(rebin='0.5'):
    a=rebin
    b = a.replace(".", "")
    notitie = 'rebin:' +rebin+ r'A, all darks, same night flat flatfields and same night bias, aD_csfF_snB_'+b
    return notitie
def run_cddo(snr, df,sf,note=0,filename_prefix='',rebin_size=0.5):
    # datafolder = str(Data_folder)+r'\Demetra\Individual\\'
    datafolder =df
    # savefolder = str(converted_Data_folder)+r'\demetra\with_orders\Individual\\'
    savefolder =sf
    linelist = str(converted_Data_folder)+r'\linelists\linelist_apo_v_cor_3.txt'
    create_datafiles_demetra_orders(datafolder,savefolder,linelist_file=linelist,snrtreshhold=snr,vshift=True,note=note,filename_prefix=filename_prefix,rebin_size=rebin_size)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin01\single_obs\\', note=makenote(rebin='0.1'),rebin_size=0.1)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\combined\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin01\combined\\', note=makenote(rebin='0.1'),rebin_size=0.1)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin02\single_obs\\', note=makenote(rebin='0.2'),rebin_size=0.2)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\combined\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin02\combined\\', note=makenote(rebin='0.2'),rebin_size=0.2)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin05\single_obs\\', note=makenote(rebin='0.5'),rebin_size=0.5)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\combined\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin05\combined\\', note=makenote(rebin='0.5'),rebin_size=0.5)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin025\single_obs\\', note=makenote(rebin='0.25'),rebin_size=0.25)
# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Zet_Ori_Data_Zet_Ori_Response3\final_spectra\combined\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\all_darks\rebin025\combined\\', note=makenote(rebin='0.25'),rebin_size=0.25)


# run_cddo(snr=None,df=str(Data_folder)+r'\Demetra\Individual\\',sf=str(converted_Data_folder)+r'\demetra\with_orders\Individual\\')
# run_cddo(snr=90)
def run_test_do():
    datafolder = str(Data_folder)+r'\Demetra\spectra_with_orders\\'
    linelist = str(converted_Data_folder)+r'\linelists\linelist_v_cor_2.txt'
    create_test_version_datafiles_demetra_orders(datafolder,linelist_file=linelist)
# a=run_test_do()
def make_testfiles_do():
    datafolder = str(Data_folder)+r'\Demetra\spectra_with_orders\\'
    savefolder = str(converted_Data_folder)+r'\demetra\with_orders\test\\'
    linelist = str(converted_Data_folder) + r'\linelists\linelist_v_cor_2.txt'
    create_datafiles_demetra_orders(datafolder, savefolder, linelist_file=linelist,snrtreshhold=100,vshift=False)
# make_testfiles_do()

def run_cdm():
    filelist = sortedfl_lapalma
    print('attention',filelist)
    savefolder = str(converted_Data_folder)+r'\mercator\ll_apo_vcor_2\\'
    linelist = str(converted_Data_folder)+r'\linelists\linelist_v_cor_2.txt'
    create_datafiles_lapalma(filelist=filelist,save_folder=savefolder,linelist_file=linelist)

def run_cdm_omar():
    filelist = fl_dataset_omar
    print('attention',filelist)
    savefolder = str(converted_Data_folder)+r'\dataset_omar\masterfiles_with_snr\\'
    linelist = str(converted_Data_folder)+r'\linelists\linelist_merc_incl_Hy.txt'
    create_datafiles_lapalma_omar(filelist=filelist,save_folder=savefolder,linelist_file=linelist)
run_cdm_omar()
# run_cdm()
# print(fl_apo_audela_all[7:10])
def run_cda():
    linelist = str(converted_Data_folder) + r'\linelists\linelist_apo_v_cor_2.txt'
    create_datafiles_audela(filelist=fl_apo_audela_all, savefolder=datafile_folder_audela_all,

                                linelist_file_path=linelist, vshift=True)
def run_mdg():
    filelist = open_masterfiles.mercator(str(converted_Data_folder) + r'\dataset_omar\\')
    linelist = open_masterfiles.open_linelist(str(converted_Data_folder) + r'\linelists\linelist_merc_incl_Hy.txt')
    savefolder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\no_degredation\rebin_05\\'
    rb=0.5
    for i,line in enumerate(linelist):
        vmin = -800
        vmax = 800
        linekey = 'line' + str(int(line[k]))
        print(i+1)
        data_grid = make_data_grid(filelist,linekey, vmin,vmax,rebin_size=rb)
        print(data_grid["paraminfo"])
        print(data_grid["snrlist"])
        savename = savefolder+'data_grid_'+line[0]+'_'+str(int(line[1]))+'rebin'+str(rb)+'_vlim'+str(vmin)+'_'+str(vmax)+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(data_grid, workfileresource)
        workfileresource.close()
def run_mdg_deg(R=10000,snr_desired = 1000):
    filelist = open_masterfiles.mercator(str(converted_Data_folder)+r'\dataset_omar\\')
    linelist = open_masterfiles.open_linelist(str(converted_Data_folder)+r'\linelists\linelist_merc_incl_Hy.txt')
    parent_path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\degraded\rebin_05'
    subforder = r'\R'+str(int(R))+'_snr'+str(int(snr_desired*5))+r'\\'
    savefolder = parent_path+subforder
    Path(savefolder).mkdir(parents=True, exist_ok=True)
    for i,line in enumerate(linelist):
        vmin = -800
        vmax = 800
        linekey = 'line' + str(int(line[k]))

        data_grid = make_data_grid_with_degradation(filelist,linekey, vmin,vmax,R=R,snr_desired=snr_desired,rebin_size=0.5)
        savename = savefolder+'data_grid_'+line[0]+'_'+str(int(line[1]))+str(vmin)+'_'+str(vmax)+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(data_grid, workfileresource)
        workfileresource.close()


def run_mdg_deg_spectra_removed(R=10000,snr_desired = 1000):
    filelist = open_masterfiles.mercator(str(converted_Data_folder)+r'\dataset_omar\\')
    linelist = open_masterfiles.open_linelist(str(converted_Data_folder)+r'\linelists\linelist_merc_incl_Hy.txt')
    parent_path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\degraded\rebin_05'
    subforder = r'\R'+str(int(R))+'_snr'+str(int(snr_desired*5))+r'\\'
    savefolder = parent_path+subforder
    Path(savefolder).mkdir(parents=True, exist_ok=True)
    for i,line in enumerate(linelist):
        vmin = -800
        vmax = 800
        linekey = 'line' + str(int(line[k]))

        data_grid = make_data_grid_with_degradation(filelist,linekey, vmin,vmax,R=R,snr_desired=snr_desired,rebin_size=0.5)
        savename = savefolder+'data_grid_'+line[0]+'_'+str(int(line[1]))+str(vmin)+'_'+str(vmax)+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(data_grid, workfileresource)
        workfileresource.close()


def run_mlb():
    input_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\no_degredation\rebin_01\\'
    output_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\original\custom_f_array\rebin01\\'
    fl = glob.glob(input_folder + r'*.txt')
    f_array = airmass.make_frequency_array(2,12,6.5,7,smallstep=1/100000)
    for filepath in tqdm.tqdm(fl):
        make_ls_brick(filepath,output_folder,frequencyarray=f_array)
def run_mlb_deg():
    input_base_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\degraded\rebin_05\\'
    output_base_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\degraded\rebin_05\\'
    #  --------- pick below to turn on folder snr selection
    input_folder_list = glob.glob(input_base_folder + r'\*')
    # folders = glob.glob(input_base_folder + r'\*')
    # input_folder_list = []
    # for folder in folders:
    #     match = re.search(r"snr(\d+)", folder)
    #     if match:
    #         x = int(match.group(1))
    #         if x < 50:
    #             input_folder_list.append(folder)
    # ---------
    print(input_folder_list)
    for folderpath in tqdm.tqdm(input_folder_list):
        subfolder = os.path.basename(folderpath)
        savefolder = output_base_folder+subfolder
        Path(savefolder).mkdir(parents=True, exist_ok=True)
        filelist = glob.glob(folderpath+'\*.txt')
        for filepath in filelist:
            lombscl_dict,output_filename = make_ls_brick(filepath)
            output_filepath = savefolder+r'\\'+output_filename
            workfileresource = open(output_filepath, 'wb')
            pickle.dump(lombscl_dict, workfileresource)
            workfileresource.close()

# run_mlb()
# run_mdg()
# snr_list = np.divide([1,2,5,10,20,30,40],5)
# # run_mdg_deg(R=10000,snr_desired=20)
# for snr in tqdm.tqdm(snr_list):
#     run_mdg_deg(R=10000,snr_desired=snr)

# run_mlb_deg()

# run_cda()
# create_datafiles_demetra(filelist=fl_demetra_good_alt,savefolder=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\altair_good\\',linelist_file_path=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\linelists\linelist_apo.txt')
# testfile_apo = fl_clean[12]
# a = Datafile_class.Datafile_apo(testfile_apo)
# testfile_merc = filelist_lapalma[3]
# a = Datafile_class.Datafile_mercator(testfile_merc, i=1)
# dl,dl2 = airmass.split_date(a.header['DATE-OBS'])
# testname = a.observatory+'{num:02d}'.format(num=a.i)+'_'+dl[0]+dl[1]+dl[2]+dl[3]+'.txt'
# print testname
#
# print a.i
# print a.phase
# print a.snr
# print a.line6562.ew
# for item in a.header:
#     print item
# print a.header['DATE-OBS']
# dl,dl2 = airmass.split_date(a.header['DATE-OBS'])
# print dl,dl2
# plt.plot(a.wl_original,a.flux_original,label = 'original')
# plt.plot(a.wl_rebin,a.flux_rebin, label = 'rebinned')
# plt.legend()
# plt.show()

# a = Datafile_class.Datafile_mercator()
# # make_line_data_dict(a.linelist,a.wl_rebin,a.flux_rebin,a.observatory)
# print a.filename
# savename = datafile_folder+a.filename+'.txt'
# pickle.dump(a, open(savename,'w'))

# line = ll_lapalma[3]
# wl = a.wl_rebin
# flux = a.flux_rebin
# obs = a.observatory
# line2,lw,lf,v,nf,vsini = line_data(line,wl,flux,obs)
# c = Line(line2,lw,lf,v,nf,vsini)


#
# b = pickle.load(datafile_folder+'test.txt','r')
# testfile = open(testsave,'r')
# b=pickle.load(testfile)
#
# print b.filename

# fl = glob.glob(r'D:\Peter\School\Master Thesis\Data\masterfiles\test\*.txt')
#
# testfile = open(fl[0],'r')
# b=pickle.load(testfile)
#
# print b.snr_original,b.snr

# plt.plot(b.line4026.v,b.line4026.flux)
# plt.show()
