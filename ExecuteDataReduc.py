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
import datareduc
import os
import Make_Plots
import open_masterfiles
from collections import defaultdict
import tqdm
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
ll_TVS_eshel = [['He', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
ll_nieuwe_lijnen = [['He_I', 35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922']]

# ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, r'He I 6678']
fl_eshel_clean_folder = r'D:\Peter\Master Thesis\Data\eShelData\data\clean' #Zonder twee spectra met rare Halpha spike
filelist_lapalma_folder = r'D:\Peter\Master Thesis\Data\LaPalmaData'
fl_eshel_all_folder = r'D:\Peter\Master Thesis\Data\eShelData\data\AlleSpectra'
fl_eshel_goodSNR_folder = r'D:\Peter\Master Thesis\Data\eShelData\data'
fl_all = glob.glob(fl_eshel_all_folder+r'\*.fit')
fl_clean = glob.glob(fl_eshel_clean_folder+r'\*.fit')
fl_goodSNR = glob.glob(fl_eshel_goodSNR_folder+r'\*.fit')
filelist_lapalma = glob.glob(filelist_lapalma_folder+r'\*.fits')

apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
             'line4861', 'line4921', 'line6678', 'line4471']
line_ha ='line6562'
apo_lines2 = ['line6562']
apo_lines3 =[ 'line5875','line6562','line4861',  'line6678' ]
ha_hb_linelist = ['line6562','line4861']
lines_III = ['line6562','line4861', 'line4685']
# esb=0
# est=-0.05
# # shw,sv = 'on', 'off'
# shw,sv = 'off', 'on'

import Path_check

import pickle
# import seaborn as sns
folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)

[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
# ----------------------
#
#
# TVS La Palma
# datareduc.plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma\test',[ll_lapalma[-1]],show='on', save = 'off',sg='on',oneline='on')
# datareduc.plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma',ll_lapalma,      show='off', save = 'on',sg='on',oneline='on')
# datareduc.plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma\smoothed',         ll_lapalma,show='off', save = 'on', sg='on',  oneline='off')
# datareduc.plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma\reference_line',   ll_lapalma,show='off', save = 'on', sg='off', oneline='on')
# datareduc.plot_TVS_Lapalma('D:\Peter\Master Thesis\Data\LaPalmaData',r'D:\Peter\Master Thesis\figures\TVS\LaPalma\both',             ll_lapalma,show='off', save = 'on', sg='on',  oneline='on')

#
#
# TVS eShel
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data\clean',r'D:\Peter\Master Thesis\figures\TVS\eShel\cleanonly\smoothed',         ll_TVS_eshel,show='off', save = 'on', sg='on',  oneline='off')
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data\clean',r'D:\Peter\Master Thesis\figures\TVS\eShel\cleanonly\reference_line',   ll_TVS_eshel,show='off', save = 'on', sg='off', oneline='on')
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data\clean',r'D:\Peter\Master Thesis\figures\TVS\eShel\cleanonly\both',             ll_TVS_eshel,show='off', save = 'on', sg='on',  oneline='on')
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\TVS\eShel\every_good_snr\smoothed',         ll_TVS_eshel,show='off', save = 'on', sg='on',  oneline='off')
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\TVS\eShel\every_good_snr\reference_line',   ll_TVS_eshel,show='on', save = 'off', sg='off', oneline='on')
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\TVS\eShel\every_good_snr\both',             ll_TVS_eshel,show='off', save = 'on', sg='on',  oneline='on')

# TVS eShel Demetra
# datareduc.plot_TVS_eShel_masterfile(apo_lines, plot_save_folder=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\apo\demetra_altair_snr100',show='off',save='on',sg='off',oneline='on', siglvlline=0.01,datafilefolder=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\altair_good\snr100\\',datareductionprogram='Demetra', norm_boundaries='on')
# datareduc.plot_TVS_orders(apo_lines, plot_save_folder=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\apo\from_orders',show='off',save='on',sg='off',oneline='on', siglvlline=0.01, norm_boundaries='on')

# pf1=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_new\all'
# pf2=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_new\snr_100'
# pf3=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_old\all'
# pf4=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_old\snr_100'
#
# df1 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\v_cor\\'
# df2 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\v_cor\snr_100\\'
# df3 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\ll_oud\\'
# df4 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\ll_oud\snr_100\\'



pf1=str(Plots_folder)+r'\TVS\demetra_from_orders\ll_new\all'
pf2=str(Plots_folder)+r'\TVS\demetra_from_orders\ll_new\snr_100'
pf3=str(Plots_folder)+r'\TVS\demetra_from_orders\ll_old\all'
pf4=str(Plots_folder)+r'\TVS\demetra_from_orders\ll_old\snr_100'


pf_lp=str(Plots_folder)+r'\TVS\mercator\cropped\vrange1000\ll_apo_vcor_2'
pf_lp_final = str(Plots_folder)+r'\TVS\mercator\final'
pf_dem_final_all = str(Plots_folder)+r'\TVS\demetra_from_orders\final\all\\'
pf_dem_final_90 = str(Plots_folder)+r'\TVS\demetra_from_orders\final\snr90\\'
pf_dem_final_100 = str(Plots_folder)+r'\TVS\demetra_from_orders\final\snr100\\'
pf_full_night_100 = r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\full_night\SNR100\\'
pf_full_night_110 = r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\full_night\SNR110\\'
pf_full_night_all = r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\full_night\SNRALL\\'
df1 = str(converted_Data_folder)+r'\demetra\with_orders\v_cor_3\\'
df2 = str(converted_Data_folder)+r'\demetra\with_orders\v_cor_3\snr_90\\'
df3 = str(converted_Data_folder)+r'\demetra\with_orders\v_cor_3\snr_100\\'
# df4 = str(converted_Data_folder)+r'\demetra\with_orders\ll_oud\snr_100\\'
df_lp = str(converted_Data_folder)+r'\mercator\ll_apo_vcor_2\\'
df_test= str(converted_Data_folder)+r'\demetra\with_orders\test\snr_100\\'
data_full_night_all= str(converted_Data_folder)+r'\demetra\with_orders\all_darks\full_night\\'
# data_full_night_all= str(converted_Data_folder)+r'\demetra\with_orders\full_night\\'
data_full_night_110= str(converted_Data_folder)+r'\demetra\with_orders\full_night\snr_110\\'
data_full_night_100= str(converted_Data_folder)+r'\demetra\with_orders\full_night\snr_100\\'
# data_individual = str(converted_Data_folder)+r'\demetra\with_orders\Individual\\'
data_individual = str(converted_Data_folder)+r'\demetra\with_orders\all_darks\single_obs\\'
# data_individual_list = open_masterfiles.apo_demetra_orders(path = data_individual,manual_filelist=None,sort_data_files='on')
# data_full_list = open_masterfiles.apo_demetra_orders(path = data_full_night_all,manual_filelist=None,sort_data_files='on')
# pfs=[pf_dem_final_all,pf_dem_final_90,pf_dem_final_100]
# dfs=[df1,df2,df3]
dfs_full_night = [data_full_night_110,data_full_night_100]
pfs_full_night = [pf_full_night_110,pf_full_night_100]
audela_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\AudeLA\all\\'
snr_comp_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\SNR_comp\comp\\'

def get_day_data_list(data_individual_folder,data_full_night_path):

    data_full_list = open_masterfiles.apo_demetra_orders(path = data_full_night_path,manual_filelist=None,sort_data_files='on')
    data_individual_list= open_masterfiles.apo_demetra_orders(path = data_individual_folder,manual_filelist=None,sort_data_files='on')
    groups = defaultdict(list)
    for obj in data_individual_list:
        # print(obj.time_and_date)
        groups[obj.time_and_date[0:5]].append(obj)
    new_list = groups.values()
    list_of_day_data = list(new_list)
    # for file in data_full_list:
    #     print(file.time_and_date[0:6])
    daydatalist = []
    for day in list_of_day_data:
        tad= day[0].time_and_date[0:6]
        full_day_file=[x for x in data_full_list if x.time_and_date[0:6] == tad][0]
        daydatalist.append([day,full_day_file])
    return daydatalist
    # files_snr_test= open_masterfiles.apo_demetra_orders(path=snr_comp_folder)
    # for day in list_of_day_data:
    #     li=day[0].line6562.lineinfo
    #     datareduc.plot_snr_test(day,li[2:6] )


# ls_databrick_original = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\original\\'
# ls_databrick_original = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\original\rebin01\\'
ls_databrick_original = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\original\custom_f_array\rebin01\\'
# ls_databrick_original = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\original\rebin05\\'
ls_databrick_filelist = glob.glob(ls_databrick_original+r'*.txt')
# datareduc.ls_sum_plotter(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\original\rebin01\selection\\',-500,500,show='on',save='off')
datareduc.ls_sum_plotter(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\original\rebin05\\',-500,500,show='on',save='off',SG=False,SGwindowsize=31)

quit()
print(ls_databrick_filelist)
for filepath in tqdm.tqdm(ls_databrick_filelist):
    datareduc.ls_brick_plotter(filepath,-500,500,plotsavefolder=r'D:\peter\Master_Thesis\Datareduction\Plots\LS_periodogram\mercator_original\rebin01\\',show='on',save='off' )

quit()
print(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin01\single_obs\\')
rebin_base_path= r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin'
rebin_size = '01'
snr110=False
single_obs_suffix = r'\single_obs\\'
if snr110 is True:
    full_night_suffix =r'\combined\high_snr\\'
else:
    full_night_suffix = r'\combined\\'
di_path = rebin_base_path+rebin_size+single_obs_suffix
fn_path = rebin_base_path+rebin_size+full_night_suffix
ddl= get_day_data_list(di_path,fn_path)
mercator_files = open_masterfiles.mercator(df_lp)
for day in ddl:
    indiv = day[1]
    fn = day[0]
    datareduc.rebin_and_overplot_demetra_orders(fn,indiv,mercator_files[0],rebin_size=0.5,boundaries=[5310,5370])
# datareduc.plot_SNR_orders(['line4861'],indiv, file_full_night = fn,plot_avg=True,rebin=0.5,plot_save_folder=r'D:\peter\Master_Thesis\Datareduction\Plots\SNR\aD_snfF',show='on',save='off', norm_boundaries='on',vrange=1000,subplotylim=[0.97,1.03])

# for line in apo_lines:
#     datareduc.LS_periodogram_merc(df_lp,line,searchrange=[1,8])
# ews,hjds,phases,errs = datareduc.equivalent_width_array_mercator(df_lp,'line4861')
# pr,pdg = airmass.ls_periodogram(hjds,ews,searchrange=[1,8])
# plt.plot(pr,pdg)
# plt.show()
# plt.close()
# snr_orders(di_path,fn_path,show='on',save='off')

# data_individual_list = open_masterfiles.apo_demetra_orders(path=di_path, manual_filelist=None, sort_data_files='on')
# datareduc.plot_order_stack(data_individual_list,wlpiece= [5315, 5365],rebinstep=0.5,day=0,from_order=False)

# datareduc.plot_TVS_orders_lines_together(lines_III,r'D:\peter\Master_Thesis\Datareduction\Plots\test\demetra_tvs\\',show='off',save='on', norm_boundaries='off', datafilefolder=fn_path,vrange=800)

# datareduc.plot_TVS_together(linelist=apo_lines3,filefolder_apo=data_full_night_100,filefolder_merc=df_lp,show='on',save='off')
# for i in range(len(dfs_full_night)):
# datareduc.plot_TVS_orders(apo_lines3, plot_save_folder=pf_full_night_all, show='on', save='off', sg='off', oneline='on', siglvlline=0.01,datafilefolder=data_full_night_all, norm_boundaries='on',vrange=1000,from_order=False)

# datareduc.plot_TVS_orders(apo_lines, plot_save_folder=pf2+r'\vrange1000', show='on', save='off', sg='off', oneline='on', siglvlline=0.01,datafilefolder=df2, norm_boundaries='on',style = None,vrange=1000, from_order=True)
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data\clean',r'D:\Peter\Master Thesis\figures\TVS\eShel\cleanonly\reference_line',ll_TVS_eshel,show='off', save = 'on',sg='on',oneline='on')
# datareduc.plot_TVS_Lapalma_masterfile(apo_lines2,plot_save_folder=pf_lp_final,datafilefolder=df_lp,show='off', save='on', sg='off', oneline='on', siglvlline=0.01, norm_boundaries='on',style = None,vrange=1000)

# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\TVS\eShel\every_good_snr',ll_TVS_eshel,show='off', save = 'on')
# datareduc.plot_TVS_eShel_masterfile(linelist=apo_lines3, plot_save_folder=pf1,show='on',save='off',sg='off',oneline='off', siglvlline=0.01,datafilefolder=audela_folder,datareductionprogram='AudeLA',norm_boundaries='on',vrange=1000)
#
# Quotient eShel
# datareduc.plot_quotient_eShel('D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\test',[ll_TVS_eshel[1]],overplot='off',show='off',save='on',sg = 'off',oneline = 'on')
# for file in fl_clean:
#     print file
# print '--------------'
# for i in range(len(fl_clean)-1):
#     print fl_clean[i],fl_clean[i+1]

# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\TVS\eShel\every_good_snr\reference_line',   ll_nieuwe_lijnen,show='on', save = 'off', sg='off', oneline='on')

# Make_Plots.plot_EW_demetra(obs='MERCATOR', orders=True, figsavefolder=r'D:\peter\Master_Thesis\Datareduction\Plots\EW_from_orders\\',
#             custom_lines=apo_lines, custom_files=None, Chunklength=3,
#             datafile_folder_demetra =r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\snr_100\\',
#             datafile_folder_mercator=r'D:\peter\Master_Thesis\Datareduction\Converted_Data\mercator\test\\',
#                     save=True, show=False)

# datareduc.create_JoO_apo_demetra(data_individual,stacked=False)
# datareduc.create_JoO_mercator(df_lp)