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
import Make_Plots
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
apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
             'line4861', 'line4921', 'line6678', 'line4471']
# ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, r'He I 6678']
fl_eshel_clean_folder = r'D:\Peter\Master Thesis\Data\eShelData\data\clean' #Zonder twee spectra met rare Halpha spike
filelist_lapalma_folder = r'D:\Peter\Master Thesis\Data\LaPalmaData'
fl_eshel_all_folder = r'D:\Peter\Master Thesis\Data\eShelData\data\AlleSpectra'
fl_eshel_goodSNR_folder = r'D:\Peter\Master Thesis\Data\eShelData\data'
fl_all = glob.glob(fl_eshel_all_folder+r'\*.fit')
fl_clean = glob.glob(fl_eshel_clean_folder+r'\*.fit')
fl_goodSNR = glob.glob(fl_eshel_goodSNR_folder+r'\*.fit')
filelist_lapalma = glob.glob(filelist_lapalma_folder+r'\*.fits')
apo_lines2 = ['line6562','line4861']
apo_lines3 =[  'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
              'line4921', 'line6678', 'line4471']
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

pf1=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_new\all'
pf2=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_new\snr_100'
pf3=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_old\all'
pf4=r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\demetra_from_orders\ll_old\snr_100'

df1 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\\'
df2 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\snr_100\\'
df3 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\ll_oud\\'
df4 = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\ll_oud\snr_100\\'
pfs=[pf1,pf2,pf3,pf4]
dfs=[df1,df2,df3,df4]

# for i in range(len(pfs)):
#     datareduc.plot_TVS_orders(apo_lines, plot_save_folder=pfs[i], show='off', save='on', sg='off', oneline='on', siglvlline=0.01,datafilefolder=dfs[i], norm_boundaries='on')

datareduc.plot_TVS_orders(apo_lines3, plot_save_folder=pf2+r'\vrange1000', show='off', save='on', sg='off', oneline='on', siglvlline=0.01,datafilefolder=df2, norm_boundaries='on',style = None,vrange=1000)
# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data\clean',r'D:\Peter\Master Thesis\figures\TVS\eShel\cleanonly\reference_line',ll_TVS_eshel,show='off', save = 'on',sg='on',oneline='on')


# datareduc.plot_TVS_eShel('D:\Peter\Master Thesis\Data\eShelData\data',r'D:\Peter\Master Thesis\figures\TVS\eShel\every_good_snr',ll_TVS_eshel,show='off', save = 'on')
#
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