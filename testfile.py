from __future__ import division
import glob
import pickle
from pathlib import Path
import numpy as np
import os
import airmass
import open_masterfiles
from astropy.time import Time
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
# os.environ['PYSYN_CDBS'] = 'C:\Users\Peter\Anaconda2\envs\p27\Lib\site-packages\pysynphot\data\cbds\grp\hst\cdbs'
import scipy
# filelist = open_masterfiles.mercator(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\')
# fl_apo = open_masterfiles.apo_demetra_orders()
# fi = fl_apo[0]
# ll = fi.available_lines
# textfilelist = os.listdir(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\no_degredation\rebin_05')
# tss_apo = []
# for line in ll:
#     li = getattr(fi,line).lineinfo
#     textstring = r'data_grid_' +li[0]+'_'+str(int(li[1]))+r'rebin0.5_vlim-800_800.txt'
#     tss_apo.append(textstring)
#
# for item in textfilelist:
#     if item not in tss_apo:
#         print(item)
# print(os.listdir(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\no_degredation\rebin_05'))

