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
import csv
import os
import open_masterfiles
import Path_check
from collections import defaultdict

folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)

[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
datafolder_omar = str(converted_Data_folder)+r'\dataset_omar\\'
def make_ew_file_normal(folder,target_folder_path,linelist,binsize='05',vlim=[-500, 500], ha_vlim=[-500,500]):

    for line in linelist:
        if binsize == 'original':
            lk = line+'_original'
        elif binsize in ['005','01','02','025','03','05']:
            lk=line+'_rebin'+binsize
        else:
            lk = line+'_original'
        ews, hjds, phases, ers = datareduc.equivalent_width_array_mercator(folder, line, vlim=vlim,ha_vlim= ha_vlim,binsize=binsize)
        datastructure = ['Header','Equivalent Widths','HJDs','Phases','Errors']
        header = [line, folder, binsize, str(vlim),datastructure]
        data = [ews,hjds,phases,ers]
        with open(target_folder_path+lk+'_EW.csv', 'w',newline='') as f:
            write = csv.writer(f)
            write.writerow(header)
            write.writerows(data)

def open_ew_file(filepath):
    with open(filepath) as fp:
        reader = csv.reader(fp, delimiter=",", quotechar='"')
        # next(reader, None)  # skip the headers
        data_read = [row for row in reader]
        header = data_read[0]
        data = data_read[1:]
        x = np.array(data)
        y = x.astype(float)
        converted_data_list = y.tolist()

    return header,converted_data_list

# r'D:\peter\Master_Thesis\Datareduction\Converted_Data\test\ew\test.csv'

# make_ew_file_normal(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\mercator\ll_apo_vcor_2',r'D:\peter\Master_Thesis\Datareduction\Converted_Data\test\ew\\',linelist= ['line6562'],binsize='original')
make_ew_file_normal(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar',r'D:\peter\Master_Thesis\Datareduction\Converted_Data\test\ew\\omar_',linelist= ['line4340'],binsize='original')



data = open_ew_file(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\test\ew\omar_line6562_original_EW.csv')

# for item in data:
#     print(item)
ews = data[1][0]
hjds = data[1][1]

datareduc.LS_periodogram_from_EW(hjds,ews,r'H\alpha',searchrange=[1,20])