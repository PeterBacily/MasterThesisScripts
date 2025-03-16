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
def make_ew_file_normal(folder,target_file_path,linelist):

    for line in linelist:
        header = [line,folder,]
        ews, hjds, phases, ers = datareduc.equivalent_width_array_mercator(folder, line, vlim=[-500, 500])
        data = [ ]
        with open(target_file_path, 'w') as f:
            write = csv.writer(f)
            write.writerow(header)
            write.writerows(data)