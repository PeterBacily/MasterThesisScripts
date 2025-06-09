from __future__ import division
import glob
import pickle
from pathlib import Path
import numpy as np
import os
import open_masterfiles
from astropy.time import Time
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
# os.environ['PYSYN_CDBS'] = 'C:\Users\Peter\Anaconda2\envs\p27\Lib\site-packages\pysynphot\data\cbds\grp\hst\cdbs'
import scipy
# beginning = 1
# end=12
# bigstep = 1
# smallstep = 1/100
# ss_start = 8
# ss_end=9
# a = np.arange(beginning,ss_start,bigstep)
# b=np.arange(ss_start,ss_end,smallstep)
# c= np.append(a,b)
# d =np.arange(ss_end,end+bigstep,bigstep)
# e=np.append(c,d)
# print(e)
#
# apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
#                  'line4861', 'line4921', 'line6678', 'line4471']
# mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
#                       'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']

def group_observations(objects, gap_days=183):
    # Sort objects by HJD
    sorted_objects = sorted(objects, key=lambda obj: obj.HJD)

    groups = []
    current_group = []

    previous_hjd = None

    for obj in sorted_objects:
        if previous_hjd is None:
            # Start the first group
            current_group = [obj]
        else:
            if obj.HJD - previous_hjd >= gap_days:
                # Gap is large enough to start a new group
                groups.append(current_group)
                current_group = [obj]
            else:
                current_group.append(obj)
        previous_hjd = obj.HJD

    # Append the last group
    if current_group:
        groups.append(current_group)

    return groups

# def group(files)
#     complete_file_list = open_masterfiles.mercator(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\')
#     hjdlist = []
#     for file in complete_file_list:
#         hjdlist.append(file.HJD)
#     t=Time(hjdlist,format='jd')
#     dtt=t.datetime
#     dtt.sort()
#
#     # Group the datetimes based on the 6-month gap rule
#     groups = []
#     current_group = [dtt[0]]
#
#     for i in range(1, len(dtt)):
#         # Check the difference between current datetime and the last datetime in the group
#         if dtt[i] - current_group[-1] >= timedelta(days=183):  # 6 months = 183 days
#             groups.append(current_group)
#             current_group = [dtt[i]]
#         else:
#             current_group.append(dtt[i])
#
#     # Add the last group
#     groups.append(current_group)
complete_file_list = open_masterfiles.mercator(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\')
groups = group_observations(complete_file_list)
# Print the results
for group in groups:
    # for file in group:
    #     dto = Time(file.HJD, 'jd').datetime
    #     print(dto.strftime("%Y-%m-%d"))
    print([Time(file.HJD,format='jd').datetime.strftime("%Y-%m-%d") for file in group])
# plt.scatter(dtt,np.zeros(len(dtt)))
# plt.show()
# plt.close()
# a = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator\degraded\rebin_05\*')
# for folder in a:
#     filepaths = glob.glob(folder+r'\*.txt')
#     fp = filepaths[6]
#     testfile = open(fp, 'rb')
#     datadict = pickle.load(testfile)
#     powerarray = datadict['powerarray']
#     fq= datadict['frequency']
#     v= datadict['v']
#     BJD = datadict['BJD']
#     header = datadict['header']
#     snrlist = datadict['snrlist']
#     lineinfo = datadict['li']
#     paraminfo = datadict['paraminfo']
#     print(lineinfo)
#     print(paraminfo[-1],str(paraminfo[-1][1]*6), np.average(snrlist))

    # for filepath in filepaths:
    #     path = Path(filepath)
    #     subdirectory = path.parent.name
    #     print(subdirectory)
# a,b,c = os.walk(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\degraded')[0]
# print(a)
# print(b)
# print(c)