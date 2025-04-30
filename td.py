from __future__ import division
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import Path_check
import os
import pickle
import calendar
import numpy as np
# import airmass
from scipy.optimize import *
from scipy.stats import chi2
from PyAstronomy import pyasl
from astropy.timeseries import LombScargle
from matplotlib import ticker, cm
import airmass
import tqdm
folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)
[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
datafile_folder_omar = str(converted_Data_folder)+r'\dataset_omar\\'
# fl_dataset_omar = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Data\Dataset_Omar\fits\*.fits')
# converted_datafiles_omar = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\'+r'*.txt')
# datafiles = []
# for filepath in converted_datafiles_omar:
#     a = open(filepath, 'rb')
#     b = pickle.load(a)
#     datafiles.append(b)
#     a.close()
# i=0
# dd = airmass.make_data_grid(datafiles,'line4861',-150,150,rebin_size=0.1)
# for f in datafiles:
#     htd = f.header['DATE-OBS']
#     bjd = f.BJD
#     bjd_new = airmass.bjd_lapalma_from_date_zet_ori_omar(htd)
#     bjd_zelf = airmass.bjd_lapalma_from_date_zet_ori(htd)
#     print('file:',bjd,' Omar:',np.abs(bjd_new[0]-bjd),' Zelf:',np.abs(bjd_zelf-bjd))
# for file in datafiles:
#     print(i, file.filename, file.header['DATE-OBS'])
#     print(file.HJD)
#     print('-------------')
#     i+=1
# testfile=pf.open(fl_dataset_omar[0])
# header=testfile[0].header
# time = header['DATE-OBS']
# bjd1=airmass.bjd_lapalma_from_date_zet_ori(time)
# bjd2 =header['BJD']
# dif = (bjd1-bjd2)*24*60*60
# a=[1,2,3,4]
# b=str(a)
# print(b[0])
# testfile.close()
# print(bjd1)
# print(bjd2)
# print(dif)
# rng = np.random.default_rng()
#
# A = 2.  # amplitude
# c = 2.  # offset
#
# prot = 4
# w0 = 2*np.pi/prot  # rad/sec
# nin = 150
# nout = 1002
# x = rng.uniform(1, 20, nin)
# w = (2*np.pi)/np.linspace(1,10,1000)
# t= (2*np.pi)/w
# print(np.linspace(1,10,1000),t)
# y = A * np.cos(w0*x) + c
#
# # w = np.linspace(0.25, 10, nout)
# from scipy.signal import lombscargle
# pgram_power = lombscargle(x, y, w, normalize=False)
# plt.plot(t,pgram_power)
# plt.show()
# plt.close()
#
# a = np.linspace(1,10,10)
# b=1/a
# print(b)

# print header['RA'],header['DEC']
# def suck():
#     # header = pf.open(file)[0].header
#     JD = 2457306.585396
#     DEC = (-1.94)* (2*math.pi/360)
#     RA = 85.18* (2*math.pi/360)
#     LAT = 28.7636* (2*math.pi/360)
#     LON = -17.8947* (2*math.pi/360)
#     D = JD - 2451545.0
#     GMST = ((18.697374558 + 24.06570982441908*D)%24) * (360/24)* (2*math.pi/360)
#     sinalt = math.sin(DEC)*math.sin(LAT) + math.cos(DEC)*math.cos(LAT)*math.cos((GMST + LON - RA))
#     alt = np.arcsin(sinalt)*360/(2*math.pi)
#     airmass = 1/sinalt
#     return airmass,alt,JD
# am,alt,jd = suck()
# print alt
# datadict =  dict(flux = fluxarraylist, wl = wl_bound, v = speed_bound,BJD= bjdlist, header = headerlist)

file = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\data_grids\vlim-800_800\data_grid_Hy_4340-800_800.txt'
a = open(file, 'rb')
b = pickle.load(a)
a.close()
fluxgrid = np.array(b['flux'])
fluxgrid_T = np.transpose(fluxgrid)
BJDlist = np.array(b['BJD'])
v = np.array(b['v'])
header = b['header']
lineinfo = b['li']
power_ls_list = []
min_freq=1/10
max_freq=1
# looping over single observation per v bin
# Here use Lomb Scargle
for single_vbin in fluxgrid_T:
    # print(single_vbin)
    # print(len(BJDlist), len(single_vbin), len(single_vbin_err))

    # autopower code
    frequency_LS, power_LS = LombScargle(BJDlist, single_vbin).autopower(
        minimum_frequency=min_freq,
        maximum_frequency=max_freq)

    power_ls_list.append(power_LS)
power_ls_array = np.asarray(power_ls_list)
lombscl_dict = [power_ls_array, frequency_LS, v, BJDlist]


def period_plotter(line_period_info):

    # line_number = star_pickle['roundline'].values
    # plot_titles = star_pickle['plottitle'].values
    #
    # infile_LS = open(pickle_location, "rb")
    # period_dict = pickle.load(infile_LS)
    #
    # line_period_info = period_dict[line]

    power_ls = line_period_info[0]
    frequency_ls = line_period_info[1]
    wave_grid = line_period_info[2]
    BJDlist = line_period_info[3]

    x_wave, y_freq = np.meshgrid(wave_grid, 1 / frequency_ls)
    power_ls_trans = power_ls.T

    # print(min(list(map(min, power_ls))))
    # print(max(list(map(max, power_ls))))

    som_frequency = np.sum(power_ls_trans, axis=1)
    som_wave_place = np.sum(power_ls_trans, axis=0)
    som_wave = som_wave_place / np.max(som_wave_place)

    fig5 = plt.figure(figsize=(15, 13.333))
    widths = [15, 5]
    heights = [5, 20]
    spec5 = fig5.add_gridspec(ncols=2, nrows=2, width_ratios=widths, height_ratios=heights, \
                              wspace=0.0, hspace=0.0)

    ax1 = fig5.add_subplot(spec5[0, 0])
    # ax1.set_title(fr'${plot_titles[index]} \ {line_number[index]}$', fontsize=26, pad=10)

    ax1.plot(wave_grid, som_wave, lw=1.0, color='k')
    ax1.set_xlim(np.min(wave_grid), np.max(wave_grid))
    ax1.tick_params(labelsize=23, bottom=False, left=False, right=True)
    ax1.tick_params(axis='y', which='major', length=8)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")

    ax2 = fig5.add_subplot(spec5[1, 0])
    ax2.set_yscale('log')
    cs = ax2.contourf(x_wave, y_freq, power_ls_trans, cmap=cm.YlOrBr)
    ax2.set_ylabel("Period (d)", fontsize=25)
    ax2.set_xlabel("Velocity (km/s)", fontsize=25)
    ax2.tick_params(labelsize=23)
    ax2.tick_params(which='major', length=8)
    ax2.tick_params(which='minor', length=4)
    # ax2.yaxis.set_minor_formatter(NullFormatter())

    ax3 = fig5.add_subplot(spec5[1, 1])
    ax3.set_yscale('log')
    ax3.plot(som_frequency, 1 / frequency_ls, lw=0.5, color='k')
    ax3.set_ylim(np.min(1 / frequency_ls), np.max(1 / frequency_ls))
    cbar = fig5.colorbar(cs)

    # set ticks left and right, turn ylabels off
    ax3.tick_params(axis='y', which='both', direction='in', left=True, right=True)
    ax3.tick_params(axis='y', which='major', length=8)
    ax3.tick_params(axis='y', which='minor', length=4)
    # ax3.yaxis.set_minor_formatter(NullFormatter())
    ax3.tick_params(axis='x', direction='in', labelsize=23)
    ax3.set_yticklabels([])
    plt.show()
    # plt.savefig(fname=source_location + '/rel_flux/' + line + '_' + star_name + '.png', format='png', \
    #             bbox_inches='tight', dpi=200)
    plt.clf()
    del (cs, power_ls_trans)

    return

period_plotter(lombscl_dict)