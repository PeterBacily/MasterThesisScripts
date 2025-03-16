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
import pickle
import os
import Datafile_class
matplotlib.style.use('classic')
import open_masterfiles
import Path_check

folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)

[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
filelist_lapalma_folder = str(Data_folder)+r'\LaPalmaData'
filelist_lapalma = glob.glob(filelist_lapalma_folder+r'\*.fits')

# merc_linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
# apo_linelist = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
merc_linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4016, 4020, 4032, 4036, 'He I 4026'], ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'], ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5884.6, 5885.5, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
apo_linelist = [['Ha', 35, 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]


Sil_Omar_normline_list= [['Ha', 6562.819, [[-1300, -1000], [1700, 1850]], r'H$\alpha$ 6563'], ['Hb', 4861.333, [[1387, 1927], [-1000, -700]], r'H$\beta$ 4861'], ['Hy', 4340.472, [[-1000, -700], [1000, 1500]], r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4016, 4020, 4032, 4036, 'He I 4026'], ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'], ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5884.6, 5885.5, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]

ll_new = open_masterfiles.open_linelist(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\linelists\linelist_apo.txt')


Sil_Omar_normline_list_apo = [['Ha', 35, 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563', [[-1300, -1000], [1700, 1850]]], ['Hb', 35, 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861',[[-1886, -1252], [1387, 1927]]], ['He_I', 35, 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713',[[-1400, -1100], [1380, 1900]]], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876',[[-1870, -1300], [1810, 1891]]], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542',[[-1444, -1303], [1305, 1500]]], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686',[[-430, -167], [237, 535]]], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412',[[-900, -500], [700, 1000]]], ['He_I', 35,4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471',[[-1967, -1694], [1545, 1936]]]]

Sil_Omar_normline_list_apo_new_short= [['Ha', 35, 6562.819, 6549, 6550.7, 6576.0, 6578.0, r'H$\alpha$ 6563', [[-1300, -1000], [1700, 1850]]], ['Hb', 35, 4861.333, 4847.3, 4848.3, 4876.3, 4877.3, r'H$\beta$ 4861',[[-1886, -1252], [1387, 1927]]],['Hy',4340.472,4326.0, 4330.3, 4355.0, 4362.2,r'H$\gamma$ 4861', [[-1000, -700], [1000, 1500]]]]
testlist = [['Ha', 35, 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563', [[-1300, -1000], [1700, 1850]]]]


    # [ ['Hb', 4861.333, [[-1000, -700],[1387, 1927]], r'H$\beta$ 4861'], ['Ha', 6562.819, [[-1300, -1000], [1700, 1850]], r'H$\alpha$ 6563'],['Hy', 4340.472, [[-1000, -700], [1000, 1500]], r'H$\gamma$ 4340']]
Sil_Omar_normline_list_test = [ ['Hb', 4861.333, [[-1000, -700],[1387, 1927]], r'H$\beta$ 4861']]

def normalize(wave,flux,a,b,c,d,startwl,endwl):
    normwave = np.hstack((wave[(wave>a)&(wave<b)],wave[(wave>c)&(wave<d)]))
    normflux = np.hstack((flux[(wave>a)&(wave<b)],flux[(wave>c)&(wave<d)]))
    #fit line trough slice
    slope,height = np.polyfit(normwave,normflux,1)
    # print 'slope and height are', slope, height
    fit = np.poly1d([slope,height])
    linewave = wave[(wave>startwl)&(wave<endwl)]
    # print wave[0],startwl,wave[-1],endwl
    # print wave
    lineflux = flux[(wave>(startwl))&(wave<(endwl))]
    normlineflux = []
    # print linewave
    for i,j in enumerate(linewave):
        normlineflux.append(lineflux[i]/fit(j))
    fluxarray = np.array(normlineflux)
    nnf = []
    for k, nwl in enumerate(normwave):
        nnf.append(normflux[k]/fit(nwl))
    return linewave, fluxarray, nnf, lineflux, fit
# apo_files = open_masterfiles.apo()
# # print apo_files
# merc_files = open_masterfiles.mercator()
# merc_file =apo_files[1]
# line = ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861']
def plot_norm_sil_omar(linelist,obs='apo',show='on',save = 'off',datafilefolder_apo=None, plot_save_folder = r'D:\peter\Master_Thesis\Datareduction\Plots\normalization\test'):
    if datafilefolder_apo ==None:
        apo_files = open_masterfiles.apo_demetra_orders()
    else:
        apo_files = open_masterfiles.apo_demetra_orders(datafilefolder_apo)
    # print apo_files
    merc_files = open_masterfiles.mercator()
    if obs == 'apo':
        k = 8
        filelist=apo_files
    elif obs == 'mercator':
        k = 8
        filelist=merc_files
    else:
        raise(Exception("obs needs to be 'apo' or 'mercator'"))
    i=0
    for file in filelist:
        i+=1
        for line in linelist:
            # line_oud = apo_linelist[i]
            a = line[k][0][0]
            b = line[k][0][1]
            c = line[k][1][0]
            d = line[k][1][1]
            spacing = 200
            start = a-spacing
            stop = d+spacing

            wl = file.wl_rebin
            flux = file.flux_rebin
            print(wl[0])

            lw,nlf,nnf,lf,fit = airmass.normalize(wl,flux,a,b,c,d,start,stop,xtype='velo',linecenter=line[2])
            cond1 = (lw>a)&(lw<b)
            cond2 = (lw>c)&(lw<d)
            normwave1 = lw[cond1]
            normflux1= lf[cond1]
            normwave2 = lw[cond2]
            normflux2= lf[cond2]
            normflux3 = lf[lw<a]
            normwave3 = lw[lw<a]
            normflux4 = lf[(lw>b)&(lw<c)]
            normwave4 = lw[(lw>b)&(lw<c)]
            normflux5 = lf[lw>d]
            normwave5 = lw[lw>d]
            # lw_oud,nlf_oud,nnf_oud,lf_oud,fit_oud = normalize(wl,flux,line_oud[k+1],line_oud[1][0][1],line_oud[1][1][0],line_oud[1][1][1],line_oud[k+1]-20,line_oud[1][1][1]+20)
            # cond1_oud = (lw>line_oud[k+1])&(lw<line_oud[1][0][1])
            # cond2_oud = (lw>line_oud[1][1][0])&(lw<line_oud[1][1][1])
            # normwave1_oud = lw[cond1_oud]
            # normflux1_oud= lf[cond1_oud]
            # normwave2_oud = lw[cond2_oud]
            # normflux2_oud= lf[cond2_oud]
            # ((lw>c)&(lw<d))]
            # normflux = lf[((lw>line[k+1])&(lw<b)) or ((lw>c)&(lw<d))]
            f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True,layout="constrained")
            # ax1.plot(lw,lf)
            ax1.plot(normwave1,normflux1, c='r', label = 'continuum data slice')
            ax1.plot(normwave2,normflux2, c='r')
            ax1.plot(normwave3,normflux3, c='b', label = 'data')
            ax1.plot(normwave4,normflux4, c='b')
            ax1.plot(normwave5,normflux5, c='b')
            # ax1.plot(normwave1_oud,normflux1_oud,c='g')
            # ax1.plot(normwave2_oud, normflux2_oud, c='g')
            # ax1.plot(lw, fit_oud(lw), c='g', label='fitted continuum')
            ax1.plot(lw,fit(lw),c='red', linestyle='--', label = 'fitted continuum')
            # ax1.legend(prop={'size': 10}, loc='best')
            # ax2.plot
            ax2.plot(lw,nlf,c='r')
            ax1.set_title('Sil en omar \n Raw Flux')
            ax2.set_title('Normalized Flux')
            [a,b,c,d],vsini = airmass.wl_to_velocity([line[3],line[4],line[5],line[6]],line[2])
            # spacing = 200
            # start = a - spacing
            # stop = d + spacing

            wl = file.wl_rebin
            flux = file.flux_rebin
            print(wl[0])

            lw, nlf, nnf, lf, fit = airmass.normalize(wl, flux, a, b, c, d, start, stop, xtype='velo', linecenter=line[2])
            cond1 = (lw > a) & (lw < b)
            cond2 = (lw > c) & (lw < d)
            normwave1 = lw[cond1]
            normflux1 = lf[cond1]
            normwave2 = lw[cond2]
            normflux2 = lf[cond2]
            normflux3 = lf[lw < a]
            normwave3 = lw[lw < a]
            normflux4 = lf[(lw > b) & (lw < c)]
            normwave4 = lw[(lw > b) & (lw < c)]
            normflux5 = lf[lw > d]
            normwave5 = lw[lw > d]

            # lw_oud,nlf_oud,nnf_oud,lf_oud,fit_oud = normalize(wl,flux,line_oud[k+1],line_oud[1][0][1],line_oud[1][1][0],line_oud[1][1][1],line_oud[k+1]-20,line_oud[1][1][1]+20)
            # cond1_oud = (lw>line_oud[k+1])&(lw<line_oud[1][0][1])
            # cond2_oud = (lw>line_oud[1][1][0])&(lw<line_oud[1][1][1])
            # normwave1_oud = lw[cond1_oud]
            # normflux1_oud= lf[cond1_oud]
            # normwave2_oud = lw[cond2_oud]
            # normflux2_oud= lf[cond2_oud]

            # ((lw>c)&(lw<d))]

            # normflux = lf[((lw>line[k+1])&(lw<b)) or ((lw>c)&(lw<d))]
            # ax1.plot(lw,lf)
            f.legend(loc='lower right', bbox_to_anchor=[1.1, 0.04], prop={'size': 10})
            ax3.plot(normwave1, normflux1, c='r', label='continuum data slice')
            ax3.plot(normwave2, normflux2, c='r')
            ax3.plot(normwave3, normflux3, c='b', label='data')
            ax3.plot(normwave4, normflux4, c='b')
            ax3.plot(normwave5, normflux5, c='b')
            # ax1.plot(normwave1_oud,normflux1_oud,c='g')
            # ax1.plot(normwave2_oud, normflux2_oud, c='g')
            # ax1.plot(lw, fit_oud(lw), c='g', label='fitted continuum')
            ax3.plot(lw, fit(lw), c='red', linestyle='--', label='fitted continuum')

            # ax2.plot
            ax4.plot(lw, nlf, c='r')

            ax3.set_title('Peter \n Raw Flux')
            ax4.set_title('Normalized Flux')





            # ax2.plot(lw_oud,nlf_oud,c='g')
            plt.suptitle(line[-2])
            if save == 'on':
                plt.savefig(plot_save_folder + r'\\' + obs +r'\\' + obs+'_' + line[0] + str(int(line[2])) + 'file_'+str(i) + '_normalization.pdf' , format='pdf', dpi=1200 )
            if show == 'on':
                plt.show()
            plt.close()


def normalize_and_plot_from_raw_file(line,filelist):
    v_rad = -18.5
    k=1
    for file in filelist:
        fn = os.path.basename(file)
        data = pf.open(file)
        header = data[0].header
        time_and_date = airmass.timeanddate2(header['DATE-OBS'])
        bjd = airmass.bjd_lapalma_from_date_zet_ori(header['DATE-OBS'])
        print(bjd)
        if 'BVCOR' in header:
            baricentric_correction = float(header['BVCOR'])
            velo_cor = baricentric_correction + v_rad
        naxis1 = header['NAXIS1']
        crval1 = header['CRVAL1']
        cdelt1 = header['CDELT1']
        wl_original = np.exp(np.arange(naxis1) * cdelt1 + crval1 + v_rad / 299792.458)
        flux_original = data[0].data
        lw, lf, nf, _, _ = airmass.normalize(wl_original, flux_original, line[k + 1], line[k + 2], line[k + 3], line[k + 4],
                                             line[k + 1] - 20, line[k + 4] + 20)
        plt.plot(lw,lf,label=time_and_date)
        data.close()
    plt.axvline(line[k + 1])
    plt.axvline(line[k + 2])
    plt.axvline(line[k + 3])
    plt.axvline(line[k + 4])
    plt.ylabel('Relative Flux')
    plt.xlabel(r'Wavelength (Å)')
    plt.show()
    plt.close()
def plot_norm_sil_omar_from_order(linelist,obs='apo',show='on',save = 'off',datafilefolder_apo=None, plot_save_folder = r'D:\peter\Master_Thesis\Datareduction\Plots\normalization\test'):
    if datafilefolder_apo ==None:
        apo_files = open_masterfiles.apo_demetra_orders()
    else:
        apo_files = open_masterfiles.apo_demetra_orders(datafilefolder_apo)
    # print apo_files
    merc_files = open_masterfiles.mercator()
    if obs == 'apo':
        k = 8
        filelist=apo_files
    elif obs == 'mercator':
        k = 8
        filelist=[merc_files[0]]
    else:
        raise(Exception("obs needs to be 'apo' or 'mercator'"))
    i=0
    for file in filelist:
        i+=1
        for line in linelist:
            # line_oud = apo_linelist[i]
            a = line[k][0][0]
            b = line[k][0][1]
            c = line[k][1][0]
            d = line[k][1][1]
            spacing = 200
            start = a-spacing
            stop = d+spacing
            order = airmass.find_order(line[2],file)
            wl = order.wl_rebin
            flux = order.flux_rebin
            print(wl[0])

            lw,nlf,nnf,lf,fit = airmass.normalize(wl,flux,a,b,c,d,start,stop,xtype='velo',linecenter=line[2])
            cond1 = (lw>a)&(lw<b)
            cond2 = (lw>c)&(lw<d)
            normwave1 = lw[cond1]
            normflux1= lf[cond1]
            normwave2 = lw[cond2]
            normflux2= lf[cond2]
            normflux3 = lf[lw<a]
            normwave3 = lw[lw<a]
            normflux4 = lf[(lw>b)&(lw<c)]
            normwave4 = lw[(lw>b)&(lw<c)]
            normflux5 = lf[lw>d]
            normwave5 = lw[lw>d]
            # lw_oud,nlf_oud,nnf_oud,lf_oud,fit_oud = normalize(wl,flux,line_oud[k+1],line_oud[1][0][1],line_oud[1][1][0],line_oud[1][1][1],line_oud[k+1]-20,line_oud[1][1][1]+20)
            # cond1_oud = (lw>line_oud[k+1])&(lw<line_oud[1][0][1])
            # cond2_oud = (lw>line_oud[1][1][0])&(lw<line_oud[1][1][1])
            # normwave1_oud = lw[cond1_oud]
            # normflux1_oud= lf[cond1_oud]
            # normwave2_oud = lw[cond2_oud]
            # normflux2_oud= lf[cond2_oud]
            # ((lw>c)&(lw<d))]
            # normflux = lf[((lw>line[k+1])&(lw<b)) or ((lw>c)&(lw<d))]
            f, (ax1, ax2, ax3, ax4) = plt.subplots(4, sharex=True,layout="constrained")
            # ax1.plot(lw,lf)
            ax1.plot(normwave1,normflux1, c='r', label = 'continuum data slice')
            ax1.plot(normwave2,normflux2, c='r')
            ax1.plot(normwave3,normflux3, c='b', label = 'data')
            ax1.plot(normwave4,normflux4, c='b')
            ax1.plot(normwave5,normflux5, c='b')
            # ax1.plot(normwave1_oud,normflux1_oud,c='g')
            # ax1.plot(normwave2_oud, normflux2_oud, c='g')
            # ax1.plot(lw, fit_oud(lw), c='g', label='fitted continuum')
            ax1.plot(lw,fit(lw),c='red', linestyle='--', label = 'fitted continuum')
            # ax1.legend(prop={'size': 10}, loc='best')
            # ax2.plot
            ax2.plot(lw,nlf,c='r')
            ax1.set_title('Sil en omar \n Raw Flux')
            ax2.set_title('Normalized Flux')
            [a,b,c,d],vsini = airmass.wl_to_velocity([line[3],line[4],line[5],line[6]],line[2])
            # spacing = 200
            # start = a - spacing
            # stop = d + spacing

            wl = order.wl_rebin
            flux = order.flux_rebin
            print(wl[0])

            lw, nlf, nnf, lf, fit = airmass.normalize(wl, flux, a, b, c, d, start, stop, xtype='velo', linecenter=line[2])
            cond1 = (lw > a) & (lw < b)
            cond2 = (lw > c) & (lw < d)
            normwave1 = lw[cond1]
            normflux1 = lf[cond1]
            normwave2 = lw[cond2]
            normflux2 = lf[cond2]
            normflux3 = lf[lw < a]
            normwave3 = lw[lw < a]
            normflux4 = lf[(lw > b) & (lw < c)]
            normwave4 = lw[(lw > b) & (lw < c)]
            normflux5 = lf[lw > d]
            normwave5 = lw[lw > d]

            # lw_oud,nlf_oud,nnf_oud,lf_oud,fit_oud = normalize(wl,flux,line_oud[k+1],line_oud[1][0][1],line_oud[1][1][0],line_oud[1][1][1],line_oud[k+1]-20,line_oud[1][1][1]+20)
            # cond1_oud = (lw>line_oud[k+1])&(lw<line_oud[1][0][1])
            # cond2_oud = (lw>line_oud[1][1][0])&(lw<line_oud[1][1][1])
            # normwave1_oud = lw[cond1_oud]
            # normflux1_oud= lf[cond1_oud]
            # normwave2_oud = lw[cond2_oud]
            # normflux2_oud= lf[cond2_oud]

            # ((lw>c)&(lw<d))]

            # normflux = lf[((lw>line[k+1])&(lw<b)) or ((lw>c)&(lw<d))]
            # ax1.plot(lw,lf)
            f.legend(loc='lower right', bbox_to_anchor=[1.1, 0.04], prop={'size': 10})
            ax3.plot(normwave1, normflux1, c='r', label='continuum data slice')
            ax3.plot(normwave2, normflux2, c='r')
            ax3.plot(normwave3, normflux3, c='b', label='data')
            ax3.plot(normwave4, normflux4, c='b')
            ax3.plot(normwave5, normflux5, c='b')
            # ax1.plot(normwave1_oud,normflux1_oud,c='g')
            # ax1.plot(normwave2_oud, normflux2_oud, c='g')
            # ax1.plot(lw, fit_oud(lw), c='g', label='fitted continuum')
            ax3.plot(lw, fit(lw), c='red', linestyle='--', label='fitted continuum')

            # ax2.plot
            ax4.plot(lw, nlf, c='r')

            ax3.set_title('Peter \n Raw Flux')
            ax4.set_title('Normalized Flux')





            # ax2.plot(lw_oud,nlf_oud,c='g')
            plt.suptitle(line[-2])
            if save == 'on':
                plt.savefig(plot_save_folder + r'\\' + obs +r'\\' + obs+'_' + line[0] + str(int(line[2])) + 'file_'+str(i) + '_normalization.pdf' , format='pdf', dpi=1200 )
            if show == 'on':
                plt.show()
            plt.close()






normalize_and_plot_from_raw_file(Sil_Omar_normline_list_apo_new_short[2],filelist_lapalma)

df_apo = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin02\combined\high_snr\\'
# plot_norm_sil_omar(Sil_Omar_normline_list_apo_new_short,obs='mercator',save='off',show='on')
# plot_norm_sil_omar_from_order(Sil_Omar_normline_list_apo_new_short,obs='apo',datafilefolder_apo=df_apo,save='off',show='on')

