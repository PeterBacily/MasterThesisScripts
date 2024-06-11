from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import glob
import SavitzkyGolay
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
from scipy.optimize import *
from scipy.odr import *
import scipy.stats
from scipy.stats import pearsonr
from PyAstronomy import pyasl
import matplotlib.style
import datareduc
import pickle
import os
import Datafile_class
import open_masterfiles
import warnings
def chunks(l, n):
    # For item i in a range that is a length of l,
    chunklist = []
    for i in range(0, len(l), n):
        # Create an index range for l of n items:
        chunk = l[i:i+n]
        chunklist.append(chunk)
    return chunklist

def my_sin(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*x  + phase)  + offset

def my_sin2(x,  amplitude, phase, offset):
    return amplitude*np.sin(4*np.pi*x  + phase)  + offset


def my_slope(p,x):
    a,b = p
    return a*x+b

# apo_lines  = ['line6562', 'line4861', 'line4713', 'line5875', 'line4541', 'line4685', 'line5411', 'line4471', 'line4921', 'line6678', 'line5592', 'line5801']
# mercator_lines = ['line6562', 'line4861', 'line4340', 'line4026', 'line4471', 'line4713', 'line5875', 'line4541', 'line4685', 'line5411', 'line4921', 'line6678', 'line5592', 'line5801']

def sort_linelist(lines, obs='MERCATOR',apowantedmarks = None):
    if obs =='APO':
        if apowantedmarks == None:
            master_files = open_masterfiles.apo()
        else:
            master_files = open_masterfiles.apo(wantedmarks=apowantedmarks)
    elif obs == 'MERCATOR':
        master_files = open_masterfiles.mercator()
    pfs = []
    for line in lines:
        phases = []
        ews = []
        ew_error = []
        for file in master_files:
            linedata = getattr(file, line)
            lineinfo = linedata.lineinfo
            phases.append(file.phase)
            ews.append(linedata.ew)
            ew_error.append(linedata.ew_error)
        chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats2(np.array(ews), np.array(phases),
                                                                              np.array(ew_error))
        pfs.append((1/probfactor))
    print(lines)
    print (pfs)
    outlist = [x for y,x in sorted(zip(pfs,lines), reverse=True)]
    print(outlist)


def plot_EW(obs='BOTH',apowantedmarks = None, figsavefolder = r'D:\Peter\School\Master Thesis\figures\EWs\verslag\\',custom_lines = None, custom_files = None ):
    # apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']
    apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']

    # mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471', 'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
    mercator_lines = ['line4861', 'line5875', 'line4340', 'line6562', 'line6678', 'line5592', 'line4685', 'line4026', 'line4471',
     'line4921', 'line4713', 'line5801', 'line5411', 'line4541']

    # split_up_lines = list(chunks(mercator_lines,3))
    # print chunks(apo_lines,3)
    # obs = 'APO'
    # obs = 'MERCATOR'
    if obs =='APO':
        if apowantedmarks == None:
            master_files = open_masterfiles.apo()
        else:
            master_files = open_masterfiles.apo(wantedmarks=apowantedmarks)
        lines = apo_lines
    elif obs == 'MERCATOR':
        master_files = open_masterfiles.mercator()
        lines =mercator_lines
    elif obs =='BOTH':
        apo_files = open_masterfiles.apo()
        mercator_files = open_masterfiles.mercator()
        lines = apo_lines
    else:
        master_files = custom_files
        lines = custom_lines
        warnings.warn('observatory is not APO, MERCATOR or BOTH, if this is not intentional look at obs variable')

    if obs != 'BOTH':
        for chunk in chunks(lines,3):
            fig, axs = plt.subplots(nrows=len(chunk), ncols=1,sharex=True,figsize=(5, 8))
            fig.subplots_adjust(top = 0.83)
            size = fig.get_size_inches()
            savename = ''
            for i,ax in enumerate(axs):
                line = chunk[i]
                print(line)
                phases = []
                ews = []
                ew_error = []
                for file in master_files:
                    linedata = getattr(file,line)
                    lineinfo = linedata.lineinfo
                    phases.append(file.phase)
                    ews.append(linedata.ew)
                    ew_error.append(linedata.ew_error)
                chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats2(np.array(ews),np.array(phases),np.array(ew_error))
                print(lineinfo)
                p1 = [0.005, 0.5, 1]
                fit1 = curve_fit(my_sin, phases, ews, p0=p1, sigma=ew_error, absolute_sigma=True)
                fit2 = curve_fit(airmass.flat_line, phases,ews,p0=[1],sigma=ew_error,absolute_sigma=True)
                t = np.linspace(0, 1, 50)
                data_fit = my_sin(t, *fit1[0])
                ax.set_title(lineinfo[-1])
                l1 = ax.errorbar(phases, ews, yerr=ew_error, fmt='o', c='r', label='Data')
                l2, = ax.plot(t, data_fit, c='r', label='Sinusodial fit')
                l3 = ax.axhline(fit2[0][0],linestyle = '--',c='gray', label = 'Weighted Average')
                aldat = np.append(ews, data_fit)
                ax.locator_params(axis='y',nbins=7)
                lowlim,uplim =np.round(np.min(aldat),decimals=2)-0.01,np.round(np.max(aldat),decimals=2)+0.01
                l4 = ax.text(0.75, 0.95*(uplim-lowlim)+lowlim,str('%.2E' % (1./probfactor)),label = 'likelihood ratio')
                ax.set_ylim(uplim,lowlim)
                savename+= lineinfo[0]+line[-4:]+'_'
            savename += obs+'_EW.pdf'
            print(savename)
            plt.figlegend((l1,l2,l3),('Data','Sinusodial fit','Weighted Average'),bbox_to_anchor = [1., 0.92], prop={'size': 10},fancybox=True, framealpha=1)
            plt.suptitle('Normalized Equivalent width \n'+obs, size=16)
            # plt.tight_layout()
            plt.savefig(figsavefolder+savename)
            # plt.show()
            plt.close()

    if obs == 'BOTH':
        for chunk in chunks(lines,3):
            fig, axs = plt.subplots(nrows=len(chunk), ncols=1,sharex=True,figsize=(5, 8))
            fig.subplots_adjust(top = 0.83)
            size = fig.get_size_inches()
            savename = ''
            for i,ax in enumerate(axs):
                line = chunk[i]
                print(line)
                phases_apo = []
                phases_mercator = []
                phases = []
                ews_apo = []
                ews_mercator = []
                ews = []
                ew_error_apo = []
                ew_error_mercator = []
                ew_error = []
                for file in apo_files:
                    linedata = getattr(file,line)
                    lineinfo = linedata.lineinfo
                    phases.append(file.phase)
                    phases_apo.append(file.phase)
                    ews.append(linedata.ew)
                    ews_apo.append(linedata.ew)
                    ew_error.append(linedata.ew_error)
                    ew_error_apo.append(linedata.ew_error)
                for file in mercator_files:
                    linedata = getattr(file,line)
                    lineinfo = linedata.lineinfo
                    phases.append(file.phase)
                    phases_mercator.append(file.phase)
                    ews.append(linedata.ew)
                    ews_mercator.append(linedata.ew)
                    ew_error.append(linedata.ew_error)
                    ew_error_mercator.append(linedata.ew_error)
                norm_EW_apo = np.array(ews_apo)/np.average(ews_apo)
                norm_error_apo =np.array(ew_error_apo)/np.average(ews_apo)
                norm_EW_mercator = np.array(ews_mercator)/np.average(ews_mercator)
                norm_error_mercator = np.array(ew_error_mercator) / np.average(ews_mercator)
                norm_ews = np.append(norm_EW_apo,norm_EW_mercator)
                norm_ew_errors = np.append(norm_error_apo,norm_error_mercator)
                norm_phases = np.append(phases_apo,phases_mercator)
                chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats2(norm_ews,norm_phases,norm_ew_errors)
                # p1 = [0.005, 0.5, 1]
                p1 = [1,1,1]
                # print norm_ews
                fit1 = curve_fit(my_sin, phases, norm_ews, p0=p1, sigma=norm_ew_errors, absolute_sigma=True)
                fit2 = curve_fit(airmass.flat_line, phases,norm_ews,p0=[1],sigma=norm_ew_errors,absolute_sigma=True)
                fit_sin_apo = curve_fit(my_sin, phases_apo, norm_EW_apo, p0=p1, sigma=norm_error_apo, absolute_sigma=True)
                fit_sin_mercator = curve_fit(my_sin, phases_mercator, norm_EW_mercator, p0=p1, sigma=norm_error_mercator, absolute_sigma=True)
                t = np.linspace(0, 1, 50)
                data_fit = my_sin(t, *fit1[0])
                data_fit_apo = my_sin(t, *fit_sin_apo[0])
                data_fit_mercator = my_sin(t, *fit_sin_mercator[0])
                ax.set_title(lineinfo[-1])
                l1_apo = ax.errorbar(phases_apo, norm_EW_apo, yerr =norm_error_apo, fmt='o', c='r', label='APO Data')
                l1_mercator = ax.errorbar(phases_mercator, norm_EW_mercator, yerr=norm_error_mercator, fmt='o', c='b', label='MERCAROR Data')
                l2a, = ax.plot(t, data_fit_apo, c='r', label='Sinusodial fit APO')
                l2b = ax.plot(t, data_fit_mercator, c='b', label='Sinusodial fit mercator')
                l3 = ax.axhline(fit2[0][0],linestyle = '--',c='gray', label = 'Weighted Average')
                # l5 = ax.errorbar(phases, norm_ews, yerr=norm_ew_errors, fmt='o', c='g', label='combined_Data')
                aldat = np.append(norm_ews, data_fit)
                ax.locator_params(axis='y',nbins=7)
                lowlim,uplim =np.round(np.min(aldat),decimals=2)-0.01,np.round(np.max(aldat),decimals=2)+0.01
                l4 = ax.text(0.75, 0.95*(uplim-lowlim)+lowlim,str('%.2E' % (1./probfactor)),label = 'likelihood ratio')
                # ax.set_ylim(uplim,lowlim)
                savename+= lineinfo[0]+line[-4:]+'_'
            savename += obs+'_EW.pdf'
            print(savename)
            plt.figlegend((l1_apo,l1_mercator,l2a,l3),('APO Data','MERCATOR Data','Sinusodial fit','Weighted Average'),bbox_to_anchor = [1., 0.92], prop={'size': 10},fancybox=True, framealpha=1)
            plt.suptitle('Normalized Equivalent width \n'+obs, size=16)
            # plt.tight_layout()
            # plt.savefig(figsavefolder+savename)
            plt.show()
            plt.close()

def plot_TVS(obs='MERCATOR',save='off',show='on',plot_save_folder='',oneline='off',sg='off', apowantedmarks = None,customlines = None,customfiles = None):
    apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
                      'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']

    if obs == 'APO':
        if apowantedmarks == None:
            master_files = open_masterfiles.apo()
        else:
            master_files = open_masterfiles.apo(wantedmarks=apowantedmarks)
        lines = apo_lines
    elif obs == 'MERCATOR':
        master_files = open_masterfiles.mercator()
        lines = mercator_lines
    else:
        master_files = customfiles
        lines = customlines
        if master_files==None or lines == None:
            raise AttributeError('Either select observatory correctly or give custom lines and files')
    vsini = 127
    for line in lines:
        wl,TVS,v,n = airmass.TVS_masterfiles(master_files,line)
        lineinfo = getattr(master_files[0], line).lineinfo
        print(lineinfo)
        sgn = 91  # window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        # print v[sgn]-v[0]
        p = scipy.stats.distributions.chi2.ppf(0.99, n - 1) / (n - 1)
        vs, lfs = airmass.overplot_masterfiles(master_files, line)
        f, (ax1, ax2) = plt.subplots(2, sharex=True)
        spec2 =[]
        for i, spec in enumerate(lfs):
            ax1.plot(vs[i], spec, linewidth=1.0)
            spec2.append( spec[(v > -300) & (v < 300)])

        plt.suptitle(lineinfo[-1], size=18)
        ax1.set_title('Spectra ' +obs)
        ax2.set_title('TVS')
        # ax1.legend()
        # ax1.set_xlim([-600,600])
        spec3 = np.array(spec2)
        mini = np.floor(20 * 0.99 * np.amin(spec3)) / 20
        maxi = np.ceil(20 * 1.01 * np.amax(spec3)) / 20
        print(maxi, np.amax(spec3))
        ax1.set_ylim([mini, maxi])
        ax1.axvline(vsini, color='k', linestyle=':', linewidth=1)
        ax1.axvline(-vsini, color='k', linestyle=':', linewidth=1)
        # if line[2]==5875.621:
        #     TVS2 = np.array(TVS)*1.4
        ax2.plot(v, TVS, color='b')
        if sg == 'on':
            ax2.plot(v, TVS_smoothed, color='r', linestyle='dashed')
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
        TVS2 = np.array(TVS)[(v > -200) & (v < 200)]
        # print TVS2
        # print np.amax(TVS2)
        maxi2 = np.ceil(np.amax(TVS2))
        ax2.set_ylim([0, maxi2])
        ax2.set_xlabel('V (km/s)')
        ax1.set_ylabel('Normlized Flux')
        ax2.set_ylabel(r'$\sigma_{obs}$' + r' \ ' + r'$\sigma_{exp}$', size=16)
        ax2.set_xlim([-600, 600])
        print(line[4:])
        if save == 'on':
            plt.savefig(plot_save_folder + r'\\'+obs + lineinfo[0] + line[4:] + '_TVS.pdf', format='pdf',
                        dpi=1200)
        if show == 'on':
            plt.show()
        plt.close()
    # vlim= 300
    # plotv = v[(v>-vlim)& (v<vlim)]
    # plotTVS = TVS[(v>-vlim)& (v<vlim)]
    # plt.plot(plotv,plotTVS)
    # plt.show()
    # plt.close()



def plot_quotient(obs='MERCATOR',save='off',show='on',plot_save_folder='',overplot = 'off',oneline='off',sg='off', apowantedmarks = None,customlines = None,customfiles = None):
    # print datafile_folder
    savepdf =  matplotlib.backends.backend_pdf.PdfPages(r'D:\Peter\School\Master Thesis\figures\quotient_' + obs + '_overplot.pdf')

    apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
                      'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
    vsini = 127
    if obs == 'APO':
        if apowantedmarks == None:
            master_files = open_masterfiles.apo()
        else:
            master_files = open_masterfiles.apo(wantedmarks=apowantedmarks)
        lines = apo_lines
    elif obs == 'MERCATOR':
        master_files = open_masterfiles.mercator()
        lines = mercator_lines
    else:
        master_files = customfiles
        lines = customlines
        if master_files==None or lines == None:
            raise AttributeError('Either select observatory correctly or give custom lines and files')
    # filelist = glob.glob(datafile_folder+'\*.fit')
    # fileinfo_header = dict_eshel['header']
    # print fileinfo_header
    # print filelist
    for i in range(len(master_files) - 1):
        file1 = master_files[i]
        file2 = master_files[i + 1]
        t1 = file1.HJD
        t2 = file2.HJD
        dt = t2 - t1
        print(dt)
        phi1 = file1.phase
        phi2 = file2.phase
        # print phi1,phi2
        dphi=(phi2 - phi1 + 1) % 1
        # file1_id = file1_info[1]+str(file1_info[0])
        # file2_id = file2_info[1] + str(file2_info[0])
        # fi=[file1_info,file2_info]
        for line in lines:
            lineobject1 = getattr(file1, line)
            lineobject2 = getattr(file2, line)
            lineinfo1 = lineobject1.lineinfo
            lineinfo2 = lineobject2.lineinfo
            v1 = np.array(lineobject1.v_cor)
            v2 = np.array(lineobject2.v_cor)
            flux1 =np.array(lineobject1.flux)
            flux2 = np.array(lineobject2.flux)
            qf = np.array(airmass.quotient(flux1,flux2))
            # swl = line[3]-40
            # ewl = line[6]+40
            # lw,v,qf =airmass.quotient_eShel(file1,file2,line,swl,ewl, v_rad=18.5)

            sgn = 101  # window size for SavitzkyGolay (must be odd integer)
            qf_smoothed = SavitzkyGolay.savitzky_golay(qf, sgn, 4)
            # print v[sgn] - v[0]
            # p = chi2.ppf(0.99, n-1)/(n-1)
# Maak textbox met info en legenda en shit.
            if overplot == 'on':
                # vs,lfs = airmass.overplot_masterfiles([file1,file2],line)
                f, (ax1, ax2) = plt.subplots(2, sharex=True)
                lb1 = obs + ' '+ str(file1.i)
                lb2 = obs + ' '+ str(file2.i)
                ax1.plot(v1, flux1, linewidth=1.0, label=lb1)
                ax1.plot(v2, flux2, linewidth=1.0, label=lb2)

                # for i,spec in enumerate(lfs):
                #     lb = fi[i][1]+str(fi[i][0])
                #     ax1.plot(vs[i],spec,linewidth=1.0,label = lb)
                # ax1.set_title(line[7])
                # ax1.legend()
                # ax1.set_xlim([-600,600])
                spec2 = np.append(flux1[(v1>-300)& (v1<300)],flux2[(v1>-300)& (v1<300)])
                mini = np.floor(10*0.99*np.amin(spec2))/10
                maxi = np.ceil(10*1.01*np.amax(spec2))/10
                print(mini, maxi)
                # ax1.set_ylim([mini,maxi])
                ax1.axvline(vsini, color='k', linestyle='dashed', linewidth=1)
                ax1.axvline(-vsini, color='k', linestyle='dashed', linewidth=1)
                ax1.set_ylabel('Normalized Flux', size=14)
                ax1.legend(loc=4,fancybox=True,framealpha=0.5)
                # if line[2]==5875.621:
                #     TVS2 = np.array(TVS)*1.4
            else:
                f, (ax2) = plt.subplots(1, sharex=True)
            plt.suptitle('Quotient spectrum ' + lineinfo1[-1] +'\n' + obs + str(file2.i)+'/'+ obs +str(file1.i), size=14)
            ax2.plot(v1, qf, color='b')
            informationtext = '$\Delta$t = '+"{:.3f}".format(dt) +'d'+ '\n'+'$\Delta\phi$ = ' +"{:.3f}".format(dphi)
            ax2.text(0.814,0.951,informationtext, transform=ax2.transAxes, verticalalignment='top', bbox = dict(boxstyle='round', facecolor='white', alpha=0.5))
            if sg=='on':
                ax2.plot(v1,qf_smoothed,color='r',linestyle='dashed')
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
            # TVS2 = np.array(qf)[(v>-200)& (v<200)]
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
                plt.savefig(plot_save_folder +r'\\' + obs+'_' + lineinfo1[0]+line[4:] +str(file2.i)+'_' +str(file1.i)+'_Quotient_Overplot.pdf',format='pdf', dpi=1200)
                savepdf.savefig(f)
            elif save =='on' and overplot == 'off':
                plt.savefig(plot_save_folder +r'\\' + obs+'_' + lineinfo1[0]+line[4:] +str(file2.i)+'_' +str(file1.i)+ '_Quotient_bare.pdf',format='pdf', dpi=1200)
            if show =='on':
                plt.show()
            plt.close()
    savepdf.close()

def plot_TVS_together(save='off',show='on',plot_save_folder='',oneline='off',sg='on', siglvlline=0.01,apowantedmarks = None,customlines = None,customfiles = None):
    apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
                      'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']


    apo_master_files = open_masterfiles.apo()

    merc_master_files = open_masterfiles.mercator()

    vsini = 127

    # for line in [apo_lines[1]]:
    for line in apo_lines:
        f, axarr = plt.subplots(2, 2, sharex=True, figsize = (8,6))
        for i in [0, 1]:
            if i == 0:
                master_files = apo_master_files
                obs = 'APO'
                k=2
            elif i == 1:
                master_files = merc_master_files
                obs = 'MERCATOR'
                k=1
            wl, TVS, v, n = airmass.TVS_masterfiles(master_files, line)
            lineinfo = getattr(master_files[0], line).lineinfo

            normwl_edges = [lineinfo[k+1],lineinfo[k+2],lineinfo[k+3],lineinfo[k+4]]
            norm_v_edges =airmass.wl_to_velocity(normwl_edges,lineinfo[k])[0]
            sgn = 91  # window size for SavitzkyGolay (must be odd integer)
            TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
            # print v[sgn]-v[0]
            # p = chi2.ppf(0.99, n - 1) / (n - 1)
            vs, lfs = airmass.overplot_masterfiles(master_files, line)
            # f, (axarr[0, i], axarr[0, i]) = plt.subplots(2, sharex=True)
            spec2 = []
            for j, spec in enumerate(lfs):
                axarr[0, i].plot(vs[j], spec, linewidth=1.0)
                spec2.append(spec[(v > -300) & (v < 300)])

            plt.suptitle(lineinfo[-1], size=24)
            axarr[0, i].set_title(r'$\bf{'+obs+r'}$' +'\n Spectra ')
            axarr[1, i].set_title('TVS ')
            # axarr[0, i].legend()
            # axarr[0, i].set_xlim([-600,600])
            spec3 = np.array(spec2)
            mini = np.floor(20 * 0.99 * np.amin(spec3)) / 20
            maxi = np.ceil(20 * 1.01 * np.amax(spec3)) / 20
            # print maxi, np.amax(spec3)
            axarr[0, i].set_ylim([mini, maxi])
            axarr[0, i].axvline(vsini, color='k', linestyle=':', linewidth=1)
            axarr[0, i].axvline(-vsini, color='k', linestyle=':', linewidth=1)
            for edge in norm_v_edges:
                axarr[0, i].axvline(edge, color='green', linestyle=':', linewidth=1)
            # if line[2]==5875.621:
            #     TVS2 = np.array(TVS)*1.4
            axarr[1, i].plot(v, TVS, color='b')
            if sg == 'on':
                axarr[1, i].plot(v, TVS_smoothed, color='r', linestyle='dashed')
            if oneline == 'on':
                axarr[1, i].axhline(y=1, color='gray', linestyle='--')
            if isinstance(siglvlline, float):
                Nfiles = len(master_files)
                p = siglvlline
                siglvl = airmass.TVS_significance_level(Nfiles, p)
                axarr[1, i].axhline(y=siglvl, color='red', linestyle='--')
            # else:
            #     axarr[1, i].plot(v,TVS)
            axarr[1, i].axvline(vsini, color='k', linestyle=':', linewidth=1)
            axarr[1, i].axvline(-vsini, color='k', linestyle=':', linewidth=1)
            # axarr[1, i].plot([v[0],v[-1]],[1,1],linestyle='dashed', linewidth=1,c='g')
            # axarr[1, i].plot([v[0],v[-1]],[p,p], linewidth=1,c='g')
            # print len(v)
            # print len(TVS)
            # print v
            TVS2 = np.array(TVS)[(v > -200) & (v < 200)]
            # print TVS2
            # print np.amax(TVS2)
            maxi2 = np.ceil(np.amax(TVS2))
            axarr[1, i].set_ylim([0, maxi2])
            # axarr[1,1].set_ylim([0,3])
            axarr[1, i].set_xlabel('V (km/s)',size = 14)
            axarr[0, 0].set_ylabel('Normlized Flux', size = 14)
            axarr[1, 0].set_ylabel(r'$\sigma_{obs}$' + r' \ ' + r'$\sigma_{exp}$', size=20)
            # axarr[1, i].set_xlim([-600, 600])
            # print line[4:]
            # plt.tight_layout()
            plt.subplots_adjust(left=0.08, bottom=None, right=0.98, top=None,
                            wspace=None, hspace=0.15)
        if save == 'on':
            plt.savefig(plot_save_folder + r'\\TVS_' +  lineinfo[0] + line[4:] + '_TVS.pdf', format='pdf',
                            dpi=1200)
        if show == 'on':
            plt.show()
        plt.close()

def plot_EW_together(save='off',show='on',plot_save_folder='',oneline='off',sg='off', apowantedmarks = None,customlines = None,customfiles = None, bad = 'off',fx2='off'):
    # fig, axs = plt.subplots(nrows=len(chunk), ncols=1, sharex=True, figsize=(5, 8))
    # fig.subplots_adjust(top=0.83)
    # size = fig.get_size_inches()
    apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    apo_master_files = open_masterfiles.apo()

    merc_master_files = open_masterfiles.mercator()

    vsini = 127

    for line in apo_lines:
    # for line in [apo_lines[1]]:
        f, axarr = plt.subplots(2, sharex=True, figsize=(9, 6))
        savename = ''
        for i in [0,1]:
            if i == 0:
                master_files = apo_master_files
                obs = 'APO'
            elif i == 1:
                master_files = merc_master_files
                obs = 'MERCATOR'
            print(line)
            phases = []
            ews = []
            ew_error = []
            for file in master_files:
                linedata = getattr(file, line)
                lineinfo = linedata.lineinfo
                phases.append(file.phase)
                ews.append(linedata.ew)
                ew_error.append(linedata.ew_error)
            chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats3(np.array(ews), np.array(phases),
                                                                                  np.array(ew_error),fx2=fx2)
            if bad =='on':
                chisq2, red_chisq2, AIC_12, AIC_22, AIC_32, probfactor2 = airmass.EW_stats2(np.array(ews), np.array(phases),
                                                                                      np.array(ew_error))
                add = '\n'+ str('LR: ' + '%.2E' % (1. / probfactor2))
            else:
                add = ''
            print(lineinfo)
            p1 = [0.005, 0.5, 1]
            t = np.linspace(0, 1, 50)
            if fx2 == 'on':
                fit1 = curve_fit(my_sin2, phases, ews, p0=p1, sigma=ew_error, absolute_sigma=True)
                data_fit = my_sin2(t, *fit1[0])
            else:
                fit1 = curve_fit(my_sin, phases, ews, p0=p1, sigma=ew_error, absolute_sigma=True)
                data_fit = my_sin(t, *fit1[0])
            fit2 = curve_fit(airmass.flat_line, phases, ews, p0=[1], sigma=ew_error, absolute_sigma=True)
            axarr[i].set_title(obs, size=14, weight='light')
            l1 = axarr[i].errorbar(phases, ews, yerr=ew_error, fmt='o', c='r', label='Data')
            l2, = axarr[i].plot(t, data_fit, c='r', label='Sinusodial fit')
            l3 = axarr[i].axhline(fit2[0][0], linestyle='--', c='gray', label='Weighted Average')
            aldat = np.append(ews, data_fit)
            # ax.locator_params(axis='y', nbins=7)
            lowlim, uplim = np.round(np.min(aldat), decimals=2) - 0.01, np.round(np.max(aldat), decimals=2) + 0.01
            l4 = axarr[i].text(0.78, (0.01 * (uplim - lowlim)) + lowlim, str('LR: ' + '%.2E' % (1. / probfactor))+add, label='likelihood ratio',size=14)
            # l4 = axarr[i].text( str('%.2E' % (1. / probfactor)),
            #                    label='likelihood ratio')
            axarr[i].set_ylim( lowlim,uplim)
        savename += lineinfo[0] + line[-4:] + '_'
        savename += '_EW.pdf'
        print(savename)
        plt.figlegend((l1, l2, l3), ('Data', 'Sinusodial fit', 'Weighted Average'), bbox_to_anchor=[0.87, 0.85],
                      prop={'size': 11}, fancybox=True, framealpha=0.8)
        plt.suptitle('Normalized Equivalent Width \n' + lineinfo[-1], size=18)
        plt.subplots_adjust(top=0.85)
        plt.xlabel('Phase (P = 6.83d)',fontsize = 13)
        # plt.tight_layout()
        if save =='on':
            plt.savefig(plot_save_folder + savename)
        if show == 'on':
            plt.show()
        plt.close()

def EW_peak_table():
    # fig, axs = plt.subplots(nrows=len(chunk), ncols=1, sharex=True, figsize=(5, 8))
    # fig.subplots_adjust(top=0.83)
    # size = fig.get_size_inches()
    # apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
    #              'line4861', 'line4921', 'line6678', 'line4471']
    apo_lines = ['line6562', 'line4861', 'line4471', 'line4713', 'line4921', 'line5875', 'line6678', 'line4541',
              'line4685', 'line5411', 'line5801', 'line5592']
    apo_master_files = open_masterfiles.apo()

    merc_master_files = open_masterfiles.mercator()

    vsini = 127

    for line in apo_lines:
        # for line in [apo_lines[1]]:
        for i in [0, 1]:
            if i == 0:
                master_files = apo_master_files
                obs = 'APO'
            elif i == 1:
                master_files = merc_master_files
                obs = 'MERCATOR'
            phases = []
            ews = []
            ew_error = []
            for file in master_files:
                linedata = getattr(file, line)
                lineinfo = linedata.lineinfo
                phases.append(file.phase)
                ews.append(linedata.ew)
                ew_error.append(linedata.ew_error)
            chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats3(np.array(ews), np.array(phases),np.array(ew_error))
            p1 = [0.005, 0.5, 1]
            popt, pcov = curve_fit(my_sin, phases, ews, p0=p1, sigma=ew_error, absolute_sigma=True)
            fit2 = curve_fit(airmass.flat_line, phases, ews, p0=[1], sigma=ew_error, absolute_sigma=True)

            # fm = lambda x: -my_sin(x, *popt)
            # r = minimize_scalar(fm, bounds=(0, 1))
            t = np.linspace(0, 1, 50000)
            y = my_sin(t, *popt)
            imax = np.argmax(y)
            # data_fit = my_sin(t, *popt)
            # aldat = np.append(ews, data_fit)
            pf = str('%.2E' % (1. / probfactor))
            if 1/probfactor> 100:
                det = 'Yes'
            elif 1/probfactor>20:
                det = 'Pos'
            else:
                det = 'No'
            if 1/probfactor>20:
                print(obs, lineinfo[-1], '   &',pf, '   &', '%.2f' % t[imax],r'\\')



def EW_fit_table():
    # fig, axs = plt.subplots(nrows=len(chunk), ncols=1, sharex=True, figsize=(5, 8))
    # fig.subplots_adjust(top=0.83)
    # size = fig.get_size_inches()
    # apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
    #              'line4861', 'line4921', 'line6678', 'line4471']
    apo_lines = ['line6562', 'line4861', 'line4471', 'line4713', 'line4921', 'line5875', 'line6678', 'line4541',
              'line4685', 'line5411', 'line5801', 'line5592']
    apo_master_files = open_masterfiles.apo()

    merc_master_files = open_masterfiles.mercator()

    vsini = 127

    for line in apo_lines:
        # for line in [apo_lines[1]]:
        for i in [0, 1]:
            if i == 0:
                master_files = apo_master_files
                obs = 'APO'
            elif i == 1:
                master_files = merc_master_files
                obs = 'MERCATOR'
            phases = []
            ews = []
            ew_error = []
            for file in master_files:
                linedata = getattr(file, line)
                lineinfo = linedata.lineinfo
                phases.append(file.phase)
                ews.append(linedata.ew)
                ew_error.append(linedata.ew_error)
            chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats3(np.array(ews), np.array(phases),np.array(ew_error))
            if i == 0:
                print(lineinfo[-1], r'  & \ref{fig:EW_' + lineinfo[-1][-4:] +r'} &'),
            p1 = [0.005, 0.5, 1]
            fit1 = curve_fit(my_sin, phases, ews, p0=p1, sigma=ew_error, absolute_sigma=True)
            fit2 = curve_fit(airmass.flat_line, phases, ews, p0=[1], sigma=ew_error, absolute_sigma=True)
            t = np.linspace(0, 1, 50)
            data_fit = my_sin(t, *fit1[0])
            aldat = np.append(ews, data_fit)
            pf = str('%.2E' % (1. / probfactor))
            if 1/probfactor> 100:
                det = 'Yes'
            elif 1/probfactor>20:
                det = 'Pos'
            else:
                det = 'No'

            print(pf, '   &', det,'   &',)
        print(r'\\')

def get_mfs(obs):
    if obs =='APO':
        mfs = open_masterfiles.apo()
    elif obs =='MERC' or obs =='MERCATOR':
        mfs =open_masterfiles.mercator()
    return mfs

def EW_EW_plots(line1,obs,line2,save = 'off',plot_save_folder = '',show='off',invert='off'):
    apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    master_files1 = get_mfs(obs)
    master_files2 = get_mfs(obs)

    phases1 = []
    ews1 = []
    ew_error1 = []
    if invert=='on':
        # print 'onononon'
        line_1 = line2
        line_2 = line1
    else:
        line_1 = line1
        line_2 = line2

    for file in master_files1:
        linedata1 = getattr(file, line_1)
        lineinfo1 = linedata1.lineinfo
        phases1.append(file.phase)
        ews1.append(linedata1.ew)
        ew_error1.append(linedata1.ew_error)
    phases2 = []
    ews2 = []
    ew_error2 = []
    for file in master_files2:
        linedata2 = getattr(file, line_2)
        lineinfo2 = linedata2.lineinfo
        phases2.append(file.phase)
        ews2.append(linedata2.ew)
        ew_error2.append(linedata2.ew_error)
    pearsonr = scipy.stats.pearsonr(ews1,ews2)
    plt.errorbar(x= ews1, y =ews2,xerr = ew_error1,yerr = ew_error2, linestyle= 'None',fmt='o')
    # fitparams = np.polyfit(x= ews1, y =ews2,deg=1)
    axes = plt.gca()
    xmin,xmax = axes.get_xlim()

    # y=x*fitparams[0]+fitparams[1]
    slope_model = Model(my_slope)
    data = RealData(np.array(ews1), np.array(ews2), sx=np.array(ew_error1), sy=np.array(ew_error2))
    odr = ODR(data, slope_model, beta0=[0.5, 0.5])
    out = odr.run()
    x = np.linspace(xmin, xmax,25)
    y=my_slope(out.beta, x)
    plt.plot(x,y,c='r')
    plt.title('Equivalent Width correlation of '+lineinfo1[-1] +' '+ ' and\n'+lineinfo2[-1]+' (' +obs +')' ,size=14, weight='light')
    # print pearsonr
    savename =  'EW_CORRELATION_'+'%.4s' % (obs,)+'_'+lineinfo1[0]+line1[-4:]+ '_'+lineinfo2[0]+line2[-4:]+'.pdf'
    plt.text(0.7, 0.01, 'Pearson R: ' + '%.2f' % pearsonr[0] +'\np-value:' +
             '%.2E' % pearsonr[1]    , size=14, transform=axes.transAxes)
    plt.xlabel('EW '+lineinfo1[-1] , fontsize = 13)
    plt.ylabel('EW '+lineinfo2[-1] , fontsize = 13)
    if save == 'on':
        plt.savefig(plot_save_folder + savename)
    if show == 'on':
        plt.show()
    plt.close()


def EW_CORRELATION(line1, obs, line2):
    apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    master_files1 = get_mfs(obs)
    master_files2 = get_mfs(obs)

    phases1 = []
    ews1 = []
    ew_error1 = []
    for file in master_files1:
        linedata1 = getattr(file, line1)
        lineinfo1 = linedata1.lineinfo
        phases1.append(file.phase)
        ews1.append(linedata1.ew)
        ew_error1.append(linedata1.ew_error)
    phases2 = []
    ews2 = []
    ew_error2 = []
    for file in master_files2:
        linedata2 = getattr(file, line2)
        lineinfo2 = linedata2.lineinfo
        phases2.append(file.phase)
        ews2.append(linedata2.ew)
        ew_error2.append(linedata2.ew_error)
    pearsonr = scipy.stats.pearsonr(ews1, ews2)
    return [obs, lineinfo1, lineinfo2, pearsonr]


def loop_EW_EW(obs = 'APO', save = 'on',show = 'off',plot_save_folder = r'D:\Peter\School\Master Thesis\figures\EWs\EW-EW\goodfit\\'):
    linelist  = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    j=1
    for i in range(len(linelist)-1):
        line1 = linelist[i]
        for k in range(len(linelist)-1-i):
            line2 = linelist[k+1+i]
            EW_EW_plots(line1,obs,line2,save=save,plot_save_folder = plot_save_folder,show=show, invert='off')
            print(obs, j)
            j+=1
def loop_EW_COR():
    linelist  = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    # linelist  = ['line6562', 'line4713', 'line5801']
    print('OBS  line 1  line 2  pearsonr    p-value')
    obses = ['APO','MERC']
    ew_cor_list = []
    for obs in obses:
        for i in range(len(linelist)-1):
            line1 = linelist[i]
            for k in range(len(linelist)-1-i):
                print(i,k)
                line2 = linelist[k+1+i]
                ewcor = EW_CORRELATION(line1,obs,line2)
                ew_cor_list.append(ewcor)
    with open(r"D:\Peter\School\Master Thesis\Data\EW\EW_correlation.txt", "wb") as fp:
        pickle.dump(ew_cor_list, fp)

def open_EW_cor_list(path = r"D:\Peter\School\Master Thesis\Data\EW\EW_correlation.txt"):
    with open(path, "rb") as fp:  # Unpickling
        ewcorlist = pickle.load(fp)
    return ewcorlist

def sorted_ew_correlation_lists():
    cl = open_EW_cor_list()
    cllist_apo = []
    cllist_merc = []
    for item in cl:
        if item[0]=='APO':
            cllist_apo.append(item)
        elif item[0]=='MERC':
            cllist_merc.append(item)
        else:
            print('error')

    sorted_apo = sorted(cllist_apo, key=lambda x: x[3][0], reverse=True)
    sorted_merc = sorted(cllist_merc, key=lambda x: x[3][0], reverse=True)
    return sorted_apo,sorted_merc

def compare_average(save = 'off',show = 'on',plot_save_folder = r'D:\peter\Master_Thesis\Datareduction\Plots\Average_comparison\\'):
    apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
                      'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']

    apo_master_files = open_masterfiles.apo()
    mercator_master_files = open_masterfiles.mercator()
    for line in apo_lines:
        lineinfo = getattr(apo_master_files[0], line).lineinfo
        wl, avg, v, n = airmass.average_masterfiles(apo_master_files, line)
        wl2, avg2, v2, n2 = airmass.average_masterfiles(mercator_master_files, line)
        # vsini = 127
        plt.title('Average of normalized spectra of the \n ' + lineinfo[-1] + 'line',
                  size=14, weight='light')

        # savename = 'EW_CORRELATION_' + '%.4s' % (obs,) + '_' + lineinfo1[0] + line1[-4:] + '_' + lineinfo2[0] + line2[-4:] + '.pdf'
        plt.plot(v, avg,label='APO')
        plt.plot(v2, avg2,label = 'Mercator')
        plt.xlim([-600, 600])
        plt.ylabel('Normalized flux')
        plt.xlabel('velocity (km/s)')
        plt.legend()
        if save == 'on':
            plt.savefig(plot_save_folder + r'\\Average_' + lineinfo[0] + line[4:] + '_line.pdf', format='pdf',
                        dpi=1200)
        if show == 'on':
            plt.show()
        plt.close()

# a = sorted_ew_correlation_lists()
# i=1
# for line in a[1][:15]:
#     print i, line[1][-1], '       ',line[2][-1]
#     print '---------------------------------'
#     i+=1

# apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']
# mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471', 'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
# sort_linelist(apo_lines, obs='APO',apowantedmarks = None)
#2

# plot_quotient(obs='APO',overplot='on', oneline='on', show='off', save='on', plot_save_folder=r'D:\Peter\School\Master Thesis\figures\Quotients\apo\overplot')
# plot_quotient(obs='MERCATOR',overplot='on', oneline='on', show='off', save='on', plot_save_folder=r'D:\Peter\School\Master Thesis\figures\Quotients\mercator\overplot')
# plot_TVS(obs = 'APO',oneline='on',plot_save_folder=r'D:\Peter\School\Master Thesis\figures\TVS\eShel\from_masterfiles',save='on',show='off')
# plot_TVS(obs = 'MERCATOR',oneline='on',plot_save_folder=r'D:\Peter\School\Master Thesis\figures\TVS\LaPalma\from_masterfiles',save='on',show='off')
# plot_EW(obs='MERCATOR',apowantedmarks = None)
# plot_EW(obs='APO',apowantedmarks = None)
plot_TVS_together(plot_save_folder= r'D:\peter\Master_Thesis\Datareduction\Plots\TVS\combined\no_sg', save = 'off', show='on',sg='off', oneline='on')
# plot_EW_together(plot_save_folder= r'D:\Peter\School\Master Thesis\figures\EWs\verslag\half_period\EW_', save = 'off', show='off', oneline='on', bad='off',fx2='on')
# compare_average(save = 'on',show = 'off')


# EW_EW_plots('line6562','APO','line4861',save = 'off',show = 'on', plot_save_folder = r'D:\Peter\School\Master Thesis\figures\EWs\EW-EW\inverted\\' )
# loop_EW_EW(obs= 'MERCATOR', save='on',show='off')
# loop_EW_EW(obs= 'APO', save='on',show='off')
# loop_EW_COR()
# loop_EW_COR(obs='MERC')
# EW_peak_table()

# ew_cor_list = open_EW_cor_list()
# linewavelengths = [6563,4713,5801,4685,5875,4541,5411,5592,4861,4921,6678,4471]
# # for lw in linewavelengths:
# #     a = next(obj for obj in ew_cor_list if int(obj[2][-1][-4:])==lw)
# #     print a[1][-1]
# linenamelist = [ew_cor_list[0][1][-1]]
# linelist = [ew_cor_list[0][1]]
# for item in ew_cor_list[:11]:
#     linelist.append(item[2])
# #     if int(item[1][-1][-4:]) in linewavelengths and int(item[2][-1][-4:]) in linewavelengths:
#     linenamelist.append(item[2][-1])
# print 'Line',
# for line in linelist:
#     print ' & ' + line[-1][:-4],
# print r' \\'
# print ' ',
# for line in linelist:
#     print ' & ' + line[-1][-4:],
# print r' \\'
# obs = 'MERC'
# for line1 in linenamelist:
#     yline = line1[-4:]
#     # print line1,
#     halves = ['top','bottom']
#     for half in halves:
#         if half == 'top':
#             print line1[:-4],
#
#         else:
#             print line1[-4:],
#         for line2 in linenamelist:
#             xline = line2[-4:]
#             xylines = [xline,yline]
#             # print xylines
#             if xline == yline:
#                 print '&--',
#             for item in ew_cor_list:
#                 if item[1][-1][-4:] in xylines and item[2][-1][-4:] in xylines and item[0]==obs:
#                     if half =='top':
#                         print  '&','%.2f' % item[3][0],
#                     else:
#                         print  '&','%.1e' % item[3][1],
#                     # print item[1][-1], item[2][-1],
#         print r'\\'
#     print r'\hline'
#     else:
#         print 'missing'