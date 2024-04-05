from __future__ import division
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import glob
import SavitzkyGolay
import pyfits as pf
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
    print lines
    print pfs
    outlist = [x for y,x in sorted(zip(pfs,lines), reverse=True)]
    print outlist


def plot_EW(obs='MERCATOR',apowantedmarks = None):
    # apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']
    apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']

    # mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471', 'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
    mercator_lines = ['line4861', 'line5875', 'line4340', 'line6562', 'line6678', 'line5592', 'line4685', 'line4026', 'line4471',
     'line4921', 'line4713', 'line5801', 'line5411', 'line4541']

    custom_lines = []
    custom_files = []

    figsavefolder = r'D:\Peter\School\Master Thesis\figures\EWs\verslag\\'
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

    if obs is not 'BOTH':
        for chunk in chunks(lines,3):
            fig, axs = plt.subplots(nrows=len(chunk), ncols=1,sharex=True,figsize=(5, 8))
            fig.subplots_adjust(top = 0.83)
            size = fig.get_size_inches()
            savename = ''
            for i,ax in enumerate(axs):
                line = chunk[i]
                print line
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
                print lineinfo
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
            print savename
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
                print line
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
            print savename
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
        print lineinfo
        sgn = 91  # window size for SavitzkyGolay (must be odd integer)
        TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
        # print v[sgn]-v[0]
        p = chi2.ppf(0.99, n - 1) / (n - 1)
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
        print maxi, np.amax(spec3)
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
        print line[4:]
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
        print dt
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
                print mini, maxi
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

def plot_TVS_together(save='off',show='on',plot_save_folder='',oneline='off',sg='off', apowantedmarks = None,customlines = None,customfiles = None):
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
            elif i == 1:
                master_files = merc_master_files
                obs = 'MERCATOR'
            wl, TVS, v, n = airmass.TVS_masterfiles(master_files, line)
            lineinfo = getattr(master_files[0], line).lineinfo
            print lineinfo
            sgn = 91  # window size for SavitzkyGolay (must be odd integer)
            TVS_smoothed = SavitzkyGolay.savitzky_golay(TVS, sgn, 4)
            # print v[sgn]-v[0]
            p = chi2.ppf(0.99, n - 1) / (n - 1)
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
            print maxi, np.amax(spec3)
            axarr[0, i].set_ylim([mini, maxi])
            axarr[0, i].axvline(vsini, color='k', linestyle=':', linewidth=1)
            axarr[0, i].axvline(-vsini, color='k', linestyle=':', linewidth=1)
            # if line[2]==5875.621:
            #     TVS2 = np.array(TVS)*1.4
            axarr[1, i].plot(v, TVS, color='b')
            if sg == 'on':
                axarr[1, i].plot(v, TVS_smoothed, color='r', linestyle='dashed')
            if oneline == 'on':
                axarr[1, i].axhline(y=1, color='gray', linestyle='--')
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
            axarr[1, i].set_xlim([-600, 600])
            print line[4:]
            # plt.tight_layout()
            plt.subplots_adjust(left=0.08, bottom=None, right=0.98, top=None,
                            wspace=None, hspace=0.15)
        if save == 'on':
            plt.savefig(plot_save_folder + r'\\TVS_' +  lineinfo[0] + line[4:] + '_TVS.pdf', format='pdf',
                            dpi=1200)
        if show == 'on':
            plt.show()
        plt.close()

def plot_EW_together(save='off',show='on',plot_save_folder='',oneline='off',sg='off', apowantedmarks = None,customlines = None,customfiles = None, bad = 'off'):
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
            print line
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
                                                                                  np.array(ew_error))
            if bad =='on':
                chisq2, red_chisq2, AIC_12, AIC_22, AIC_32, probfactor2 = airmass.EW_stats2(np.array(ews), np.array(phases),
                                                                                      np.array(ew_error))
                add = '\n'+ str('LR: ' + '%.2E' % (1. / probfactor2))
            else:
                add = ''
            print lineinfo
            p1 = [0.005, 0.5, 1]
            fit1 = curve_fit(my_sin, phases, ews, p0=p1, sigma=ew_error, absolute_sigma=True)
            fit2 = curve_fit(airmass.flat_line, phases, ews, p0=[1], sigma=ew_error, absolute_sigma=True)
            t = np.linspace(0, 1, 50)
            data_fit = my_sin(t, *fit1[0])
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
        print savename
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


def plot_EW_together2(save='off',show='on',plot_save_folder='',oneline='off',sg='off', apowantedmarks = None,customlines = None,customfiles = None, bad = 'off'):
    # fig, axs = plt.subplots(nrows=len(chunk), ncols=1, sharex=True, figsize=(5, 8))
    # fig.subplots_adjust(top=0.83)
    # size = fig.get_size_inches()
    apo_lines = ['line6562', 'line4713', 'line5801', 'line4685', 'line5875', 'line4541', 'line5411', 'line5592',
                 'line4861', 'line4921', 'line6678', 'line4471']
    apo_master_files = open_masterfiles.apo()
    apo_testlist = glob.glob(r'D:\Peter\School\Master Thesis\Data\masterfiles\apo\test\\' + r'*.txt')
    apo_test_files = open_masterfiles.apo(manual_filelist=apo_testlist)
    merc_master_files = open_masterfiles.mercator()
    merc_testlist = glob.glob(r'D:\Peter\School\Master Thesis\Data\masterfiles\mercator\test\\' + r'*.txt')
    merc_test_files = open_masterfiles.mercator(manual_filelist=merc_testlist)
    vsini = 127

    for line in apo_lines:
    # for line in [apo_lines[1]]:
        f, axarr = plt.subplots(2, sharex=True, figsize=(9, 6))
        savename = ''
        for i in [0,1]:
            if i == 0:
                master_files = apo_master_files
                master_files2 = apo_test_files
                obs = 'APO'
            elif i == 1:
                master_files = merc_master_files
                master_files2 = merc_test_files
                obs = 'MERCATOR'
            print line
            phases = []
            ews = []
            ew_error = []
            for file in master_files:
                linedata = getattr(file, line)
                lineinfo = linedata.lineinfo
                phases.append(file.phase)
                ews.append(linedata.ew)
                ew_error.append(linedata.ew_error)
            ew1_avg = np.average(ews)
            norm_ews = np.array(ews)/ew1_avg
            norm_ew_error = np.array(ew_error)/ew1_avg
            phases2 = []
            ews2 = []
            ew_error2 = []
            for file in master_files2:
                linedata2 = getattr(file, line)
                lineinfo2 = linedata.lineinfo
                phases2.append(file.phase)
                ews2.append(linedata2.ew)
                ew_error2.append(linedata2.ew_error)
            ew2_avg = np.average(ews2)
            norm_ews2 = np.array(ews2)/ew2_avg
            norm_ew_error2 = np.array(ew_error2)/ew2_avg
            chisq, red_chisq, AIC_1, AIC_2, AIC_3, probfactor = airmass.EW_stats3(np.array(norm_ews), np.array(phases),
                                                                                  np.array(norm_ew_error))
            chisq2, red_chisq2, AIC_12, AIC_22, AIC_32, probfactor2 = airmass.EW_stats3(np.array(norm_ews2), np.array(phases2),
                                                                                  np.array(norm_ew_error2))
            if bad =='on':
                chisq2, red_chisq2, AIC_12, AIC_22, AIC_32, probfactor2 = airmass.EW_stats2(np.array(norm_ews), np.array(phases),
                                                                                      np.array(norm_ew_error))
                add = '\n'+ str('LR: ' + '%.2E' % (1. / probfactor2))
            else:
                add = ''
            print lineinfo
            p1 = [0.005, 0.5, 1]
            fit1 = curve_fit(my_sin, phases, norm_ews, p0=p1, sigma=norm_ew_error, absolute_sigma=True)
            fit2 = curve_fit(airmass.flat_line, phases, norm_ews, p0=[1], sigma=norm_ew_error, absolute_sigma=True)
            t = np.linspace(0, 1, 50)
            data_fit = my_sin(t, *fit1[0])
            axarr[i].set_title(obs, size=14, weight='light')
            axarr[i].errorbar(phases, norm_ews, yerr=norm_ew_error, fmt='o', c='r', label='Data')
            axarr[i].plot(t, data_fit, c='r', label='Sinusodial fit')
            axarr[i].axhline(fit2[0][0], linestyle='--', c='gray', label='Weighted Average')

            fit12 = curve_fit(my_sin, phases2, norm_ews2, p0=p1, sigma=norm_ew_error2, absolute_sigma=True)
            fit22 = curve_fit(airmass.flat_line, phases2, norm_ews2, p0=[1], sigma=norm_ew_error2, absolute_sigma=True)
            data_fit2 = my_sin(t, *fit12[0])
            # axarr[i].set_title(obs, size=14, weight='light')
            axarr[i].errorbar(phases2, norm_ews2, yerr=norm_ew_error2, fmt='o', c='b', label='Data2')
            axarr[i].plot(t, data_fit2, c='b', label='Sinusodial fit 2')
            axarr[i].axhline(fit22[0][0], linestyle='--', c='gray', label='Weighted Average 2')

            aldat = np.append(norm_ews, data_fit)
            # ax.locator_params(axis='y', nbins=7)
            lowlim, uplim = np.round(np.min(aldat), decimals=2) - 0.01, np.round(np.max(aldat), decimals=2) + 0.01
            # l4 = axarr[i].text(0.78, (0.01 * (uplim - lowlim)) + lowlim, str('LR: ' + '%.2E' % (1. / probfactor))+add, label='likelihood ratio',size=14)
            l4 = axarr[i].text(0.78, (0.01 * (uplim - lowlim)) + lowlim, str('LR: ' + '%.2E' % (1. / probfactor)) + '\n'+ str('LR: ' + '%.2E' % (1. / probfactor2)),
                               label='likelihood ratio', size=14)
            # l4 = axarr[i].text( str('%.2E' % (1. / probfactor)),
            #                    label='likelihood ratio')
            # axarr[i].set_ylim( lowlim,uplim)
        savename += lineinfo[0] + line[-4:] + '_'
        savename += '_EW.pdf'
        print savename
        # plt.figlegend((l1, l2, l3), ('Data', 'Sinusodial fit', 'Weighted Average'), bbox_to_anchor=[0.87, 0.85],
        #               prop={'size': 11}, fancybox=True, framealpha=0.8)
        plt.suptitle('Normalized Equivalent Width \n' + lineinfo[-1], size=18)
        plt.subplots_adjust(top=0.85)
        plt.xlabel('Phase (P = 6.83d)',fontsize = 13)
        # plt.tight_layout()
        if save =='on':
            plt.savefig(plot_save_folder + savename)
        if show == 'on':
            plt.show()
        plt.close()



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
                print lineinfo[-1], r'  & \ref{fig:EW_' + lineinfo[-1][-4:] +r'} &',
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

            print pf, '   &', det,'   &',
        print r'\\'
# apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592', 'line4861', 'line4921', 'line6678', 'line4471']
# mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471', 'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']
# sort_linelist(apo_lines, obs='APO',apowantedmarks = None)


# plot_quotient(obs='APO',overplot='on', oneline='on', show='off', save='on', plot_save_folder=r'D:\Peter\School\Master Thesis\figures\Quotients\apo\overplot')
# plot_quotient(obs='MERCATOR',overplot='on', oneline='on', show='off', save='on', plot_save_folder=r'D:\Peter\School\Master Thesis\figures\Quotients\mercator\overplot')
# plot_TVS(obs = 'APO',oneline='on',plot_save_folder=r'D:\Peter\School\Master Thesis\figures\TVS\eShel\from_masterfiles',save='on',show='off')
# plot_TVS(obs = 'MERCATOR',oneline='on',plot_save_folder=r'D:\Peter\School\Master Thesis\figures\TVS\LaPalma\from_masterfiles',save='on',show='off')
# plot_EW(obs='MERCATOR',apowantedmarks = None)
# plot_EW(obs='APO',apowantedmarks = None)
# plot_TVS_together(plot_save_folder= r'D:\Peter\School\Master Thesis\figures\TVS\verslag', save = 'on', show='off', oneline='on')
# plot_EW_together(plot_save_folder= r'D:\Peter\School\Master Thesis\figures\EWs\verslag\EW_', save = 'off', show='on', oneline='on', bad='off')
plot_EW_together2(plot_save_folder= r'D:\Peter\School\Master Thesis\figures\EWs\verslag\EW_', save = 'off', show='on', oneline='on', bad='off')
# EW_fit_table()