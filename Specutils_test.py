# __author__ = 'PeterBacily'
from __future__ import division
import warnings
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
import open_masterfiles
import scipy.stats as ss
from scipy.optimize import *
from PyAstronomy import pyasl
import Path_check
import os
import specutils
from specutils.fitting import fit_generic_continuum
from specutils import analysis
import astropy.units as ap_unit
import pickle

AA = ap_unit.AA
linelist = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
             'line4861', 'line4921', 'line6678', 'line4471']
apo_eshel_files_combined = open_masterfiles.apo_demetra_orders(r"D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin01\combined\\")
apo_eshel_files= open_masterfiles.apo_demetra_orders(r"D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin01\combined\high_snr\\")
apo_eshel_files_single_obs= open_masterfiles.apo_demetra_orders(r"D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin01\single_obs\\")
mercator_files = open_masterfiles.mercator(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\mercator\ll_apo_vcor_2\\')



def slice_and_norm(wl,flux,start,end,rebin=None):
    slice_flux = flux[(wl > start) & (wl < end)]
    slice_wl = wl[(wl > start) & (wl < end)]
    slice_flux_norm = slice_flux / np.average(slice_flux)
    if rebin == None:
        return slice_wl,slice_flux_norm
    else:
        slice_wl_rebinned,slice_flux_rebinned = airmass.rebin2(slice_wl,slice_flux_norm,step=rebin)
        return slice_wl_rebinned,slice_flux_rebinned


def show_specutils_fit(file):
    line = linelist[0]
    list_of_orders = file.orders
    ha_order = list(filter(lambda x: x.order_number_demetra == '34', list_of_orders))[0]
    wl = ha_order.wl_rebin
    flux = ha_order.flux_rebin
    spec1d = specutils.Spectrum1D(spectral_axis=wl * AA, flux=flux * ap_unit.Jy)
    linekey = line + '_order_rebin'
    line_instance_order = getattr(file, linekey)
    lineinfo = line_instance_order.lineinfo
    norm_boundaries = [6538,6546,6575,6589]
    region = [(norm_boundaries[0] * AA, norm_boundaries[1] * AA), (norm_boundaries[2] * AA, norm_boundaries[3] * AA)]
    g1_fit = specutils.fitting.fit_generic_continuum(spec1d)
    continuum_fitted = g1_fit(wl * ap_unit.AA)

    plt.plot(wl,flux, label = 'Data')
    plt.plot(wl,continuum_fitted, label = 'Fitted continuum')
    plt.xlabel('Wavelength (AA)')
    plt.ylabel('Relative Flux')
    plt.legend(loc='lower right')
    plt.ylim(0.8, 1.15)
    plt.xlim(6530,6590)
    # for bound in norm_boundaries:
    #     plt.axvline(bound, color='black',linestyle='--',linewidth=1)
    plt.title('Specutils Continuum fit')
    plt.show()
    plt.close()

def plot_normalization_test(apo_eshel_files,mercator_files):
    f, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, sharex=True)
    line = linelist[0]
    for file in apo_eshel_files:
        linekey = line+'_order_rebin'
        line_instance_order = getattr(file, linekey)
        lineinfo = line_instance_order.lineinfo
        norm_boundaries = lineinfo[2:6]
        ax1.axvspan(norm_boundaries[0] , norm_boundaries[1], facecolor='0.95', edgecolor='0', linestyle='--', alpha=1)
        ax1.axvspan(norm_boundaries[2] , norm_boundaries[3], facecolor='0.95', edgecolor='0', linestyle='--', alpha=1)
        wl = line_instance_order.wl
        flux = line_instance_order.flux
        ax1.plot(wl,flux)
        start = wl[0]
        stop = wl[-1]


    for file in apo_eshel_files:
        list_of_orders = file.orders
        ha_order = list(filter(lambda x: x.order_number_demetra == '34', list_of_orders))[0]
        wl = ha_order.wl_rebin
        flux=ha_order.flux_rebin
        spec1d = specutils.Spectrum1D(spectral_axis=wl* ap_unit.AA, flux=flux*ap_unit.Jy)
        ax2.plot(spec1d.spectral_axis, spec1d.flux)
        norm_boundaries = [6538,6546,6575,6589]
        region = [(norm_boundaries[0] * AA, norm_boundaries[1] * AA), (norm_boundaries[2] * AA, norm_boundaries[3] * AA)]
        g1_fit = specutils.fitting.fit_continuum(spec1d)
        continuum_fitted = g1_fit(wl * AA)
        g1_fit = fit_generic_continuum(spec1d)
        y_continuum_fitted = g1_fit(wl * AA)
        spec_normalized = spec1d/ y_continuum_fitted
        ax3.plot(spec_normalized.spectral_axis, spec_normalized.flux)
        ax3.axvspan(norm_boundaries[0] , norm_boundaries[1], facecolor='0.95', edgecolor='0', linestyle='--', alpha=1)
        ax3.axvspan(norm_boundaries[2] , norm_boundaries[3], facecolor='0.95', edgecolor='0', linestyle='--', alpha=1)

    for file in mercator_files:
        linekey = line+'_rebin'
        line_instance_order = getattr(file, linekey)
        lineinfo = line_instance_order.lineinfo
        norm_boundaries = lineinfo[2:6]
        ax4.axvspan(norm_boundaries[0] , norm_boundaries[1], facecolor='0.95', edgecolor='0', linestyle='--', alpha=1)
        ax4.axvspan(norm_boundaries[2] , norm_boundaries[3], facecolor='0.95', edgecolor='0', linestyle='--', alpha=1)
        wl = line_instance_order.wl
        flux = line_instance_order.flux
        ax4.plot(wl,flux)


    ax1.set_title('APO, Normalized with Linear Fit')
    ax2.set_title('APO, Not Normalized')
    ax3.set_title('APO, Normalized with Specutils')
    ax4.set_title('Mercator, normalized with Linear Fit')
    f.supxlabel('Wavelength (AA)')
    f.supylabel('Relative Flux')
    f.suptitle('Effect of changes in normalization method')
    ax1.set_xlim(start,stop)
    ax2.set_xlim(start,stop)
    ax3.set_xlim(start,stop)
    ax4.set_xlim(start,stop)
    ax1.set_ylim(0.8,1.15)
    ax2.set_ylim(0.8,1.15)
    ax3.set_ylim(0.8,1.15)
    ax4.set_ylim(0.8,1.15)
    plt.show()
    plt.close()


bd = [5170, 5190]
bd_ha = [6614.0, 6625.0]


def SNR_3(wl,flux,boundaries='Halpha',rebin=False,separate=False):
    if boundaries == 'Halpha':
        bd = [6614.0, 6625.0]
    elif boundaries == 'flat_continuum':
        bd = [5170, 5190]
    elif (np.array(boundaries).shape ==(4,) or np.array(boundaries).shape ==(2,)):
        bd = boundaries
    else:
        raise TypeError('Boundaries needs to be \'Halpha\', \'flat_continuum\', or a list of 2 or 4 boundaries')
    snr = airmass.snr_2(wl, flux, boundaries=bd, rebin=rebin, rebin_size=0.1, separate=separate)
    return snr

def SNR_merc(masterfile):
    m_wl = masterfile.wl_rebin2
    m_flux = masterfile.flux_rebin2
    snr_ha = SNR_3(m_wl,m_flux,boundaries='Halpha',rebin=False,separate=False)
    snr_straight = SNR_3(m_wl,m_flux,boundaries='flat_continuum',rebin=False,separate=False)
    return snr_ha,snr_straight

def SNR_merc_degen(wl,flux):
    m_wl = wl
    m_flux = flux
    snr_ha = SNR_3(m_wl,m_flux,boundaries='Halpha',rebin=False,separate=False)
    snr_straight = SNR_3(m_wl,m_flux,boundaries='flat_continuum',rebin=False,separate=False)
    return snr_ha,snr_straight


def SNR_apo_orders(file):
    bd_ha = [6614, 6625]
    bd_straight = [5170, 5190]
    apo_order_ha = airmass.find_order(bd_ha, file)
    apo_order_straight_line = airmass.find_order(bd_straight, file)
    apo_sl_wl = apo_order_straight_line.wl_rebin[2:-3]
    apo_sl_flux = apo_order_straight_line.flux_rebin[2:-3]
    apo_ha_wl = apo_order_ha.wl_rebin[2:-3]
    apo_ha_flux = apo_order_ha.flux_rebin[2:-3]
    snr_ha = SNR_3(apo_ha_wl,apo_ha_flux,boundaries='Halpha',rebin=False,separate=False)
    snr_straight = SNR_3(apo_sl_wl,apo_sl_flux,boundaries='flat_continuum',rebin=False,separate=False)
    return snr_ha,snr_straight

mercfile = mercator_files[1]
m_wl = mercfile.wl_rebin2
m_flux = mercfile.flux_rebin2
print(m_wl[1001]-m_wl[1000])
m_wl_rebin,m_flux_rebin = airmass.rebin2(m_wl,m_flux)
m_wl_deg,m_flux_deg = airmass.degrade_spectrum(m_wl,m_flux,pre_rebin=0.1,spectral_resolution=10000,desired_snr=120)
m_spec1d_deg = specutils.Spectrum1D(spectral_axis=m_wl_deg * AA, flux=m_flux_deg * ap_unit.Jy)
m_spec1d_normal = specutils.Spectrum1D(spectral_axis=m_wl_rebin * AA, flux=m_flux_rebin * ap_unit.Jy)





m_snr_straight = airmass.snr_2(m_wl_deg,m_flux_deg,boundaries=bd,rebin=False,rebin_size=0.1,separate=False)
m_snr_ha = airmass.snr_2(m_wl_deg,m_flux_deg,boundaries=bd_ha,rebin=False,rebin_size=0.1,separate=False)
m_snr_specutils_ha = analysis.snr_derived(m_spec1d_deg,specutils.SpectralRegion(bd_ha[0]*AA, bd_ha[1]*AA))
m_snr_specutils_straight = analysis.snr_derived(m_spec1d_deg,specutils.SpectralRegion(bd[0]*AA, bd[1]*AA))

# single_obs_inspect_list = [apo_eshel_files_single_obs[2],apo_eshel_files_single_obs[16],apo_eshel_files_single_obs[25],apo_eshel_files_single_obs[33]]
# single_obs_inspect_list = [apo_eshel_files_single_obs[2],apo_eshel_files_single_obs[33]]
i=0
snrs_ha=[]
snrs_str=[]

# for file in mercator_files:
#     # snr_ha,snr_str= SNR_merc(file)
#     m_wl = file.wl_rebin2
#     m_flux = file.flux_rebin2
#     m_spec1d = specutils.Spectrum1D(spectral_axis=m_wl * AA, flux=m_flux * ap_unit.Jy)
#     snr_ha=m_snr_specutils_ha = analysis.snr_derived(m_spec1d,specutils.SpectralRegion(bd_ha[0]*AA, bd_ha[1]*AA))
#     snr_str=analysis.snr_derived(m_spec1d,specutils.SpectralRegion(bd[0]*AA, bd[1]*AA))
#     snrs_ha.append(snr_ha)
#     snrs_str.append(snr_str)

deg_list=  np.arange(50, 200, 10)
m_wl = mercfile.wl_rebin2
m_flux = mercfile.flux_rebin2
for deg in deg_list:
    # snr_ha,snr_str= SNR_merc(file)
    deg_wl,deg_flux= airmass.degrade_spectrum(m_wl,m_flux,desired_snr=deg)
    # m_spec1d = specutils.Spectrum1D(spectral_axis=deg_wl * AA, flux=deg_flux * ap_unit.Jy)
    # snr_ha=m_snr_specutils_ha = analysis.snr_derived(m_spec1d,specutils.SpectralRegion(bd_ha[0]*AA, bd_ha[1]*AA))
    # snr_str=analysis.snr_derived(m_spec1d,specutils.SpectralRegion(bd[0]*AA, bd[1]*AA))
    snr_ha,snr_str = SNR_merc_degen(deg_wl,deg_flux)
    snrs_ha.append(snr_ha)
    snrs_str.append(snr_str)
i=0
apo_eshel_files
for apo_test_file in apo_eshel_files_combined:
# for apo_test_file in [apo_eshel_files_single_obs[19]]:
# for apo_test_file in [apo_eshel_files[3]]:
    apo_order_straight_line = airmass.find_order(bd,apo_test_file)
    apo_order_ha = airmass.find_order(bd_ha,apo_test_file)
    apo_sl_wl = apo_order_straight_line.wl_rebin[2:]
    apo_sl_flux = apo_order_straight_line.flux_rebin[2:]/np.average(apo_order_straight_line.flux_rebin[2:])
    apo_spec1d_sl = specutils.Spectrum1D(spectral_axis=apo_sl_wl * AA, flux=apo_sl_flux * ap_unit.Jy)
    apo_ha_wl = apo_order_ha.wl_rebin[2:]
    apo_ha_flux = apo_order_ha.flux_rebin[2:]/np.average(apo_order_ha.flux_rebin[2:])
    apo_spec1d_ha = specutils.Spectrum1D(spectral_axis=apo_ha_wl * AA, flux=apo_ha_flux * ap_unit.Jy)
    snr_ha,snr_straight = SNR_apo_orders(apo_test_file)
    print('No.',i, 'SNR:',str(np.round(snr_ha,1)))
    i+=1
exit()
deg_merc_wl,deg_merc_flux = slice_and_norm(m_wl_deg,m_flux_deg,start=apo_ha_wl[2],end=apo_ha_wl[-3])

plt.plot(apo_ha_wl,apo_ha_flux,label = 'APO, SNR = '+str(np.round(snr_ha,0)))
plt.plot(deg_merc_wl,deg_merc_flux,label='Mercator degraded, SNR = 120')
plt.xlim(6614,6625)
plt.ylim(0.9,1.1)
plt.xlabel('Wavelength (Å)')
plt.ylabel('Relative Flux')
plt.legend()
plt.show()
plt.close()
#
#
#     a_snr_straight = airmass.snr_2(apo_sl_wl,apo_sl_flux,boundaries=bd,rebin=False,rebin_size=0.1,separate=False)
#     a_snr_ha = airmass.snr_2(apo_ha_wl,apo_ha_flux,boundaries=bd_ha,rebin=False,rebin_size=0.1,separate=False)
#     a_snr_specutils_ha = analysis.snr_derived(apo_spec1d_ha,specutils.SpectralRegion(bd_ha[0]*AA, bd_ha[1]*AA))
#     a_snr_specutils_straight = analysis.snr_derived(apo_spec1d_sl,specutils.SpectralRegion(bd[0]*AA, bd[1]*AA))
#     # print('Mercator degraded:')
#     # print('specutils snr straight: ',m_snr_specutils_straight,' ha: ',m_snr_specutils_ha)
#     # print('my snr straight: ',m_snr_straight,' ha: ',m_snr_ha)
#     snr_ha_2,snr_str_2 = SNR_apo_orders(apo_test_file)
#     # snrs_ha.append((snr_ha_2))
#     # snrs_str.append((snr_str_2))
#     snrs_ha.append((a_snr_specutils_ha))
#     snrs_str.append((a_snr_specutils_straight))
#     print(i,' APO:')
#     print('specutils snr:    straight:',a_snr_specutils_straight,' ha:',a_snr_specutils_ha)
#     print('fu snr:    straight: ', snr_str_2, ' ha: ', snr_ha_2, )
#     print('my snr:    straight: ',a_snr_straight,' ha: ',a_snr_ha,'\n')
#
#     m_slice_wl_sl, m_slice_flux_sl = slice_and_norm(m_wl_deg,m_flux_deg,apo_sl_wl[0],apo_sl_wl[-1])
#     m_slice_wl_ha, m_slice_flux_ha = slice_and_norm(m_wl_deg, m_flux_deg, apo_ha_wl[0], apo_ha_wl[-1])
#     i+=1

# plt.scatter(snrs_ha,snrs_str)
# plt.title('SNR from a degenerated Mercator spectra derived with a straight line fit')
# plt.xlabel('SNR Hα 6614Å')
# plt.ylabel('SNR flat continuum 5180Å')
# # plt.title('SNR from Mercator spectra in flat continuum 5180Å')
# # plt.xlabel('SNR from straight line fit')
# # plt.ylabel('SNR from Specutils')
# plt.show()
# plt.close()
    # plt.plot(m_slice_wl_ha,m_slice_flux_ha, label='mercator degraded to 120SNR')
#     plt.plot(apo_ha_wl,apo_ha_flux,label = str(i))
#     plt.ylim(0.8,1.1)
# plt.axvline(bd_ha[0])
# plt.axvline(bd_ha[1])
# plt.legend()
# plt.show()
# plt.close()

# print('Mercator degraded:')
# print('specutils snr: straight: ',m_snr_specutils_straight,' ha: ',m_snr_specutils_ha)
# print('my snr: straight: ',m_snr_straight,' ha: ',m_snr_ha)











# show_specutils_fit(apo_eshel_files[0])
# plot_normalization_test(apo_eshel_files,mercator_files)