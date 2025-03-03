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

apo_eshel_files= open_masterfiles.apo_demetra_orders(r"D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\all_darks\rebin01\combined\high_snr\\")
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
mercfile = mercator_files[1]
m_wl = mercfile.wl_rebin
m_flux = mercfile.flux_rebin
m_wl_rebin,m_flux_rebin = airmass.rebin2(m_wl,m_flux)
m_wl_deg,m_flux_deg = airmass.degrade_spectrum(m_wl,m_flux,pre_rebin=0.1,desired_snr=120)

m_spec1d_deg = specutils.Spectrum1D(spectral_axis=m_wl_deg * AA, flux=m_flux_deg * ap_unit.Jy)
print('specutils snr',analysis.snr_derived(m_spec1d_deg,specutils.SpectralRegion(bd_ha[0]*AA, bd_ha[1]*AA)))

m_snr_straight = airmass.snr_2(m_wl_deg,m_flux_deg,boundaries=bd,rebin=False,rebin_size=0.1,separate=False)
m_snr_ha = airmass.snr_2(m_wl_deg,m_flux_deg,boundaries=bd_ha,rebin=False,rebin_size=0.1,separate=False)

apo_test_file = apo_eshel_files[0]
apo_order_straight_line = airmass.find_order(bd,apo_test_file)

#write code to measure snr in apo spectra
apo_sl_wl = apo_order_straight_line.wl_rebin[2:]
apo_sl_flux = apo_order_straight_line.flux_rebin[2:]/np.average(apo_order_straight_line.flux_rebin[2:])
m_slice_wl, m_slice_flux = slice_and_norm(m_wl_deg,m_flux_deg,apo_sl_wl[0],apo_sl_wl[-1])
print(m_snr_straight,m_snr_ha)
plt.plot(apo_sl_wl,apo_sl_flux)
plt.plot(m_slice_wl,m_slice_flux)
plt.ylim(0.8,1.1)
plt.axvline(bd[0])
plt.axvline(bd[1])
plt.show()
plt.close()












# show_specutils_fit(apo_eshel_files[0])
# plot_normalization_test(apo_eshel_files,mercator_files)