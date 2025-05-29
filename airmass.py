from __future__ import division

import time
import matplotlib.pyplot as plt
import glob
import pickle
import warnings
import astropy.io.fits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
from scipy import interpolate
from scipy.optimize import *

from PyAstronomy import pyasl
import scipy.stats as ss
import scipy.signal as scisig
from collections import defaultdict

warnings.simplefilter('ignore')
from SavitzkyGolay import savitzky_golay

from pysynphot import observation
from pysynphot import spectrum
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
import astropy.units as u
from astropy.timeseries import LombScargle
warnings.resetwarnings()
from linedict import linedict
import Datafile_class
c_light = 299792.458
# filelist = glob.glob('C:/peter/School/Master Scriptie/Data/data/*.fit')
# print filelist
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=PendingDeprecationWarning)
import Omar_functies_to_use
import bjd_from_jd

def my_sin(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*x  + phase)  + offset
def my_sin2(x,  amplitude, phase, offset):
    return amplitude*np.sin(4*np.pi*x  + phase)  + offset

def my_line(x,a,b):
    return a*x+b

def flat_line(x,a):
    return a

# this comment is to test if pushing work

def redchisqg(ydata,ymod,deg=2,sd=None):
      """
 Returns the reduced chi-square error statistic for an arbitrary model,
 chisq/nu, where nu is the number of degrees of freedom. If individual
 standard deviations (array sd) are supplied, then the chi-square error
 statistic is computed as the sum of squared errors divided by the standard
 deviations. See http://en.wikipedia.org/wiki/Goodness_of_fit for reference.

 ydata,ymod,sd assumed to be Numpy arrays. deg integer.

 Usage:
 # >>> chisq=redchisqg(ydata,ymod,n,sd)
 where
  ydata : data
  ymod : model evaluated at the same x points as ydata
  n : number of free parameters in the model
  sd : uncertainties in ydata

 Rodrigo Nemmen
 http://goo.gl/8S1Oo
       """
      # Chi-square statistic
      if sd==None:
           chisq=np.sum((ydata-ymod)**2)
      else:
           chisq=np.sum( ((ydata-ymod)/sd)**2 )

      # Number of degrees of freedom assuming 2 free parameters
      nu=ydata.size-1-deg

      return chisq/nu

def aphase(HJD):
    # datafile = pf.open(file)
    period = 6.829621
    jdstart = 2454391.
    # jdstart = 2452734.2
    # header =  datafile[0].header
    # filetime = Time(header['MJD-OBS'], format='jd').jd
    phase = ((HJD- jdstart)% period )/ period
    return phase

def barcor(file, JDOBS=None):
    datafile = pf.open(file)
    if JDOBS == None:
        JD = 2400000.5 + datafile[0].header['MJD-OBS']
    else:
        JD = JDOBS
    DEC = (-1.9425736)
    RA = 85.18
    LAT = 52.354094
    LON = 4.955252
    ALT = 20
    bcor = pyasl.helcorr(LON, LAT, ALT, RA, DEC, JD, debug=False)
    datafile.close()
    return bcor[0], bcor[1]

def airmass(file, JDOBS=None):
    datafile = pf.open(file)
    if JDOBS == None:
        JD =  2400000.5+datafile[0].header['MJD-OBS']
    else:
        JD = JDOBS
    DEC = (-1.94)* (2*math.pi/360)
    RA = 85.18* (2*math.pi/360)
    LAT = 52.354094* (2*math.pi/360)
    LON = 4.955252* (2*math.pi/360)
    D = JD - 2451545.0
    GMST = ((18.697374558 + 24.06570982441908*D)%24) * (360/24)* (2*math.pi/360)
    sinalt = math.sin(DEC)*math.sin(LAT) + math.cos(DEC)*math.cos(LAT)*math.cos((GMST + LON - RA))
    alt = np.arcsin(sinalt)*360/(2*math.pi)
    airmass = 1/sinalt
    datafile.close()
    return airmass,alt,JD

def bjd_lapalma_from_date_zet_ori(header_time_and_date):
    coordinates = SkyCoord("05h40m45.50s -01d56m34.30s", frame=Galactic)
    timeobject = Time(header_time_and_date, format='isot')
    bjd= bjd_from_jd.bjd_tdb(jd_utc=timeobject.jd,sky_position=coordinates)
    bjd_float=float(bjd.to_value('jd', subfmt='str'))
    return bjd_float
def bjd_lapalma_from_date_zet_ori_omar(header_time_and_date):
    coord_zet_ori = "05h40m45.50s -01d56m34.30s"
    (bjd,bccor) = Omar_functies_to_use.BJD_calculator(header_time_and_date,coord_zet_ori)
    return bjd,bccor
def timeanddate(file):
    datafile = pf.open(file)
    header =  datafile[0].header
    dt = header['DATE-OBS']
    MJD = header['MJD-OBS']
    yr = dt[0:4]
    m = int(dt[5:7])
    d = dt[8:10]
    t = dt[11:16]
    m_n = calendar.month_abbr[m]
    datafile.close()
    return str(d) + ' ' + str( m_n) + ' ' + str(t),MJD

def timeanddate2(DATE_OBS):
    dt = DATE_OBS
    yr = dt[0:4]
    m = int(dt[5:7])
    d = dt[8:10]
    t = dt[11:16]
    m_n = calendar.month_abbr[m]
    return str(d) + ' ' + str( m_n) + ' ' + str(t)

def date(file):
    datafile = pf.open(file)
    header =  datafile[0].header
    dt = header['DATE-OBS']
    MJD = header['MJD-OBS']
    yr = dt[0:4]
    m = dt[5:7]
    d = dt[8:10]
    t = dt[11:16]
    # m_n = calendar.month_abbr[m]
    datafile.close()
    return yr+m+d

def timeanddatelp(file):
    datafile = pf.open(file)
    header =  datafile[0].header
    dt = header['DATE-OBS']
    # MJD = header['MJD-OBS']
    yr = dt[0:4]
    m = int(dt[5:7])
    d = dt[8:10]
    t = dt[11:16]
    m_n = calendar.month_abbr[m]
    datafile.close()
    return str(d) + ' ' + str( m_n) + ' ' + str(t)


def split_date(DATE_OBS):
    dt = DATE_OBS
    # MJD = header['MJD-OBS']
    yr = int(dt[0:4])
    yr_str = str(yr)
    m = int(dt[5:7])
    m_str = '{num:02d}'.format(num=m)
    d = int(dt[8:10])
    d_str = '{num:02d}'.format(num=d)
    # print dt
    # print dt[11:13],dt[14:16]
    t = int(dt[11:13]+dt[14:16])
    t_str = '{num:04d}'.format(num=t)
    # m_n = calendar.month_abbr[m]
    return [yr_str,m_str,d_str,t_str],[yr,m,d,t]


def phase(file):
    datafile = pf.open(file)
    period = 6.829
    jdstart = Time(2454391., format='jd').jd
    header =  datafile[0].header
    filetime = Time(header['MJD-OBS'], format='jd').jd
    phase = ((filetime- jdstart)% period )/ period
    datafile.close()
    return phase


def extractdata(j,file,return_header = 'off'):
    datafile = pf.open(file)
    data = datafile[j].data
    header = datafile[j].header
    currentwl = header['CRVAL1']
    step =header['CDELT1']
    wl=[]
    for i in range(len(data)):
        wl.append(currentwl)
        currentwl+=step
    wl2 = np.array(wl)
    datafile.close()
    if return_header == 'on':
        return wl2, data, header
    else:
        return wl2, data
    pf.close(file)

def grab_norm_area(wave,flux,a,b,c,d ):
    normwave = np.hstack((wave[(wave > a) & (wave < b)], wave[(wave > c) & (wave < d)]))
    normflux = np.hstack((flux[(wave > a) & (wave < b)], flux[(wave > c) & (wave < d)]))
    return normwave,normflux

def group_by_day(list_of_masterfile_objects):
    groups = defaultdict(list)
    for obj in list_of_masterfile_objects:
        groups[obj.time_and_date[0:5]].append(obj)
    new_list = groups.values()
    return new_list

def AIC_rel_likelihood(AIC_1, AIC_2):
    if AIC_1 < AIC_2:
        AIC_min = AIC_1
        AIC_max =AIC_2
        return np.exp((AIC_min-AIC_max)/2)
    else:
        AIC_min = AIC_2
        AIC_max =AIC_1
        print('second AIC is lower, so they are reversed')
        return 0


def exposuretime(file):
    datafile = pf.open(file)
    header = datafile[0].header
    datafile.close()
    return header['EXPOSURE']


def HJD_rounded(file):
    BCCor, HJD = barcor(file)
    return round(HJD - 2457000, 3)

def normalize(wls,flux,a,b,c,d,startwl,endwl,xtype = 'wave',linecenter = None):
    if xtype == 'wave':
        wave = wls
    if xtype == 'velo':
        wave,vsini = wl_to_velocity(wls,linecenter)
    normwave = np.hstack((wave[(wave>a)&(wave<b)],wave[(wave>c)&(wave<d)]))
    normflux = np.hstack((flux[(wave>a)&(wave<b)],flux[(wave>c)&(wave<d)]))
    #fit line trough slice
    # print('ss', startwl, endwl)
    slope,height = np.polyfit(normwave,normflux,1)
    # print 'slope and height are', slope, height
    fit = np.poly1d([slope,height])
    linewave = wave[(wave>startwl)&(wave<endwl)]


    # print('linewave=', linewave)
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

    return linewave, fluxarray,nnf,lineflux, fit


def normalize_single_part(x,y,start=None,stop=None):
    slope, height = np.polyfit(x, y, 1)
    # print 'slope and height are', slope, height
    fit = np.poly1d([slope, height])
    if start == None and stop ==None:
        xslice=x
        yslice=y
    else:
        xslice = x[(x > start) & (x < stop)]
        yslice= y[(x > start) & (x < stop)]
    norm_y = []
    # print linewave
    for i, j in enumerate(xslice):
        norm_y.append(yslice[i] / fit(j))
    fluxarray = np.array(norm_y)
    return xslice, fluxarray, fit


def snr(wl,flux):

    # wl,data = extractdata(35, file)
    # print wl
    start =5375
    stop = 5385
    # print wl[(wl>start-5) & (wl<stop+5)]
    slice = flux[(wl>start) & (wl<stop)]
    wlslice = wl[(wl>start) & (wl<stop)]

    l = int((len(wlslice)-1)/2)
    lw,lf, _,_,_ = normalize(wlslice,slice,wlslice[0],wlslice[l],wlslice[l+1],wlslice[-1],wlslice[0],wlslice[-1])
    # plt.plot(lw,lf)
    # plt.pause(0.5)
    avgcounts = np.average(lf)
    stand_dev = np.std(lf)
    stnr = avgcounts/stand_dev
    return stnr




def snr_2(wl,flux,boundaries=[6549, 6550.7, 6576.0, 6578.0],rebin=True,rebin_size=0.1,separate=False):
    if len(boundaries)==2:
        [a,d]=boundaries
        separate = False
    elif len(boundaries)==4:
        [a,b,c,d]=boundaries
    else:
        raise Exception('boundaries  needs to be a list of either 2 or 4 values')
    start =a-2
    stop = d+2
    slice = flux[(wl>start) & (wl<stop)]
    wlslice = wl[(wl>start) & (wl<stop)]

    if len(boundaries)==2:
        l = int((len(wlslice) - 1) / 2)
        b=l
        c=l+1
    if rebin is True:
        wl_rebin = np.arange(wlslice[10], wlslice[-10], rebin_size)
        flux_rebin = rebin_spec(wlslice, slice, wl_rebin)
        lw,lf, _,_,_ = normalize(wl_rebin,flux_rebin,a,b,c,d,wl_rebin[0],wl_rebin[-1])
    else:
        lw,lf, _,_,_ = normalize(wlslice,slice,a,b,c,d,a,d)
    if separate is False:
        normwave = np.hstack((lw[(lw > a) & (lw < b)], lw[(lw > c) & (lw < d)]))
        normflux = np.hstack((lf[(lw > a) & (lw < b)], lf[(lw > c) & (lw < d)]))
        avgcounts = np.average(normflux)
        stand_dev = np.std(normflux)
        stnr = avgcounts/stand_dev
    if separate is True:
        left_part_wl = lw[(lw > a) & (lw < b)]
        left_part_flux = lf[(lw > a) & (lw < b)]
        right_part_wl = lw[(lw > c) & (lw < d)]
        right_part_flux =  lf[(lw > c) & (lw < d)]
        left_part_normwave,left_part_normflux,left_fit  = normalize_single_part(left_part_wl,left_part_flux)
        right_part_normwave, right_part_normflux, right_fit = normalize_single_part(right_part_wl, right_part_flux)
        stds = []
        weights = []
        for normflux in [left_part_normflux,right_part_normflux]:
            avgcounts = np.average(normflux)
            stand_dev = np.std(normflux)
            sr_part = avgcounts/stand_dev
            stds.append(sr_part)
            weights.append(len(normflux))
        stnr = np.average(stds,weights=weights)

    return stnr


def SNR_3(wl,flux,boundaries='Halpha',rebin=False,separate=False):
    if boundaries == 'Halpha':
        bd = [6614.0, 6625.0]
    elif boundaries == 'flat_continuum':
        bd = [5170, 5190]
    elif (np.array(boundaries).shape ==(4,) or np.array(boundaries).shape ==(2,)):
        bd = boundaries
    else:
        raise TypeError('Boundaries needs to be \'Halpha\', \'flat_continuum\', or a list of 2 or 4 boundaries')
    snr = snr_2(wl, flux, boundaries=bd, rebin=rebin, rebin_size=0.1, separate=separate)
    return snr

def SNR_merc(masterfile,binsize = '01'):
    m_wl = getattr(masterfile,'wl'+'_rebin'+binsize)
    m_flux = getattr(masterfile,'flux'+'_rebin'+binsize)
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
    apo_order_ha = find_order(bd_ha, file)
    apo_order_straight_line = find_order(bd_straight, file)
    apo_sl_wl = apo_order_straight_line.wl_rebin[2:-3]
    apo_sl_flux = apo_order_straight_line.flux_rebin[2:-3]
    apo_ha_wl = apo_order_ha.wl_rebin[2:-3]
    apo_ha_flux = apo_order_ha.flux_rebin[2:-3]
    snr_ha = SNR_3(apo_ha_wl,apo_ha_flux,boundaries='Halpha',rebin=False,separate=False)
    snr_straight = SNR_3(apo_sl_wl,apo_sl_flux,boundaries='flat_continuum',rebin=False,separate=False)
    return snr_ha,snr_straight






def snr_ha(datafile_class_file,return_only_snr=False,boundaries=[6549, 6550.7, 6576.0, 6578.0]):
    file=datafile_class_file
    [a, b, c, d] = boundaries
    wave = file.line6562_order_rebin.wl
    flux = file.line6562_order_rebin.flux
    normwave = np.hstack((wave[(wave > a) & (wave < b)], wave[(wave > c) & (wave < d)]))
    normflux = np.hstack((flux[(wave > a) & (wave < b)], flux[(wave > c) & (wave < d)]))
    # fit line trough slice
    # print('ss', startwl, endwl)
    slope, height = np.polyfit(normwave, normflux, 1)
    # print 'slope and height are', slope, height
    fit = np.poly1d([slope, height])
    nnf = []
    for k, nwl in enumerate(normwave):
        nnf.append(normflux[k] / fit(nwl))
    snr_ha = 1 / np.std(nnf)
    snr_per_60 = snr_ha * 60 / file.exptime
    count_per_60 = (snr_ha ** 2) * 60 / file.exptime
    if return_only_snr is True:
        return snr_ha
    else:
        return snr_ha,snr_per_60,count_per_60

def rebin_spec(wl, flux, wavnew):
    wave = np.array(wl)
    flux = np.array(flux)
    spec = spectrum.ArraySourceSpectrum(wave=wave, flux=flux)
    f = np.ones(len(wave))
    filt = spectrum.ArraySpectralElement(wave, f, waveunits='angstrom')
    obs = observation.Observation(spec, filt, binset=wavnew, force='taper')

    return obs.binflux


def rebin(wlarray, fluxarray):
    # print wlarray[0],wlarray[-1]
    spl = interpolate.InterpolatedUnivariateSpline(wlarray, fluxarray)
    xs = np.arange(wlarray[0], wlarray[-1], 0.05)
    return xs, spl(xs)

def rebin2(wlarray, fluxarray,step=None):
    # print wlarray[0],wlarray[-1]
    # spl = interpolate.InterpolatedUnivariateSpline(wlarray, fluxarray)
    # xs = np.arange(pars[0], pars[1], pars[2])
    if step==None:
        dx = 0.1
    else:
        dx=step
    wlarray2,fluxarray2 = remove_nan(wlarray,fluxarray)
    xs = np.arange(wlarray2[0], wlarray2[-1], dx)
    ys = rebin_spec(wlarray2,fluxarray2,xs)
    return xs,ys

def remove_nan(wl,flux):
    flux2=[]
    wl2=[]
    for i, item in enumerate(flux):
        if np.isnan(item) == False:
            wl2.append(wl[i])
            flux2.append(item)
    return wl2,flux2

def slice_spec(wl,flux,wl_start,wl_stop):
    slice = flux[(wl>wl_start) & (wl<wl_stop)]
    wlslice = wl[(wl>wl_start) & (wl<wl_stop)]
    return wlslice,slice

def quotient(lf1, lf2):
    # wl_rebin1,flux_rebin1 = reduce_spectrum(file1,18.5)
    # wl_rebin2,flux_rebin2 = reduce_spectrum(file2,18.5)
    # wl1, f1 = extractdata(line[1],file1)
    # cor_v1 = barcor(file1)[0]-18.5
    # wl_cor1 = wlcorrection(wl1,correction_velocity=cor_v1)
    # wl_rebin1,flux_rebin1 = rebin(wl_cor1,f1)
    # lw1, lf1,nf = normalize(wl_rebin1,flux_rebin1,line[3],line[4],line[5],line[6],startwl,endwl)
    # # wl2, f2 = extractdata(line[1],file2)
    # # cor_v2 = barcor(file2)[0]+18.5
    # # wl_cor2 = wlcorrection(wl2,correction_velocity=cor_v2)
    # # wl_rebin2,flux_rebin2 = rebin(wl_cor2,f2)
    # lw2, lf2, nf = normalize(wl_rebin2,flux_rebin2,line[3],line[4],line[5],line[6],startwl,endwl)
    quotientflux = []
    for i in range(len(lf1)):
        qf = lf2[i]/lf1[i]
        quotientflux.append(qf)
    return quotientflux

def quotient_eShel(file1,file2,line,startwl,endwl,v_rad=18.5):
    wl_rebin1, flux_rebin1 = reduce_spectrum(file1, radial_velocity=v_rad)
    lw1, lf1, nf = normalize(wl_rebin1, flux_rebin1, line[3], line[4], line[5], line[6], startwl, endwl)
    v, vsini = wl_to_velocity(lw1, line[2])
    wl_rebin2, flux_rebin2 = reduce_spectrum(file2, radial_velocity=v_rad)
    lw2, lf2, nf = normalize(wl_rebin2, flux_rebin2, line[3], line[4], line[5], line[6], startwl, endwl)
    qf = quotient(lf1,lf2)
    return lw1,v, qf

def velocity_to_wl(v_list, linecenter):
    wl=[]
    for v in v_list:
        l=((v/299792.5)+1)*linecenter
        wl.append(l)
    return np.array(wl)

def wl_to_velocity(wavelengths, linecenter):
    vsini = 127
    # vrad = 18.5
    v = []
    # BCCor = barcor(file)
    for l in wavelengths:
        velo = (299792.5*(l/linecenter - 1))
        v.append(velo)
    return np.array(v), vsini


def wl_to_velocity_2(wavelengths, linecenter, offset_v):
    vsini = 127
    # vrad = 18.5
    v = []
    # BCCor = barcor(file)
    for l in wavelengths:
        velo = (299792.5*(l/linecenter - 1))+offset_v
        v.append(velo)
    return np.array(v)


def wlcorrection(wlarray,correction_velocity=18.5):
    c = 299792.458
    v = correction_velocity
    wl2 = []
    for wl in wlarray:
        wl2.append(wl*(1+(v/c)))
        # wl2.append(math.exp(logwl2))
    return wl2

def reduce_spectrum(file,radial_velocity=0):
    wl,flux = extractdata(35,file)
    BCCor = barcor(file)[0]
    wl2 = wlcorrection(wl,BCCor-radial_velocity)
    # return wl2,flux
    wlrebin, fluxrebin = rebin(wl2, flux)
    return wlrebin,fluxrebin

def reduce_spectrum_lapalma(file, v_rad=18.5):
    datafile = pf.open(file)
    header = datafile[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    flux = datafile[0].data
    wl = np.exp(np.arange(naxis1) * cdelt1 + crval1 - v_rad / 299792.458)
    wl_rebin, flux_rebin = rebin2(wl, flux)
    datafile.close()
    return wl_rebin, flux_rebin


def reduce_spectrum2(file,radial_velocity=0):
    wl,flux = extractdata(35,file)
    BCCor = barcor(file)[0]
    wl2 = wlcorrection(wl,BCCor-radial_velocity)
    # return wl2,flux
    wlrebin, fluxrebin = rebin2(wl2, flux)
    return wlrebin,fluxrebin

def TVS(filelist,line,startwl,endwl,v_rad=18.5):
    lws = []
    lfs = []
    std_exp_list = []
    std_exp_weights = []
    for file in filelist:
        wl_rebin,flux_rebin = reduce_spectrum2(file,radial_velocity=v_rad)
        # lw, lf, nf = normalize(wl_rebin,flux_rebin,line[5],line[5]+((line[6]-line[5])/2),line[5]+((line[6]-line[5])/2),line[6],startwl,endwl)
        lw, lf, nf = normalize(wl_rebin,flux_rebin,line[3],line[4],line[5],line[6],startwl,endwl)
        v,vsini = wl_to_velocity(lw, line[2])
        lws.append(lw)
        lfs.append(lf)
        s = np.std(nf)
        std_exp_list.append(s)
        std_exp_weights.append(1/s**2)
    std_expected = np.average(std_exp_list, weights=std_exp_weights)
    TVS  = []
    for datapoint in np.transpose(lfs):
        TVS.append(np.std(datapoint)/std_expected)

    return lws[0],TVS,v,len(filelist)


def TVS_LaPalma(filelist,line,startwl,endwl,v_rad=18.5):
    lws = []
    lfs = []
    std_exp_list = []
    std_exp_weights = []
    for file in filelist:
        datafile = pf.open(file)
        header = datafile[0].header
        naxis1 = header['NAXIS1']
        crval1 = header['CRVAL1']
        cdelt1 = header['CDELT1']
        flux = datafile[0].data
        wl = np.exp(np.arange(naxis1)*cdelt1 + crval1 - v_rad/299792.458)
        wl_rebin, flux_rebin = rebin2(wl,flux)
        lw, lf, nf = normalize(wl_rebin,flux_rebin,line[2],line[3],line[4],line[5],startwl,endwl)
        v, vsini = wl_to_velocity(lw, line[1])
        lws.append(lw)
        lfs.append(lf)
        s = np.std(nf)
        std_exp_list.append(s)
        std_exp_weights.append(1/s**2)
        datafile.close()
    std_expected = np.average(std_exp_list, weights=std_exp_weights)
    TVS  = []
    for datapoint in np.transpose(lfs):
        TVS.append(np.std(datapoint)/std_expected)

    return lws[0],TVS,v,len(filelist)

def TVS_masterfiles(filelist,line):
    wls = []
    lfs = []
    std_exp_list = []
    std_exp_weights = []
    for file in filelist:
        linedata = getattr(file, line)
        flux = linedata.flux
        wl = linedata.wl
        v = linedata.v_cor
        v2=linedata.v
        if (v2==v).all:
            print('good')
        else:
            print('bad')
        nf = linedata.normalizationflux
        wls.append(wl)
        lfs.append(flux)
        s = np.std(nf)
        std_exp_list.append(s)
        std_exp_weights.append(1/s**2)
    std_expected = np.average(std_exp_list, weights=std_exp_weights)
    TVS  = []
    # print(len( np.transpose(lfs)))
    for datapoint in np.transpose(lfs):
        TVS.append(np.std(datapoint)/std_expected)
    return np.array(wl),np.array(TVS),np.array(v),len(filelist)

def TVS_masterfiles_order(filelist,line):
    wls = []
    lfs = []
    std_exp_list = []
    std_exp_weights = []
    wl_file_1 = getattr(filelist[0], line).wl
    wavenew = np.arange(wl_file_1[10], wl_file_1[-10], 0.05)
    lineinfo = getattr(filelist[0], line).lineinfo
    linecenter = lineinfo[1]
    vrad = -18.5
    for file in filelist:
        offset_v = file.baricentric_correction + vrad
        linedata = getattr(file, line)
        v = linedata.v_cor
        v2 = linedata.v
        if (v2 == v).all():
            print('good')
        else:
            print('bad')

        flux = linedata.flux
        wl = linedata.wl
        flux_rebin = rebin_spec(wl,flux,wavenew)
        vnew=wl_to_velocity(wavenew,linecenter)[0]
        nf = linedata.normalizationflux
        wls.append(wavenew)
        lfs.append(flux_rebin)
        s = np.std(nf)
        std_exp_list.append(s)
        std_exp_weights.append(1/s**2)
    std_expected = np.average(std_exp_list, weights=std_exp_weights)
    TVS  = []
    # print(len( np.transpose(lfs)))
    for datapoint in np.transpose(lfs):
        TVS.append(np.std(datapoint)/std_expected)
    return np.array(wavenew),np.array(TVS),np.array(vnew),len(filelist)

def average_masterfiles(filelist,line):
    wls = []
    lfs = []
    for file in filelist:
        linedata = getattr(file, line)
        flux = linedata.flux
        wl = linedata.wl
        v = linedata.v_cor
        nf = linedata.normalizationflux
        wls.append(wl)
        lfs.append(flux)
    average_lineprofile = []
    for datapoint in np.transpose(lfs):
        average_lineprofile.append(np.average(datapoint))
    return np.array(wl),np.array(average_lineprofile),np.array(v),len(filelist)



def overplot_LaPalma(filelist,line,startwl,endwl,separate_lines=False, v_rad = 18.5):
    vsini = 127
    lfs = []
    vs = []
    a = 0.0
    if separate_lines:
        ap = 0.1
    elif not separate_lines:
        ap = 0.0
    else:
        raise SyntaxError('separate_lines needs to be Boolean')
    for file in filelist:
        datafile = pf.open(file)
        header = datafile[0].header
        naxis1 = header['NAXIS1']
        crval1 = header['CRVAL1']
        cdelt1 = header['CDELT1']
        flux = datafile[0].data
        wl = np.exp(np.arange(naxis1) * cdelt1 + crval1 - v_rad / 299792.458)
        lw, lf, nf = normalize(wl, flux, line[2], line[3], line[4], line[5], startwl, endwl)
        v, vsini = wl_to_velocity(lw, line[1])
        lfs.append(np.array(lf) + a)
        vs.append(v)
        a+=ap
        datafile.close()
    return vs, lfs

def overplot_masterfiles(filelist,line,separate_lines=False):
    vsini = 127
    lfs = []
    vs = []
    a = 0.0
    if separate_lines:
        ap = 0.1
    elif not separate_lines:
        ap = 0.0
    else:
        raise SyntaxError('separate_lines needs to be Boolean')
    for file in filelist:
        linedata = getattr(file, line)
        flux = linedata.flux
        wl = linedata.wl

        v = linedata.v_cor
        lfs.append(np.array(flux) + a)
        vs.append(v)
        a += ap
    return vs, lfs

def line_separation(bl):
    if bl is True:
        return 0.1
    elif bl is False:
        return 0.0
    else:
        raise SyntaxError('separate_lines needs to be Boolean')
def overplot_masterfiles_order(filelist,line,separate_lines=False,return_wl=False,manual_rebin=True,rebin_size=0.1):
    vsini = 127
    lfs = []
    vs = []
    wls = []
    a = 0.0
    wl_file_1 = getattr(filelist[0], line).wl
    wavenew = np.arange(wl_file_1[10], wl_file_1[-10], rebin_size)
    lineinfo = getattr(filelist[0], line).lineinfo
    linecenter = lineinfo[1]
    # vrad = -18.5
    # offset_v = filelist[0].baricentric_correction + vrad
    # v = wl_to_velocity_2(wavenew, linecenter, offset_v)
    ap=line_separation(separate_lines)
    for file in filelist:
        # offset_v = file.baricentric_correction + vrad
        v_new = wl_to_velocity(wavenew, linecenter)[0]
        linedata = getattr(file, line)
        v = linedata.v_cor
        v2 = linedata.v
        flux = linedata.flux
        wl = linedata.wl
        flux_rebin = rebin_spec(wl, flux, wavenew)
        # print(flux_rebin.shape,v_new.shape)
        if manual_rebin is True:
            lfs.append(np.array(flux_rebin) + a)
            vs.append(v_new)
            wls.append(wavenew)
        elif manual_rebin is False:
            lfs.append(flux)
            vs.append(v)
            wls.append(wl)
        a += ap
    if return_wl is True:
        return wls,vs,lfs
    else:
        return vs, lfs

# def plot_TVS(filelist, linelist ):
#     for line in linelist:
#         lws,TVS = TVS(filelist,line,line[3],line[6])
#         v,vsini = wl_to_velocity(lws, line[2], filelist[1])
#         plt.title(line[0]+str(line[2]))
#         plt.axvline(x=-vsini)
#         plt.axvline(x=vsini)
#         plt.plot(v,TVS)
#         plt.show()
#         plt.close()
def overplot2(filelist,lapalmafilelist,line,lapalmaline,v_rad,startwl,endwl,together=True):

    vsini =127
    lfs = []
    vs= []
    a = 0.0
    for file in filelist:
        if together == False:
            a+=0.1
        wl_rebin,flux_rebin = reduce_spectrum2(file,radial_velocity=v_rad)
        # lw, lf, nf = normalize(wl_rebin,flux_rebin,line[5],line[5]+((line[6]-line[5])/2),line[5]+((line[6]-line[5])/2),line[6],startwl,endwl)
        lw, lf, nf = normalize(wl_rebin,flux_rebin,line[3],line[4],line[5],line[6],startwl,endwl)
        v,vsini = wl_to_velocity(lw, line[2])
        lfs.append(np.array(lf)+a)
        vs.append(v)
    return vs,lfs

def overplot(filelist,lapalmafilelist,line,lapalmaline,v_rad,startwl,endwl,together=False):
    vsini =127
    lfs = []
    vs= []
    for file in filelist:
        # print 'asdf',file
        wl_rebin,flux_rebin = reduce_spectrum2(file,radial_velocity=v_rad)
        # lw, lf, nf = normalize(wl_rebin,flux_rebin,line[5],line[5]+((line[6]-line[5])/2),line[5]+((line[6]-line[5])/2),line[6],startwl,endwl)
        lw, lf, nf = normalize(wl_rebin,flux_rebin,line[3],line[4],line[5],line[6],startwl,endwl)
        v,vsini = wl_to_velocity(lw, line[2])
        lfs.append(lf)
        vs.append(v)
    return vs,lfs
        # plt.plot(v,lf)
    # plt.title(line[0]+' '+str(round(line[2],1)))
    # plt.xlabel('Velocity (km/s)')
    # plt.ylabel('Normalized Flux')
    # plt.axvline(x=-vsini)
    # plt.axvline(x=vsini)
    # plt.xlim(-500,500)
    # if together==False:
    #     plt.show()
    #     plt.close()
    # for file in lapalmafilelist:
    #     datafile = pf.open(file)
    #     header = datafile[0].header
    #     naxis1 = header['NAXIS1']
    #     crval1 = header['CRVAL1']
    #     cdelt1 = header['CDELT1']
    #     flux = datafile[0].data
    #     wl = np.exp(np.arange(naxis1)*cdelt1 + crval1 - v_rad/299792.458)
    #     lw, lf, nf = normalize(wl,flux,lapalmaline[2],lapalmaline[3],lapalmaline[4],lapalmaline[5],startwl=lapalmaline[2],endwl=lapalmaline[5])
    #     v,vsini = wl_to_velocity(lw, lapalmaline[1])
    #     plt.plot(v,lf,color='blue')
    # if together == False:
    #     plt.title(line[0]+' '+str(round(line[2],1)))
    #     plt.xlabel('Velocity (km/s)')
    #     plt.ylabel('Normalized Flux')
    #     plt.xlim(-500,500)
    #     plt.axvline(x=-vsini)
    #     plt.axvline(x=vsini)
    # plt.show()
    # plt.close()

# def equivalent_width(wl,flux,linecenter):
#     # print wl, linecenter
#     velo,vsini = wl_to_velocity(wl,linecenter)
#     # print vsini
#     ex=50
#     a = -vsini-ex
#     b = vsini+ex
#     v_linepart = velo[(velo>a)&(velo<b)]
#     # print velo
#     flux_linepart = flux[(velo>a)&(velo<b)]
#     # print flux_linepart
#     ew =np.sum(flux_linepart)
#     return v_linepart,flux_linepart,ew

# def equivalent_width(velo,wl,flux,linecenter,snr):
#     # print wl, linecenter
#
#     # print vsini
#     vlim=200
#     if 6560>linecenter>6565:
#         vlim = 500
#     wl_linepart = wl[(velo>-vlim)&(velo<vlim)]
#     dwl = wl_linepart[-1]-wl_linepart[0]
#     v_linepart = velo[(velo>-vlim)&(velo<vlim)]
#     # print velo
#     flux_linepart = flux[(velo>-vlim)&(velo<vlim)]
#     # print flux_linepart
#     F_avg = np.average(flux_linepart)
#     # print F_avg
#     ew =dwl*(1-F_avg)
#     # print F_avg,dwl,ew,snr
#     er = np.sqrt(1+(1/F_avg))*(dwl-ew)/snr
#     # print er
#     return er, ew


def equivalent_width(velo,wl,flux,linecenter,snr,vlim=[-500,500],ha_vlim=[-380,620]):
    if 6560 <= linecenter <= 6565:
        ll = ha_vlim[0]
        ul = ha_vlim[1]
    else:
        ll=vlim[0]
        ul=vlim[1]
    # print(linecenter,vlim)
    wl_linepart = wl[(velo > ll) & (velo < ul)]
    dwl = wl_linepart[-1] - wl_linepart[0]
    v_linepart = velo[(velo > ll) & (velo < ul)]
    # print velo
    flux_linepart = flux[(velo > ll) & (velo < ul)]
    # print flux_linepart
    F_avg = np.average(flux_linepart)
    # print F_avg
    ew = dwl * (1 - F_avg)
    # print F_avg,dwl,ew,snr
    er = np.sqrt(1 + (1 / F_avg)) * (dwl - ew) / snr
    # print er
    return er, ew



def t_ew(velo,wl,flux,linecenter,snr,vlim=200):
    ll=-vlim
    ul = vlim
    if 6560 <= linecenter <=6565:
        vlim = 800
        ll = -380
        ul = 620
    # print(linecenter,vlim)
    wl_linepart = wl[(velo > ll) & (velo < ul)]
    dwl = wl_linepart[-1] - wl_linepart[0]
    v_linepart = velo[(velo > ll) & (velo < ul)]
    # print velo
    flux_linepart = flux[(velo > ll) & (velo < ul)]
    # print flux_linepart
    F_avg = np.average(flux_linepart)
    # print F_avg
    ew = dwl * (1 - F_avg)
    # print F_avg,dwl,ew,snr
    er = np.sqrt(1 + (1 / F_avg)) * (dwl - ew) / snr
    # print er
    return er, ew, wl_linepart,v_linepart,flux_linepart

def ls_periodogram(times,ews, searchrange=[2,20]):
    rotational_periods = np.linspace(searchrange[0],searchrange[1],1000)
    frequency_array =(2*np.pi)/rotational_periods

    ls = scisig.lombscargle(times,ews,frequency_array)
    return rotational_periods, ls

def double_line(x,x0,a1,b1,tau1,a2,b2,tau2):
    # return a1*np.exp(-tau1*np.exp(-((x-x0)/2*b1)**2))+a2*np.exp(-tau2*np.exp(-((x-x0-14.332)/2*b2)**2))+c
    return a1*np.exp(-tau1*np.exp(-((x-x0)/2*b1)**2))+a2*np.exp(-tau2*np.exp(-((x-x0+5.97)/2*b2)**2))

def fitfraun(file):
    apo_wl,apo_flux = extractdata(35,file)
    # print len(apo_wl), len(apo_flux)
    dat_x = apo_wl[(apo_wl>5880)&(apo_wl<5905)]
    dat_y = apo_flux[(apo_wl>5880)&(apo_wl<5905)]
    parms, pcov = curve_fit(double_line,dat_x,dat_y,p0=(5896.92,0.55,2.5,0.4,0.55,2.5,0.4))
    # wl_shift.append(parms[0])
    # x = np.linspace(5880,5905,200)
    # y = double_line(x,parms[0],parms[1],parms[2],parms[3],parms[4],parms[5],parms[6])

    return parms[0]

def fitfraunlp(file):

    datafile = pf.open(file)
    header = datafile[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    flux = datafile[0].data
    wl = np.exp(np.arange(naxis1)*cdelt1 + crval1 )
    apo_wl, apo_flux, nnf, lineflux,fit= normalize(np.array(wl),np.array(flux),5888.5,5889,5894.2,5894.75,5800,6000)
    # print len(apo_wl), len(apo_flux)
    dat_x = apo_wl[(apo_wl>5880)&(apo_wl<5905)]
    dat_y = apo_flux[(apo_wl>5880)&(apo_wl<5905)]
    parms, pcov = curve_fit(double_line,dat_x,dat_y,p0=(5896.92,0.55,2.5,0.4,0.55,2.5,0.4))
    # wl_shift.append(parms[0])
    # x = np.linspace(5880,5905,200)
    # y = double_line(x,parms[0],parms[1],parms[2],parms[3],parms[4],parms[5],parms[6])
    # plt.plot(dat_x,dat_y)
    # plt.plot(x,y)
    # plt.show()
    # plt.close()
    datafile.close()
    return parms[0]

def fitfraun_demetra(file):

    datafile = pf.open(file)
    header = datafile[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    flux = datafile[0].data
    wl = np.arange(naxis1)*cdelt1 + crval1
    apo_wl, apo_flux, nnf,lineflux, fit = normalize(np.array(wl),np.array(flux),5888.5,5889,5894.2,5894.75,5800,6000)
    # print len(apo_wl), len(apo_flux)
    dat_x = apo_wl[(apo_wl>5880)&(apo_wl<5905)]
    dat_y = apo_flux[(apo_wl>5880)&(apo_wl<5905)]
    parms, pcov = curve_fit(double_line,dat_x,dat_y,p0=(5896.92,0.55,2.5,0.4,0.55,2.5,0.4))
    # wl_shift.append(parms[0])
    # x = np.linspace(5880,5905,200)
    # y = double_line(x,parms[0],parms[1],parms[2],parms[3],parms[4],parms[5],parms[6])
    # plt.plot(dat_x,dat_y)
    # plt.plot(x,y)
    # plt.show()
    # plt.close()
    datafile.close()
    return parms[0]

# def ew(linelist,filelist):
#     vsini =127
#     ews_avg = []
#     for file in filelist:
#         ews = []
#         print file
#         for line in linelist:
#             wl_rebin,flux_rebin = reduce_spectrum(file,radial_velocity=-18.5)
#             lw, lf, nf = normalize(wl_rebin,flux_rebin,line[3],line[4],line[5],line[6],startwl=line[3],endwl=line[6])
#             # print lf
#             # print lf
#             _,__,ew = equivalent_width(lw,lf,line[2])
#             # print ew
#             ews.append(ew)
#         ews_avg.append(np.average(ews))
#     return ews_avg
#


# def sl(x,x0,tau,w,A,A2,w2,tau2):
#     return A*np.exp(-tau*np.exp(-((x-x0)/2*w)**2))+A2*np.exp(-tau2*np.exp(-((x-x0)/2*w2)**2))



def normalize_fluxarray(flux,wl=None,SG = False):
    if wl is None:
        wl = range(len(flux))
    fitparams = np.polyfit(wl,flux,1)
    fitparams2 = np.polyfit(wl,flux,3)
    x1 = fitparams[0]
    x0 = fitparams[1]
    fitted_flux = np.array(wl)*x1 +x0
    # fitted_flux_2 = np.polyval(fitparams2,wl)
    fit_SG_flux = savitzky_golay(flux,25,4)
    # normalized_flux = flux/fitted_flux
    if SG == True:
        # plt.plot(flux)
        # plt.plot(fitted_flux_2)
        # plt.show()
        # plt.close()
        normalized_flux = flux/fit_SG_flux
        # normalized_flux = flux/fitted_flux_2
    elif SG == False:
        normalized_flux = flux/fitted_flux
    else:
        raise TypeError('SG needs to be Boolean')
    return normalized_flux

def real_snr(file,linecenter,continuumpiece = None,SG=False):
    linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, [6615.0, 6618.3]], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, [4877.0, 4883.0]], ['Hy', 4340.472, 4322, 4324, 4357, 4360, [4298.0, 4302.0]], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, [4031.0, 4038.0]], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, [4458.0, 4463.0]], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, [4734.0, 4740.0]], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, [5819.0, 5833.0]], ['He_II', 4541.6, 4498, 4499, 4580, 4581, [4579.2, 4587.6]], ['He_II', 4685.804, 4679, 4680, 4690, 4691, [4718.0, 4740.0]], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, [5370.0, 5380.0]], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, [5564.3, 5575.0]], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, [5831.2, 5833.1]]]
    sl = sorted(linelist, key=lambda x: abs(x[1]-linecenter))
    # print sl


    # print type(wl),type(flux)
    if continuumpiece==None:
        closest_line = sl[0]
        start = closest_line[6][0]
        stop = closest_line[6][1]

        dif = np.min(np.absolute(np.array([start,stop])-linecenter))
        if dif > 100:
            warnings.warn('Warning: SNR calculated from continuum piece '+str(dif)+' Angstrom from linecenter, SNR might be inaccurate')
        # print start,stop
    else:
        start = continuumpiece[0]
        stop = continuumpiece[1]

    datafile = pf.open(file)
    header = datafile[0].header
    naxis1 = header['NAXIS1']
    crval1 = header['CRVAL1']
    cdelt1 = header['CDELT1']
    snroud = header['SNR50']
    rawflux = np.array(datafile[0].data)
    rawwl = np.exp(np.arange(naxis1)*cdelt1 + crval1 - 18.5/299792.458)
    wl, flux = np.array(remove_nan(rawwl,rawflux))
    continuumpiece  = flux[(wl>start)&(wl<stop)]


    # wl_piece = wl[(wl>start)&(wl<stop)]
    # print '-----------'
    # print np.shape(wl_piece),np.shape(flux_piece)
    # print '-------------'
    snr_poisson = np.average(continuumpiece) / np.sqrt(np.average(continuumpiece))
    flux_piece_norm = normalize_fluxarray(continuumpiece,SG=SG)
    snr_real = np.average(flux_piece_norm)/np.std(flux_piece_norm)
    datafile.close()
    return snr_real,snr_poisson


def EW_stats(dat):
    ew = dat[0]
    print(ew)
    phase = dat[1]
    error = dat[2]
    # error = np.std(ew)
    # lineused = dat[3][0]
    # print lineused
    # linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
    p1=[0.005,0.5, 1]
    fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
    p2 = [0,1]
    fit2 = curve_fit(my_line, phase, ew, p0=p2,  sigma=error, absolute_sigma=True)
    p3 =[1]
    fit3 = curve_fit(flat_line , phase, ew, p0=p3,  sigma=error, absolute_sigma=True)
    # t= np.linspace(0,1,50)
    data_fit_sin = my_sin(np.array(phase), *fit1[0])
    data_fit_line = my_line(np.array(phase), *fit2[0])
    data_fit_flat = flat_line(np.array(phase), *fit3[0])

    residuals = (ew - data_fit_sin)/error
    chisq = np.sum(residuals**2) / (len(ew)-1-3)
    red_chisq = redchisqg(ew,data_fit_sin,deg=3)

    k_1 = 3
    res_1 = (ew-data_fit_sin)
    SSR_1 = np.sum(res_1**2)
    N_1 = len(ew)
    s2_1 = SSR_1 / N_1
    L_1 = ( 1.0/np.sqrt(2*np.pi*s2_1) ) ** N_1 * np.exp( -SSR_1/(s2_1*2.0) )

    k_2 = 2
    res_2 = (ew-data_fit_line)
    SSR_2 = np.sum(res_2**2)
    N_2 = len(ew)
    s2_2 = SSR_2 / N_2
    L_2 = ( 1.0/np.sqrt(2*np.pi*s2_2) ) ** N_2 * np.exp( -SSR_2/(s2_2*2.0) )

    k_3 = 1
    res_3 = (ew-data_fit_flat)
    SSR_3 = np.sum(res_3**2)
    N_3 = len(ew)
    s2_3 = SSR_3 / N_3
    L_3 = ( 1.0/np.sqrt(2*np.pi*s2_3) ) ** N_3 * np.exp( -SSR_3/(s2_3*2.0) )

    AIC_1 = 2*k_1 - 2*np.log(L_1)
    AIC_2 = 2*k_2 - 2*np.log(L_2)
    AIC_3 = 2*k_3 - 2*np.log(L_3)
    probfactor = np.exp((AIC_1-AIC_2)/2)

    # print 'ln(L) =', np.log( L )

    # print chisq, chisq2, chisq-chisq2
    chi2_sin,p_sin = ss.chisquare(ew,data_fit_sin,ddof=3)
    chi2_line,p_line = ss.chisquare(ew,data_fit_line,ddof=2)
    chi2_flat,p_flat = ss.chisquare(ew,data_fit_flat,ddof=1)
    # dif1 = p_sin/p_line
    # dif2 = p_sin/p_flat
    # line_to_write = linename+' \t '+ str(chisq2)+ '\t' +str(AIC_1) + '\t' +str(AIC_2) + '\t' +str(AIC_3) + '\n'
    # print line_to_write
    # f.write(line_to_write)
    # print '#######'
    # print linename
    # print chi2_sin-chi2_line
    # print dif1
    # print '#######'
    return chisq, red_chisq, AIC_1,AIC_2,AIC_3,probfactor

def EW_stats2(ew, phase, error):
    # error = np.std(ew)
    # lineused = dat[3][0]
    # print lineused
    # linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
    p1=[0.005,0.5, 1]
    fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
    p2 = [0,1]
    fit2 = curve_fit(my_line, phase, ew, p0=p2,  sigma=error, absolute_sigma=True)
    p3 =[1]
    fit3 = curve_fit(flat_line , phase, ew, p0=p3,  sigma=error, absolute_sigma=True)
    # t= np.linspace(0,1,50)
    data_fit_sin = my_sin(np.array(phase), *fit1[0])
    data_fit_line = my_line(np.array(phase), *fit2[0])
    data_fit_flat = flat_line(np.array(phase), *fit3[0])

    residuals = (ew - data_fit_sin)/error
    chisq = np.sum(residuals**2) / (len(ew)-1-3)
    red_chisq = redchisqg(ew,data_fit_sin,deg=3)

    k_1 = 3
    res_1 = (ew-data_fit_sin)
    SSR_1 = np.sum(res_1**2)
    N_1 = len(ew)
    s2_1 = SSR_1 / N_1
    L_1 = ( 1.0/np.sqrt(2*np.pi*s2_1) ) ** N_1 * np.exp( -SSR_1/(s2_1*2.0) )

    k_2 = 2
    res_2 = (ew-data_fit_line)
    SSR_2 = np.sum(res_2**2)
    N_2 = len(ew)
    s2_2 = SSR_2 / N_2
    L_2 = ( 1.0/np.sqrt(2*np.pi*s2_2) ) ** N_2 * np.exp( -SSR_2/(s2_2*2.0) )

    k_3 = 1
    res_3 = (ew-data_fit_flat)
    SSR_3 = np.sum(res_3**2)
    N_3 = len(ew)
    s2_3 = SSR_3 / N_3
    L_3 = ( 1.0/np.sqrt(2*np.pi*s2_3) ) ** N_3 * np.exp( -SSR_3/(s2_3*2.0) )

    AIC_1 = 2*k_1 - 2*np.log(L_1)
    AIC_2 = 2*k_2 - 2*np.log(L_2)
    AIC_3 = 2*k_3 - 2*np.log(L_3)
    probfactor = np.exp((AIC_1-AIC_3)/2)

    # print 'ln(L) =', np.log( L )

    # print chisq, chisq2, chisq-chisq2
    # chi2_sin,p_sin = ss.chisquare(ew,data_fit_sin,ddof=3)
    # chi2_line,p_line = ss.chisquare(ew,data_fit_line,ddof=2)
    # chi2_flat,p_flat = ss.chisquare(ew,data_fit_flat,ddof=1)
    # dif1 = p_sin/p_line
    # dif2 = p_sin/p_flat
    # line_to_write = linename+' \t '+ str(chisq2)+ '\t' +str(AIC_1) + '\t' +str(AIC_2) + '\t' +str(AIC_3) + '\n'
    # print line_to_write
    # f.write(line_to_write)
    # print '#######'
    # print linename
    # print chi2_sin-chi2_line
    # print dif1
    # print '#######'
    return chisq, red_chisq, AIC_1,AIC_2,AIC_3,probfactor

# def general_TVS(filelist, wl_start,wl_stop,cont_wl_start, cont_wl_stop):
def calculate_likelihood(data, fitdata,error):
    # print type(data), type(fitdata), type(error)
    p_array = 1/np.sqrt(2*np.pi*(error**2))*np.exp(-((data-fitdata)**2)/(2*(error**2)))
    likelihood = np.prod(p_array)
    lp_array =  1/np.sqrt(2*np.pi*(error**2))*-((data-fitdata)**2)/(2*(error**2))

    return likelihood

def calculate_log_likelihood(data,fitdata,error):
    llh_array = -0.5*(np.log(2*np.pi*error**2)+(data-fitdata)**2/(error**2))
    llh = np.sum(llh_array)
    return llh
def EW_stats3(ew, phase, error,fx2='off'):
    # error = np.std(ew)
    # lineused = dat[3][0]
    # print lineused
    # linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
    p1=[0.005,0.5, 1]
    if fx2 == 'om':
        fit1 = curve_fit(my_sin2, phase, ew, p0=p1, sigma=error, absolute_sigma=True)
    else:
        fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
    p2 = [0,1]
    fit2 = curve_fit(my_line, phase, ew, p0=p2,  sigma=error, absolute_sigma=True)
    p3 =[1]
    fit3 = curve_fit(flat_line , phase, ew, p0=p3,  sigma=error, absolute_sigma=True)
    # t= np.linspace(0,1,50)
    data_fit_sin = my_sin(np.array(phase), *fit1[0])
    data_fit_line = my_line(np.array(phase), *fit2[0])
    data_fit_flat = flat_line(np.array(phase), *fit3[0])

    residuals = (ew - data_fit_sin)/error
    chisq = np.sum(residuals**2) / (len(ew)-1-3)
    red_chisq = redchisqg(ew,data_fit_sin,deg=3)

    k_1 = 3
    res_1 = (ew-data_fit_sin)/error
    SSR_1 = np.sum(res_1**2)
    N_1 = len(ew)
    s2_1 = SSR_1 / N_1
    L_1 = calculate_likelihood(ew,data_fit_sin,error)
    llh_1 = calculate_log_likelihood(ew,data_fit_sin,error)

    k_2 = 2
    res_2 = (ew-data_fit_line)
    SSR_2 = np.sum(res_2**2)
    N_2 = len(ew)
    s2_2 = SSR_2 / N_2
    L_2 = calculate_likelihood(ew,data_fit_line,error)
    llh_2 = calculate_log_likelihood(ew, data_fit_line, error)
    k_3 = 1
    res_3 = (ew-data_fit_flat)
    SSR_3 = np.sum(res_3**2)
    N_3 = len(ew)
    s2_3 = SSR_3 / N_3
    L_3 = calculate_likelihood(ew,data_fit_flat,error)
    llh_3 = calculate_log_likelihood(ew, data_fit_flat, error)

    AIC_1 = 2*k_1 - 2*np.log(L_1)
    AIC_1_2 =2*k_1 - 2*llh_1
    # print ' ------'
    # print AIC_1-AIC_1_2
    AIC_2 = 2*k_2 - 2*np.log(L_2)
    AIC_2_2= 2 * k_2 - 2 * llh_2
    # print AIC_2 - AIC_2_2
    AIC_3 = 2*k_3 - 2*np.log(L_3)
    AIC_3_2 = 2 * k_3 - 2 * llh_3
    # print AIC_3 - AIC_3_2
    # print' -------'
    probfactor = np.exp((AIC_1-AIC_3)/2)

    # print 'ln(L) =', np.log( L )

    # print chisq, chisq2, chisq-chisq2
    chi2_sin,p_sin = ss.chisquare(ew,data_fit_sin,ddof=3)
    chi2_line,p_line = ss.chisquare(ew,data_fit_line,ddof=2)
    chi2_flat,p_flat = ss.chisquare(ew,data_fit_flat,ddof=1)
    # dif1 = p_sin/p_line
    # dif2 = p_sin/p_flat
    # line_to_write = linename+' \t '+ str(chisq2)+ '\t' +str(AIC_1) + '\t' +str(AIC_2) + '\t' +str(AIC_3) + '\n'
    # print line_to_write
    # f.write(line_to_write)
    # print '#######'
    # print linename
    # print chi2_sin-chi2_line
    # print dif1
    # print '#######'
    return chisq, red_chisq, AIC_1,AIC_2,AIC_3,probfactor

def TVS_significance_level(Nfiles, p):
    degfree = Nfiles-1
    siglvl = np.sqrt(ss.distributions.chi2.ppf(1-p, df=degfree)/degfree)
    return siglvl


def open_linelist(filepath):
    workfileresource = open(filepath, 'rb')
    list = pickle.load(workfileresource)
    workfileresource.close()
    return list


def make_linelist(list, filepath):
    workfileresource = open(filepath, 'wb')
    pickle.dump(list, workfileresource)
    workfileresource.close()

def find_order(wl , demetra_file,return_full_order_list=False):
    orders = demetra_file.orders
    centerwl=np.average(wl)
    minwl=np.min(wl)
    maxwl=np.max(wl)
    ol = sorted(orders, key=lambda x: np.abs(x.wl_avg - centerwl))
    best_order = ol[0]
    if minwl<best_order.wl_start:
        warnings.warn('Minimum of WL array out of bounds of order, you cannot get this full interval from the same order')
    if maxwl>best_order.wl_end:
        warnings.warn('Maximum of WL array out of bounds of order, you cannot get this full interval from the same order')
    if return_full_order_list is True:
        return ol
    else:
        return best_order

def degrade_spectrum(wl,flux,spectral_resolution=10000, desired_snr=100,pre_rebin=False):
    deg_wl = wl

    deg_flux, fwhm = pyasl.instrBroadGaussFast(wl, flux, spectral_resolution,
                                         edgeHandling="firstlast", fullout=True, maxsig=5.0)
    if type(pre_rebin) is float:
        rebin_wl,loop_flux = rebin2(deg_wl,deg_flux,step=pre_rebin)
    elif pre_rebin is False or None:
        loop_flux = deg_flux
        rebin_wl = deg_wl
    else:
        print('Rebin turned off, to turn on, specify stepsize for wl array')
        loop_flux = deg_flux
        rebin_wl = deg_wl
    noise_array = []
    for f in loop_flux:
        noise_exp = f/desired_snr
        noise = np.random.normal(0, noise_exp)
        noise_array.append(noise)
    deg_flux_with_noise = np.array(loop_flux)+np.array(noise_array)
    return rebin_wl,deg_flux_with_noise

def degrade_spectrum_noise_first(wl,flux,spectral_resolution=10000, desired_snr=100,pre_rebin=0.05):
    noise_array = []
    for f in flux:
        noise_exp = f/desired_snr
        noise = np.random.normal(0, noise_exp)
        noise_array.append(noise)
    noisy_flux = np.array(flux)+noise_array
    if type(pre_rebin) == float or type(pre_rebin) == int:
        rebin_wl,rebin_noisy_flux = rebin2(wl,noisy_flux,step=pre_rebin)
    elif pre_rebin is False or pre_rebin is None:
        rebin_noisy_flux = noisy_flux
        rebin_wl, rebin_noisy_flux = rebin2(wl, noisy_flux, step=np.max(np.diff(wl)))
    else:
        print('incorrect import for pre_rebin parameter. Rebin set to 0.05, to customize either specify stepsize or set to None')
        rebin_wl, rebin_noisy_flux = rebin2(wl, noisy_flux, step=0.05)
    deg_noisy_flux, fwhm = pyasl.instrBroadGaussFast(rebin_wl, rebin_noisy_flux, spectral_resolution,
                                         edgeHandling="firstlast", fullout=True, maxsig=5.0)

    return rebin_wl,deg_noisy_flux

def make_frequency_array(start,stop,ss_start,ss_end,bigstep = 1/1000,smallstep=1/10000):
    #in days, will return 1/d frequency array
    a = np.arange(start, ss_start, bigstep)
    b = np.arange(ss_start, ss_end, smallstep)
    c = np.append(a, b)
    d = np.arange(ss_end, stop + bigstep, bigstep)
    e = np.append(c, d)
    return 1/e




