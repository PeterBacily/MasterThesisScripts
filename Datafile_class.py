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
import pickle
import os

class Line:
    def __init__(self, li,lc,wave,fl,velo,nf,vsini,snr,barcor,vrad,norm_boundaries):
        self.lineinfo = li
        self.normalization_boundaries_wl=norm_boundaries[0]
        self.normalization_boundaries_v = norm_boundaries[1]
        self.wl = wave
        self.v = velo
        self.v_cor = np.array(velo)+(barcor+vrad)
        self.flux = fl
        self.normalizationflux = nf
        self.vsini = vsini
        self.ew_error, self.ew = airmass.equivalent_width(self.v_cor,wave,fl,lc,snr)


def line_data(line,wl,flux,observatory,snr,bccor,vrad):
    if observatory == 'MERC':
        k = 1
        barcor = 0
    elif observatory == 'APO':
        k = 2
        barcor = bccor
    elif observatory == 'APO_DEMETRA':
        barcor = bccor
        k = 1
    else:
        print('observatory needs to be APO or MERC')
    center_wl = int(line[k])
    lw, lf, nf, _,_ = airmass.normalize(wl,flux,line[k+1],line[k+2],line[k+3],line[k+4],line[k+1]-20,line[k+4]+20)
    v, vsini = airmass.wl_to_velocity(lw, line[k])
    normalization_wl= [line[k+1],line[k+2],line[k+3],line[k+4]]
    normalization_v = airmass.wl_to_velocity(normalization_wl, line[k])
    return Line(line,line[k],lw,lf,v,nf,vsini,snr, barcor,vrad,[normalization_wl,normalization_v]),'line'+str(center_wl)

def linecenters(linelist,observatory):
    lcs = []
    if observatory == 'MERC':
        k=1
    elif observatory == 'APO':
        k=2
    elif observatory == 'APO_DEMETRA':
        k=1
    else:
        print('observatory needs to be APO or MERC or APO_DEMETRA')
    for line in linelist:
        lcs.append(int(line[k]))
    return lcs
# def save(object)



class Datafile_mercator:
    observatory = 'MERC'
    # linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, 'He I 4026'], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, 'He I 5876'], ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4542'], ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4686'], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5412'], ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]
    linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
                     ['Hb', 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'],
                     ['Hy', 4340.472, 4322, 4324, 4357, 4360, r'H$\gamma$ 4340'],
                     ['He_I', 4026.1914, 4016, 4020, 4032, 4036, 'He I 4026'],
                     ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'],
                     ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'],
                     ['He_I', 5875.621, 5863.0, 5864.5, 5884.6, 5885.5, 'He I 5875'],
                     ['He_II', 4541.6, 4498, 4499, 4580, 4581, 'He II 4541'],
                     ['He_II', 4685.804, 4679, 4680, 4690, 4691, 'He II 4685'],
                     ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, 'He II 5411'],
                     ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4921'],
                     ['He_I', 6678.15, 6656, 6660, 6690, 6695, r'He I 6678'],
                     ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, r'O III 5592'],
                     ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, r'C IV 5801']]

    def __init__(self, file,v_rad = 18.5,i='n/a'):
        fn = os.path.basename(file)
        data = pf.open(file)
        self.i =i
        self.filename = fn[:fn.rfind(".")]
        self.header = data[0].header
        self.time_and_date = airmass.timeanddate2(self.header['DATE-OBS'])
        self.original_filepath = file
        self.HJD = float(self.header['BJD'])
        self.phase =  airmass.aphase(self.header['BJD'])
        self.exptime = self.header['EXPTIME']
        self.altitude = float(self.header['TELALT'])
        self.airmass = 1/np.sin(2*np.pi*self.altitude/360)
        self.baricentric_correction = float(self.header['BVCOR'])
        self.fwl = airmass.fitfraunlp(file)
        self.velshift = 299792.458*(self.fwl-5895.92)/5895.92
        naxis1 = self.header['NAXIS1']
        crval1 = self.header['CRVAL1']
        cdelt1 = self.header['CDELT1']
        self.wl_original = np.exp(np.arange(naxis1) * cdelt1 + crval1 - v_rad / 299792.458)
        self.flux_original = data[0].data
        self.wl_rebin, self.flux_rebin = airmass.rebin2(self.wl_original,self.flux_original)
        self.available_lines = []
        self.snr_original =airmass.snr(self.wl_original,self.flux_original)
        self.snr = airmass.snr(self.wl_rebin,self.flux_rebin)
        for line in self.linelist:
            linedata,linekey = line_data(line,self.wl_rebin,self.flux_rebin,self.observatory,self.snr,0,0)
            linedata_original, lk = line_data(line,self.wl_original,self.flux_original,self.observatory,self.snr_original,0,0)
            setattr(self,linekey , linedata)
            setattr(self,linekey+'_original',linedata_original)
            self.available_lines.append(linekey)
        data.close()

class Datafile_apo:
    observatory = 'APO'
    # linelist = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
    linelist = [['Ha', 35, 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
     ['Hb', 35, 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'],
     ['He_I', 35, 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'],
     ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5875'],
     ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4541'],
     ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4685'],
     ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5411'],
     ['He_I', 35, 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'],
     ['He_I', 35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4921'],
     ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'],
     ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'],
     ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]

    def __init__(self, file,v_rad = 18.5,i='n/a',mark = 0):
        fn = os.path.basename(file)
        data = pf.open(file)
        self.original_filepath = file
        self.i =i
        self.mark = mark
        self.mark_explanation = '0 = no weird stuff,      1 = not usable due to very poor SNR,    2 = Not usable for EW, TVS, Quotient due to insufficient SNR,   3 =   Shows weird feature in Halpha'
        self.filename = fn[:fn.rfind(".")]
        self.header = data[0].header
        self.time_and_date = airmass.timeanddate2(self.header['DATE-OBS'])
        self.baricentric_correction, self.HJD = airmass.barcor(file)
        self.phase =  airmass.aphase(self.HJD)
        self.exptime = airmass.exposuretime(file)
        self.airmass, self.alt, JD = airmass.airmass(file)
        try:
            frwl = airmass.fitfraun(file)
        except RuntimeError:
            frwl = 5895.92
        self.fwl = frwl
        self.velshift = 299792.458*(self.fwl-5895.92)/5895.92
        self.wl_original, self.flux_original = airmass.extractdata(35, file)
        self.wl_rebin, self.flux_rebin = airmass.rebin2(self.wl_original,self.flux_original)
        self.available_lines = []
        self.snr_original =airmass.snr(self.wl_original,self.flux_original)
        self.snr = airmass.snr(self.wl_rebin,self.flux_rebin)
        for line in self.linelist:
            linedata,linekey = line_data(line,self.wl_rebin,self.flux_rebin,self.observatory,self.snr, self.baricentric_correction,-18.5)
            linedata_original, lk = line_data(line,self.wl_original,self.flux_original,self.observatory,self.snr_original,0,0)
            setattr(self,linekey , linedata)
            setattr(self,linekey+'_original',linedata_original)
            self.available_lines.append(linekey)
        data.close()


class Datafile_apo_demetra:
    observatory = 'APO_DEMETRA'
    # linelist = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
    linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
     ['Hb', 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'],
     ['He_I', 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'],
     ['He_I', 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5875'],
     ['He_II', 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4541'],
     ['He_II', 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4685'],
     ['He_II', 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5411'],
     ['He_I', 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'],
     ['He_I', 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4921'],
     ['He_I', 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'],
     ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'],
     ['C_IV', 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]

    def __init__(self, file,ll=None,v_rad = 18.5,i='n/a',mark = 0):
        if ll == None:
            pass
        else:
            linelist = 'invullen'

        fn = os.path.basename(file)
        data = pf.open(file)
        self.original_filepath = file
        self.i =i
        self.mark = mark
        self.mark_explanation = '0 = no weird stuff,      1 = not usable due to very poor SNR,    2 = Not usable for EW, TVS, Quotient due to insufficient SNR,   3 =   Shows weird feature in Halpha'
        self.filename = fn[:fn.rfind(".")]
        self.header = data[0].header
        self.time_and_date = airmass.timeanddate2(self.header['DATE-OBS'])
        self.baricentric_correction, self.HJD = airmass.barcor(file,JDOBS=self.header['JD-MID'])
        self.phase =  airmass.aphase(self.HJD)
        self.exptime = airmass.exposuretime(file)
        self.airmass, self.alt, JD = airmass.airmass(file,JDOBS=self.header['JD-MID'])
        try:
            frwl = airmass.fitfraun_demetra(file)
        except RuntimeError:
            frwl = 5895.92
        self.fwl = frwl
        self.velshift = 299792.458*(self.fwl-5895.92)/5895.92
        naxis1 = self.header['NAXIS1']
        crval1 = self.header['CRVAL1']
        cdelt1 = self.header['CDELT1']
        self.wl_original = np.arange(naxis1) * cdelt1 + crval1 - v_rad / 299792.458
        self.flux_original = data[0].data
        self.wl_rebin, self.flux_rebin = airmass.rebin2(self.wl_original,self.flux_original)
        self.available_lines = []
        self.snr_original =airmass.snr(self.wl_original,self.flux_original)
        self.snr = airmass.snr(self.wl_rebin,self.flux_rebin)
        for line in self.linelist:
            linedata,linekey = line_data(line,self.wl_rebin,self.flux_rebin,self.observatory,self.snr, self.baricentric_correction,-18.5)
            linedata_original, lk = line_data(line,self.wl_original,self.flux_original,self.observatory,self.snr_original,0,0)
            setattr(self,linekey , linedata)
            setattr(self,linekey+'_original',linedata_original)
            self.available_lines.append(linekey)
        data.close()



# class datafile_stack_apo:
#     observatory = 'APO'
#     # linelist = [['Ha', 35, 6562.819, 6554, 6556, 6578, 6579, r'H$\alpha$ 6563'], ['Hb', 35, 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, r'H$\beta$ 4861'], ['He_I', 35, 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, 'He I 4713'], ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'], ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'], ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'], ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'], ['He_I', 35,4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, 'He I 4471'] , ['He_I',35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'] , ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'] ,['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'], ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
#     linelist = [['Ha', 35, 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
#      ['Hb', 35, 4861.333, 4828.0, 4839.0, 4880.0, 4891.0, r'H$\beta$ 4861'],
#      ['He_I', 35, 4713.1457, 4701, 4703, 4718, 4720, 'He I 4713'],
#      ['He_I', 35, 5875.621, 5863.0, 5864.5, 5892.7, 5894.6, 'He I 5876'],
#      ['He_II', 35, 4541.6, 4523, 4529, 4546, 4548.5, 'He II 4542'],
#      ['He_II', 35, 4685.804, 4671.5, 4672.2, 4693.3, 4694.3, 'He II 4686'],
#      ['He_II', 35, 5411.521, 5405.2, 5406.6, 5425.0, 5428.2, 'He II 5412'],
#      ['He_I', 35, 4471.4802, 4459.0, 4462, 4475.5, 4478.5, 'He I 4471'],
#      ['He_I', 35, 4921.93, 4910, 4913, 4928.2, 4931.5, 'He I 4922'],
#      ['He_I', 35, 6678.15, 6656, 6660, 6690, 6695, 'He I 6678'],
#      ['O_III', 35, 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, 'O III 5592'],
#      ['C_IV', 35, 5801.33, 5793.8, 5796.2, 5817.1, 5819.5, 'C IV 5801']]
#
#     def __init__(self, stack, v_rad=18.5, i='n/a', mark=0):
#         middleIndex = int((len(stack) - 1) / 2)
#         self.header = stack[middleIndex].header
#         self.time_and_date = airmass.timeanddate2(self.header['DATE-OBS'])
#         HJDs = []
#         for file in stack:
#