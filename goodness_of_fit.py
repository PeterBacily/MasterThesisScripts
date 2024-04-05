from __future__ import division
import matplotlib.pyplot as plt
import glob
import pyfits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
import airmass
import scipy.stats as ss
from scipy.optimize import *
from SavitzkyGolay import savitzky_golay
from PyAstronomy import pyasl
print 'asdf'
nrow = 6
ncol = 2
# fig, axs = plt.subplots(nrows=nrow, ncols=ncol,sharex=True,figsize=(10,10))

def my_sin(x,  amplitude, phase, offset):
    return amplitude*np.sin(2*np.pi*x  + phase)  + offset

def my_line(x,a,b):
    return a*x+b

def flat_line(x,a):
    return a

rawfilelist = glob.glob('C:\peter\School\LaPalma\Data\All_data\zet Ori/*.fits')



datafolder = r'C:\peter\School\Master Scriptie\Data\EW\lp/'

datafile_list = glob.glob(datafolder + '*_sg.npy')
f = open(datafolder + 'line_Stats_snr_full_test3.txt', 'w')
f.write('line name'+' \t '+'chi2 sin' +'\t'+ 'AIC Sin' +'\t'+ 'AIC line' +'\t'+ 'AIC flat' +'\n')
###############################
def stats(dat):
    ew = dat[0]
    print ew
    phase = dat[1]
    error = dat[2]
    # error = np.std(ew)
    lineused = dat[3][0]
    # print lineused
    linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
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
    chisq2 = airmass.redchisqg(ew,data_fit_sin,deg=3)

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
    dif1 = p_sin/p_line
    dif2 = p_sin/p_flat
    line_to_write = linename+' \t '+ str(chisq2)+ '\t' +str(AIC_1) + '\t' +str(AIC_2) + '\t' +str(AIC_3) + '\n'
    # print line_to_write
    f.write(line_to_write)
    print '#######'
    print linename
    print chi2_sin-chi2_line
    print dif1
    print '#######'

for file in datafile_list:
    dat = np.load(file)
    stats(dat)

##########################

# for file in datafile_list:
#     file1 =np.load(file)
#     ew = file1[0]
#     print ew
#     phase = file1[1]
#     error = file1[2]
#     # error = np.std(ew)
#     lineused = file1[3][0]
#     # print lineused
#     linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
#     p1=[0.005,0.5, 1]
#     fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
#     p2 = [0,1]
#     fit2 = curve_fit(my_line, phase, ew, p0=p2,  sigma=error, absolute_sigma=True)
#     p3 =[1]
#     fit3 = curve_fit(flat_line , phase, ew, p0=p3,  sigma=error, absolute_sigma=True)
#     # t= np.linspace(0,1,50)
#     data_fit_sin = my_sin(np.array(phase), *fit1[0])
#     data_fit_line = my_line(np.array(phase), *fit2[0])
#     data_fit_flat = flat_line(np.array(phase), *fit3[0])
#
#     residuals = (ew - data_fit_sin)/error
#     chisq = np.sum(residuals**2) / (len(ew)-1-3)
#     chisq2 = airmass.redchisqg(ew,data_fit_sin,deg=3)
#
#     k_1 = 3
#     res_1 = (ew-data_fit_sin)
#     SSR_1 = np.sum(res_1**2)
#     N_1 = len(ew)
#     s2_1 = SSR_1 / N_1
#     L_1 = ( 1.0/np.sqrt(2*np.pi*s2_1) ) ** N_1 * np.exp( -SSR_1/(s2_1*2.0) )
#
#     k_2 = 2
#     res_2 = (ew-data_fit_line)
#     SSR_2 = np.sum(res_2**2)
#     N_2 = len(ew)
#     s2_2 = SSR_2 / N_2
#     L_2 = ( 1.0/np.sqrt(2*np.pi*s2_2) ) ** N_2 * np.exp( -SSR_2/(s2_2*2.0) )
#
#     k_3 = 1
#     res_3 = (ew-data_fit_flat)
#     SSR_3 = np.sum(res_3**2)
#     N_3 = len(ew)
#     s2_3 = SSR_3 / N_3
#     L_3 = ( 1.0/np.sqrt(2*np.pi*s2_3) ) ** N_3 * np.exp( -SSR_3/(s2_3*2.0) )
#
#     AIC_1 = 2*k_1 - 2*np.log(L_1)
#     AIC_2 = 2*k_2 - 2*np.log(L_2)
#     AIC_3 = 2*k_3 - 2*np.log(L_3)
#     probfactor = np.exp((AIC_1-AIC_2)/2)
#
#     # print 'ln(L) =', np.log( L )
#
#     # print chisq, chisq2, chisq-chisq2
#     chi2_sin,p_sin = ss.chisquare(ew,data_fit_sin,ddof=3)
#     chi2_line,p_line = ss.chisquare(ew,data_fit_line,ddof=2)
#     chi2_flat,p_flat = ss.chisquare(ew,data_fit_flat,ddof=1)
#     dif1 = p_sin/p_line
#     dif2 = p_sin/p_flat
#     line_to_write = linename+' \t '+ str(chisq2)+ '\t' +str(AIC_1) + '\t' +str(AIC_2) + '\t' +str(AIC_3) + '\n'
#     # print line_to_write
#     f.write(line_to_write)
#     print '#######'
#     print linename
#     print chi2_sin-chi2_line
#     print dif1
#     print '#######'



    # print ew
    # print data_fit_sin
    # print ew-data_fit_sin
    # print error
#################################


# # f = open(datafolder + 'llhr_manual.txt', 'w')
# # f.write('line name'+' \t '+'chi2 sin' +'\n' )
# for file in datafile_list:
#     file1 =np.load(file)
#     ew = file1[0]
#     # print ew
#     phase = file1[1]
#     error = file1[2]
#     lineused = file1[3][0]
#     # print lineused
#     linename = (lineused[0] +' '+ str(int(lineused[1]))).replace('_',' ')
#     p1=[0.005,0.5, 1]
#     fit1 = curve_fit(my_sin, phase, ew, p0=p1,  sigma=error, absolute_sigma=True)
#     p2 = [0,1]
#     fit2 = curve_fit(my_line, phase, ew, p0=p2,  sigma=error, absolute_sigma=True)
#     p3 =[1]
#     fit3 = curve_fit(flat_line , phase, ew, p0=p3,  sigma=error, absolute_sigma=True)
#     # t= np.linspace(0,1,50)
#     data_fit_sin = my_sin(np.array(phase), *fit1[0])
#     data_fit_line = my_line(np.array(phase), *fit2[0])
#     data_fit_flat = flat_line(np.array(phase), *fit3[0])
#
#     llh_sin,p_sin = ss.power_divergence(ew, f_exp=data_fit_sin, ddof=3, axis=None, lambda_="log-likelihood")
#     llh_line,p_line = ss.power_divergence(ew, f_exp=data_fit_line, ddof=2, axis=None, lambda_="log-likelihood")
#     llh_flat,p_flat = ss.power_divergence(ew, f_exp=data_fit_flat, ddof=1, axis=None, lambda_="log-likelihood")
#     print '--------------------------'
#     print linename
#     print llh_sin,llh_line,llh_flat
#     print p_sin,p_line,p_flat
#     llhr = 2*(llh_line-llh_sin)
#     print llhr
#     p_llhr = ss.chisqprob(llhr,1)
#     print p_llhr
#     print '--------------------------'


# def make_lp_TVS(filelist,line,v_rad=18.5):
    # linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, [6615.0, 6618.3]], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, [4877.0, 4883.0]], ['Hy', 4340.472, 4322, 4324, 4357, 4360, [4298.0, 4302.0]], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, [4031.0, 4038.0]], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, [4458.0, 4463.0]], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, [4734.0, 4740.0]], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, [5819.0, 5833.0]], ['He_II', 4541.6, 4498, 4499, 4580, 4581, [4579.2, 4587.6]], ['He_II', 4685.804, 4679, 4680, 4690, 4691, [4718.0, 4740.0]], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, [5370.0, 5380.0]], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, [5564.3, 5575.0]], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, [5831.2, 5833.1]]]
    # lws = []
    # lfs = []
    # std_exp_list = []
    # std_exp_weights = []
    # for file in filelist:
    #     datafile = pf.open(file)
    #     header = datafile[0].header
    #     naxis1 = header['NAXIS1']
    #     crval1 = header['CRVAL1']
    #     cdelt1 = header['CDELT1']
    #     flux = datafile[0].data
    #     wl = np.exp(np.arange(naxis1)*cdelt1 + crval1 - v_rad/299792.458)
    # wl, TVS = airmass.TVS_LaPalma(filelist,line,)
def make_large_TVS(filelist,save_folder):
    linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, [6615.0, 6618.3]], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, [4877.0, 4883.0]], ['Hy', 4340.472, 4322, 4324, 4357, 4360, [4298.0, 4302.0]], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, [4031.0, 4038.0]], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, [4458.0, 4463.0]], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, [4734.0, 4740.0]], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, [5819.0, 5833.0]], ['He_II', 4541.6, 4498, 4499, 4580, 4581, [4579.2, 4587.6]], ['He_II', 4685.804, 4679, 4680, 4690, 4691, [4718.0, 4740.0]], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, [5370.0, 5380.0]], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, [5564.3, 5575.0]], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, [5831.2, 5833.1]]]
    for line in linelist:
        wl, TVS = airmass.TVS_LaPalma(filelist,line,3800,9000)
        linename = line[0]+str(int(round(line[1])))
        dat = np.array([wl,TVS])
        np.save(save_folder+'TVS_'+linename+'.npy',dat)

# make_large_TVS(rawfilelist,r'C:\peter\School\Master Scriptie\Data\TVS\LaPalma')
#
# tvsfilelist = glob.glob(r'C:\peter\School\Master Scriptie\Data\TVS\*.npy')
# print tvsfilelist
# linelist = [['Ha', 6562.819, 6551, 6552, 6578, 6579, [6615.0, 6618.3]], ['Hb', 4861.333, 4838.0, 4839.0, 4880.0, 4881.0, [4877.0, 4883.0]], ['Hy', 4340.472, 4322, 4324, 4357, 4360, [4298.0, 4302.0]], ['He_I', 4026.1914, 4018.1914, 4022.1914, 4030.1914, 4034.1914, [4031.0, 4038.0]], ['He_I', 4471.4802, 4466.0, 4467.0, 4475.0, 4476.0, [4458.0, 4463.0]], ['He_I', 4713.1457, 4708.15, 4709.15, 4718.15, 4719.15, [4734.0, 4740.0]], ['He_I', 5875.621, 5863.0, 5864.5, 5885.0, 5885.8, [5819.0, 5833.0]], ['He_II', 4541.6, 4498, 4499, 4580, 4581, [4579.2, 4587.6]], ['He_II', 4685.804, 4679, 4680, 4690, 4691, [4718.0, 4740.0]], ['He_II', 5411.521, 5400.7, 5401.7, 5422.0, 5423.0, [5370.0, 5380.0]], ['O_III', 5592.37, 5586.0, 5587.0, 5598.0, 5599.0, [5564.3, 5575.0]], ['C_IV', 5801.33, 5794.6, 5795.6, 5807.1, 5808.1, [5831.2, 5833.1]]]
#
# for file in tvsfilelist:
#     data = np.load(file)
#     wl = data[0]
#     tvs = data[1]
#     linewl  = file.split('TVS_', 1)[1].split('.')[0]
#     central_wavelength = float(linewl[-4:])
#     line = sorted(linelist, key=lambda x: abs(x[1]-central_wavelength))[0]
#     a=line[2]
#     b=line[3]
#     c=line[4]
#     d=line[5]
#     normwave = np.hstack((wl[(wl>a)&(wl<b)],wl[(wl>c)&(wl<d)]))
#     normtvs = np.hstack((tvs[(wl>a)&(wl<b)],tvs[(wl>c)&(wl<d)]))
#     start = line[6][0]
#     stop = line[6][1]
#     wl_piece = wl[(wl>start)&(wl<stop)]
#     tvs_piece = tvs[(wl>start)&(wl<stop)]
#     tvs_avg = str(np.average(tvs_piece))
#     print np.average(normtvs), tvs_avg
#     plt.title(file.split('TVS_', 1)[1].split('.')[0])
#     plt.plot(normwave,normtvs,label=tvs_avg)
#     plt.legend()
#     plt.show()
#     plt.close()


