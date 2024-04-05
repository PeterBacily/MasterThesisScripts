from __future__ import division
import matplotlib.pyplot as plt
import glob
import pyfits as pf
from astropy.time import Time
import math
import datetime
filelist = glob.glob('C:/peter/School/Master Scriptie/Data/RAW/*')
lapalmafilelist = glob.glob(r'C:/peter/School/Master Scriptie/Data/LaPalmaData/*')
period = 6.829
print filelist
file = pf.open(filelist[0])
for item in file[0].header:
    print item
jdstart = Time(2452734.2, format='jd').jd
jds = []
for item in filelist:
    file = pf.open(item)
    jds.append(((Time(file[0].header['DATE-OBS'], format='isot', scale='utc').jd - jdstart)% period)/period)

ljds = []
for item in lapalmafilelist:
    file = pf.open(item)
    ljds.append(((Time(file[0].header['DATE-OBS'], format='isot', scale='utc').jd - jdstart)% period)/period)
monday = ((Time('2016-04-28T21:00:00.123', format='isot', scale='utc').jd - jdstart)% period)/period
wd = ((Time('2016-03-17T21:00:00.123', format='isot', scale='utc').jd - jdstart)% period)/period
y=[1]*len(jds)
y2 = [1]*len(ljds)
plt.scatter(jds,y, c='g',edgecolor='g', label='Already Obtained')
plt.scatter(monday,1,c='b', label='today 21:00')
plt.scatter(ljds,y2, c='r',edgecolor='r', label='Data LaPalma')


plt.title('Phase Coverage for Zeta Ori')
plt.xlabel('Phase')
plt.yticks([])
plt.xlim(0,1)
plt.ylim(0,2)
plt.legend()
plt.show()


