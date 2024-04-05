from __future__ import division
import pyfits as pf
import glob
import airmass
# import datareduc
import matplotlib.pyplot as plt
import numpy as np

mf = np.load(r'C:\peter\School\Master Scriptie\Data\masterfiles\\eShel.npy')
legend = mf[0]
data = mf[1]

print legend[2][1][1]
print data[8][2][1][1]
