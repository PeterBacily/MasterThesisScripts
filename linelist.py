from __future__ import division
import matplotlib.pyplot as plt
import glob
import pyfits as pf
from astropy.time import Time
import math
import calendar
import numpy as np
from PyAstronomy import pyasl
from scipy import interpolate

ll = np.loadtxt(r'C:\peter\School\Master Scriptie\linelist.txt')
# for i,l in enumerate(ll[142:207]):
# for i,l in enumerate(ll):
#     print i, l[0]
a = sorted(ll[142:207], key=lambda row: row[2])
for line in a:
    print line[0], line[2], line[1]
#
# for l in a:
#     print l[2]