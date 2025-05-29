from __future__ import division
import numpy as np
# os.environ['PYSYN_CDBS'] = 'C:\Users\Peter\Anaconda2\envs\p27\Lib\site-packages\pysynphot\data\cbds\grp\hst\cdbs'
import scipy
# beginning = 1
# end=12
# bigstep = 1
# smallstep = 1/100
# ss_start = 8
# ss_end=9
# a = np.arange(beginning,ss_start,bigstep)
# b=np.arange(ss_start,ss_end,smallstep)
# c= np.append(a,b)
# d =np.arange(ss_end,end+bigstep,bigstep)
# e=np.append(c,d)
# print(e)
#
# apo_lines = ['line6562', 'line4713', 'line5411', 'line5801', 'line4541', 'line4685', 'line5875', 'line5592',
#                  'line4861', 'line4921', 'line6678', 'line4471']
# mercator_lines = ['line5875', 'line4861', 'line4340', 'line6562', 'line6678', 'line5592', 'line4026', 'line4471',
#                       'line4685', 'line4921', 'line5801', 'line4713', 'line5411', 'line4541']

a = np.array([1,2,3,5,7,9,20,30,40])
b= np.diff(a)
print(b)