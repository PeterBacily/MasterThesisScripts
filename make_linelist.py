import pickle
import airmass
import open_masterfiles
import Datafile_class
import os
import Path_check


folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)
[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)
def open_linelist(filepath):
    workfileresource = open(filepath, 'rb')
    list = pickle.load(workfileresource)
    workfileresource.close()
    return list

def make_linelist(list,filepath):
    workfileresource = open(filepath, 'wb')
    pickle.dump(list, workfileresource)
    workfileresource.close()


linelist_apo = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
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


linelist_apo2 = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
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

linelist_apo_new = [['Ha', 6562.819, 6550.2, 6551.2, 6578, 6579, r'H$\alpha$ 6563'],
     ['Hb', 4861.333, 4848.0, 4849.0, 4877.0, 4878.0, r'H$\beta$ 4861'],
     ['He_I', 4713.1457, 4708.3, 4709.3, 4718, 4719, 'He I 4713'],
     ['He_I', 5875.621, 5867, 5868, 5881.7, 5882.7, 'He I 5875'],
     ['He_II', 4541.6, 4537.5, 4538.5, 4546.6, 4547.6, 'He II 4541'],
     ['He_II', 4685.804, 4682.3, 4683.3, 4690.5, 4691.5, 'He II 4685'],
     ['He_II', 5411.521, 5405.2, 5406.2, 5417.0, 5418.2, 'He II 5411'],
     ['He_I', 4471.4802, 4463.0, 4464, 4476.2, 4477.2, 'He I 4471'],
     ['He_I', 4921.93, 4914.2, 4915.2, 4928.2, 4929.2, 'He I 4921'],
     ['He_I', 6678.15, 6670, 6671, 6685, 6686, 'He I 6678'],
     ['O_III', 5592.37, 5588.3, 5589.3, 5597.5, 5598.5, 'O III 5592'],
     ['C_IV', 5801.33, 5793.8, 5794.8, 5817.1, 5818.1, 'C IV 5801']]

linelist_standard = [['Ha', 6562.819, 6551, 6552, 6578, 6579, r'H$\alpha$ 6563'],
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

linelist_apo_vshift = [['Ha', 6562.819, 6549.2, 6550.2, 6577.0, 6578.0, r'H$\alpha$ 6563'],
['Hb', 4861.333, 4847.3, 4848.3, 4876.3, 4877.3, r'H$\beta$ 4861'],
['He_I', 4713.1457, 4707.6, 4708.6, 4717.3, 4718.3, 'He I 4713'],
['He_I', 5875.621, 5866.1, 5867.1, 5880.8, 5881.8, 'He I 5875'],
['He_II', 4541.6, 4536.8, 4537.8, 4545.9, 4546.9, 'He II 4541'],
['He_II', 4685.804, 4681.6, 4682.6, 4689.8, 4690.8, 'He II 4685'],
['He_II', 5411.521, 5404.4, 5405.4, 5416.2, 5417.4, 'He II 5411'],
['He_I', 4471.4802, 4462.3, 4463.3, 4475.5, 4476.5, 'He I 4471'],
['He_I', 4921.93, 4913.5, 4914.5, 4927.5, 4928.5, 'He I 4921'],
['He_I', 6678.15, 6669.0, 6670.0, 6684.0, 6685.0, 'He I 6678'],
['O_III', 5592.37, 5587.5, 5588.5, 5596.7, 5597.7, 'O III 5592'],
['C_IV', 5801.33, 5792.9, 5793.9, 5816.2, 5817.2, 'C IV 5801']]

linelist_merc_vshift_revised = [['Ha', 6562.819, 6549.7, 6550.7, 6577.0, 6578.0, r'H$\alpha$ 6563'],
['Hb', 4861.333, 4847.3, 4848.3, 4876.3, 4877.3, r'H$\beta$ 4861'],
['He_I', 4713.1457, 4707.6, 4708.6, 4717.3, 4718.3, 'He I 4713'],
['He_I', 5875.621, 5866.3, 5867.3, 5884.4, 5885.3, 'He I 5875'],
['He_II', 4541.6, 4536.8, 4537.8, 4545.9, 4546.9, 'He II 4541'],
['He_II', 4685.804, 4681.6, 4682.6, 4689.8, 4690.8, 'He II 4685'],
['He_II', 5411.521, 5404.4, 5405.4, 5416.2, 5417.4, 'He II 5411'],
['He_I', 4471.4802, 4462.3, 4463.3, 4475.5, 4476.5, 'He I 4471'],
['He_I', 4921.93, 4913.5, 4914.5, 4927.5, 4928.5, 'He I 4921'],
['He_I', 6678.15, 6669.0, 6670.0, 6684.0, 6685.0, 'He I 6678'],
['O_III', 5592.37, 5587.5, 5588.5, 5596.7, 5597.7, 'O III 5592'],
['C_IV', 5801.33, 5792.9, 5793.9, 5816.2, 5817.2, 'C IV 5801']]

linelist_apo_vshift_revised = [['Ha', 6562.819, 6549, 6550.7, 6576.0, 6578.0, r'H$\alpha$ 6563'],
['Hb', 4861.333, 4847.3, 4848.3, 4876.3, 4877.3, r'H$\beta$ 4861'],
['He_I', 4713.1457, 4707.6, 4708.6, 4717.3, 4718.3, 'He I 4713'],
['He_I', 5875.621, 5866.3, 5867.3, 5883.5, 5884.5, 'He I 5875'],
['He_II', 4541.6, 4536.8, 4537.8, 4545.9, 4546.9, 'He II 4541'],
['He_II', 4685.804, 4681.6, 4682.6, 4689.8, 4690.8, 'He II 4685'],
['He_II', 5411.521, 5404.4, 5405.4, 5416.2, 5417.4, 'He II 5411'],
['He_I', 4471.4802, 4462.3, 4463.3, 4475.5, 4476.5, 'He I 4471'],
['He_I', 4921.93, 4913.5, 4914.5, 4927.5, 4928.5, 'He I 4921'],
['He_I', 6678.15, 6669.0, 6670.0, 6684.0, 6685.0, 'He I 6678'],
['O_III', 5592.37, 5587.5, 5588.5, 5596.7, 5597.7, 'O III 5592'],
['C_IV', 5801.33, 5792.9, 5793.9, 5816.2, 5817.2, 'C IV 5801']]

linelist_apo_vshift_revised = [['Ha', 6562.819, 6549, 6550.7, 6576.0, 6578.0, r'H$\alpha$ 6563'],
['Hb', 4861.333, 4847.3, 4848.3, 4876.3, 4877.3, r'H$\beta$ 4861'],
['He_I', 4713.1457, 4707.6, 4708.6, 4717.3, 4718.3, 'He I 4713'],
['He_I', 5875.621, 5866.3, 5867.3, 5883.5, 5884.5, 'He I 5875'],
['He_II', 4541.6, 4536.8, 4537.8, 4545.9, 4546.9, 'He II 4541'],
['He_II', 4685.804, 4681.6, 4682.6, 4689.8, 4690.8, 'He II 4685'],
['He_II', 5411.521, 5404.4, 5405.4, 5416.2, 5417.4, 'He II 5411'],
['He_I', 4471.4802, 4462.3, 4463.3, 4475.5, 4476.5, 'He I 4471'],
['He_I', 4921.93, 4913.5, 4914.5, 4927.5, 4928.5, 'He I 4921'],
['He_I', 6678.15, 6669.0, 6670.0, 6684.0, 6685.0, 'He I 6678'],
['O_III', 5592.37, 5587.5, 5588.5, 5596.7, 5597.7, 'O III 5592'],
['C_IV', 5801.33, 5792.9, 5793.9, 5816.2, 5817.2, 'C IV 5801'],
['Hy',4340.472,4326.0, 4330.3, 4355.0, 4362.2,r'H$\gamma$ 4861']]

boundaries = {'Ha6562':[0,0], 'Hb4861''Ha6562':[0,0], 'He_I4713''Ha6562':[0,0], 'He_I5875''Ha6562':[0,0],
              'He_II4541''Ha6562':[0,0], 'He_II4685''Ha6562':[0,0], 'He_II5411''Ha6562':[0,0], 'He_I4471''Ha6562':[0,0],
              'He_I4921''Ha6562':[0,0], 'He_I6678''Ha6562':[0,0], 'O_III5592''Ha6562':[0,0], 'C_IV5801''Ha6562':[0,0],
              'Hy4340''Ha6562':[0,0]}

a= open_linelist(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\linelists\linelist_apo_v_cor_3.txt')
print(a)
b= open_linelist("D:\peter\Master_Thesis\Datareduction\Converted_Data\linelists\linelist_merc_incl_Hy.txt")
print(b)


def deep_equal(x, y):
    if isinstance(x, list) and isinstance(y, list):
        if len(x) != len(y):
            return False
        return all(deep_equal(i, j) for i, j in zip(x, y))
    return x == y


def find_inserted_elements(a, b):
    ai = 0
    inserted = []

    for bi in range(len(b)):
        if ai < len(a) and deep_equal(a[ai], b[bi]):
            ai += 1  # Match found, move to next element in 'a'
        else:
            inserted.append(b[bi])  # Extra element in 'b'

    if ai == len(a):
        return inserted  # All elements of 'a' matched in order
    else:
        return None  # 'a' is not a subsequence of 'b'
print(find_inserted_elements(a,b))

quit()
linelist_apo_ha_test = [['Ha', 6562.819, 6549, 6550.7, 6576.0, 6578.0, r'H$\alpha$ 6563']]
# print(airmass.velocity_to_wl([-1000, -700, 1000, 1500],4340.472))
# datafiles = open_masterfiles.apo_demetra_orders(path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\demetra\with_orders\v_cor\snr_100\\',manual_filelist=None,sorted='off')
# usefile=datafiles[0]
# wl_factor=usefile.wl_offset_factor

linelist_file_name = 'linelist_merc_incl_Hy.txt'
linelist_file_path = str(converted_Data_folder)+ r'\linelists\\' + linelist_file_name
# test_file_name = str(converted_Data_folder)+r'\linelists\linelist_apo_v_cor_3.txt'
# iter = [2,3,4,5]
# for line in linelist_apo_new:
# 	for i in iter:
# 		a=line[i]
# 		b=round(a*wl_factor,1)
# 		line[i]=b
# 	print(line,', ')

# make_linelist(linelist_apo_vshift_revised,linelist_file_path)
# make_linelist(linelist_apo_vshift_revised,test_file_name)
# mylist = open_linelist(test_file_name)
# print(mylist)