import glob
import os
import pathlib
import numpy as np
import tqdm
import airmass
import create_datafiles
import datareduc
import Datafile_class
import open_masterfiles
import pickle
import Path_check
import math
import random
folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
Path_check.dir_check(folder_of_this_file)

[converted_Data_folder, Data_folder, Plots_folder, Scripts_folder] = Path_check.dir_paths(folder_of_this_file)

# filelist = open_masterfiles.mercator(str(converted_Data_folder)+r'\dataset_omar\\')

# linelist = open_masterfiles.open_linelist(str(converted_Data_folder)+r'\linelists\final_lls\linelist_no_hy.txt')
def open_linelist(filepath):
    workfileresource = open(filepath, 'rb')
    list = pickle.load(workfileresource)
    workfileresource.close()
    return list

def make_folderpath(parent_path = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator_apolines\selection',
                              R=10000,snr_desired = 1000, selectionstring = ''):
    if selectionstring == None or selectionstring =='':
        selectionpath = r'\All'
    else:
        selectionpath = '\selection_'+selectionstring
    if R is None and snr_desired is None:
        degstring = ''
        degfolderstring = r'\original\\'
    else:
        degstring = 'R'+str(int(R))+'_snr'+str(int(snr_desired*5))
        degfolderstring = r'\degraded\\'
    path_with_deg = parent_path +selectionpath+ degfolderstring
    savefolder = path_with_deg+degstring
    pathlib.Path(savefolder).mkdir(parents=True, exist_ok=True)
    return(savefolder+r'\\')

def flatten_list(nested_list):
    flat_list = []
    for item in nested_list:
        if isinstance(item, list):
            flat_list.extend(flatten_list(item))
        else:
            flat_list.append(item)
    return flat_list

def run_mdg_deg(filelist,linelist,savefolder,R=10000,snr_desired = 1000,vmin = -800,vmax = 800,k=1,selectionstring='All'):
    for i,line in enumerate(linelist):
        linekey = 'line' + str(int(line[k]))

        data_grid = create_datafiles.make_data_grid_with_degradation(filelist,linekey, vmin,vmax,R=R,snr_desired=snr_desired,rebin_size=0.5,selectionstring=selectionstring)
        savename = savefolder+'data_grid_'+line[0]+'_'+str(int(line[1]))+str(vmin)+'_'+str(vmax)+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(data_grid, workfileresource)
        workfileresource.close()

def run_mdg(filelist,linelist,datagridfolder,vmin = -800,vmax = 800,k=1,selectionstring='All'):
    savefolder = datagridfolder
    rb=0.5
    for i,line in enumerate(linelist):
        linekey = 'line' + str(int(line[k]))
        data_grid = create_datafiles.make_data_grid(filelist,linekey, vmin,vmax,rebin_size=rb,selectionstring=selectionstring)
        savename = savefolder+'data_grid_'+line[0]+'_'+str(int(line[1]))+'rebin'+str(rb)+'_vlim'+str(vmin)+'_'+str(vmax)+'.txt'
        workfileresource = open(savename, 'wb')
        pickle.dump(data_grid, workfileresource)
        workfileresource.close()


def run_mlb(data_brick_folderpath,ls_grid_folderpath,frequencyarray = None):
    savefolder = ls_grid_folderpath
    filelist = glob.glob(data_brick_folderpath+'\*.txt')
    for filepath in filelist:
        lombscl_dict,output_filename = create_datafiles.make_ls_brick(filepath,frequencyarray=frequencyarray)
        output_filepath = savefolder+r'\\'+output_filename
        workfileresource = open(output_filepath, 'wb')
        pickle.dump(lombscl_dict, workfileresource)
        workfileresource.close()


def plot_LS_grid(ls_brick_folder,LS_plot_folder,v_min_ls=-500,v_max_ls=500):
    filelist = glob.glob(ls_brick_folder + '\*.txt')
    for file in filelist:
        datareduc.ls_brick_plotter(file, v_min_ls, v_max_ls, plotsavefolder=LS_plot_folder, show='off', save='on')


def make_sumplot(LS_brick_folder,savefolder,excludehy=True):
    for bool in [True,False]:
        datareduc.ls_sum_plotter(
            LS_brick_folder, -500,
            500, plotsavefolder=savefolder, show='off', save='on', SG=bool, SGwindowsize=201,excludehy=excludehy)
        datareduc.ls_sum_plotter(
            LS_brick_folder, -500,
            500, plotsavefolder=savefolder, show='off', save='on', SG=bool, SGwindowsize=201,excludehy=excludehy)


def run_full_pipeline(filelist,linelist,datagrid_folder,LS_brick_folder,LS_plot_folder,Sumplot_folder,R=None,SNR_desired=None,vmin=-800,vmax=800,v_min_ls=-500,v_max_ls=500,frequency_array=None,selectionstring='All'):
    if R is None and SNR_desired is None:
        run_mdg(filelist,linelist,datagrid_folder,vmin=vmin,vmax=vmax,selectionstring=selectionstring)
    else:
        run_mdg_deg(filelist,linelist,datagrid_folder,R,SNR_desired,vmin=vmin,vmax=vmax,selectionstring=selectionstring)
    run_mlb(datagrid_folder,LS_brick_folder,frequencyarray=frequency_array)
    plot_LS_grid(LS_brick_folder,LS_plot_folder,v_min_ls=v_min_ls,v_max_ls=v_max_ls)
    make_sumplot(LS_brick_folder,Sumplot_folder)
def feed_selection_into_pipeline(filelist, LS_brick_folder,LS_plot_folder,datagrid_folder,Sumplot_folder,plot_selectionstring,R=None,SNR_desired=None):
    linelist = open_masterfiles.open_linelist(str(converted_Data_folder) + r'\linelists\final_lls\linelist_no_hy.txt')
    # filelistpaths = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\test\pipelinetest\masterfiles\\'+'*.txt')
    # filelist_nested = airmass.group_observations(open_masterfiles.mercator(manual_filelist=filelistpaths))
    # Wat hierboven staat als parameters meegeven
    filelist_flattened = flatten_list(filelist)
    filelist_sorted = sorted(filelist_flattened, key=lambda obj: obj.HJD)
    run_full_pipeline(filelist_sorted,linelist,datagrid_folder,LS_brick_folder,LS_plot_folder,Sumplot_folder,R=R,SNR_desired=SNR_desired,selectionstring=plot_selectionstring)

def generate_selection_strings(filelist,selection_method='Group',randomremoval_percent_removed = '0'):
    if selection_method=='Group':
        filelist_flattened = flatten_list(filelist)
        filelist_sorted = sorted(filelist_flattened, key=lambda obj: obj.HJD)
        firstfile = filelist_sorted[0]
        lastfile = filelist_sorted[-1]
        [yr_f,m_f,d_f,t_f]=airmass.split_date(firstfile.header['DATE-OBS'])[0]
        [yr_l, m_l, d_l, t_l] = airmass.split_date(lastfile.header['DATE-OBS'])[0]
        selectionstring = 'from '+yr_f+'-'+m_f+'-'+d_f+' to '+yr_l+'-' +m_l+'-'+ d_l
        foldernamestring = yr_f+m_f+'-'+yr_l+ m_l
    elif selection_method=='Random':
        selectionstring = 'Randomly removed '+randomremoval_percent_removed+'%'
        foldernamestring ='RandomRemoved'+randomremoval_percent_removed
    else:
        print('Selection method needs to be "Group" or "Random"')
    return selectionstring,foldernamestring

def run_selection(selection,R=10000,snr_desired=20,sm='Group',rm=''):
    plot_selectionstring,folder_selectionstring = generate_selection_strings(selection,selection_method=sm,randomremoval_percent_removed=rm)
    databrick_base_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\data_bricks\mercator_apolines_rebin05'
    lsbrick_base_folder = r'D:\peter\Master_Thesis\Datareduction\Converted_Data\ls_bricks\mercator_apolines_rebin05'
    ls_plot_base_folder = r'D:\peter\Master_Thesis\Datareduction\Plots\LS_periodogram\normal\mercator_apolines_rebin05'
    sum_ls_plot_base_folder = r'D:\peter\Master_Thesis\Datareduction\Plots\LS_periodogram\summed\mercator_apolines_rebin05'
    databrickfolder = make_folderpath(databrick_base_folder,R=R,snr_desired=snr_desired,selectionstring=folder_selectionstring)
    lsbrickfolder = make_folderpath(lsbrick_base_folder,R=R,snr_desired=snr_desired,selectionstring=folder_selectionstring)
    lsplotfolder = make_folderpath(ls_plot_base_folder,R=R,snr_desired=snr_desired,selectionstring=folder_selectionstring)
    sumlsplotfolder = make_folderpath(sum_ls_plot_base_folder,R=R,snr_desired=snr_desired,selectionstring=folder_selectionstring)
    feed_selection_into_pipeline(selection,LS_brick_folder=lsbrickfolder,datagrid_folder=databrickfolder,LS_plot_folder=lsplotfolder,Sumplot_folder=sumlsplotfolder,plot_selectionstring=plot_selectionstring, R=R,SNR_desired=snr_desired)
# run_test_selection(open_masterfiles.mercator(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\'))





def run_pipeline_variations(observations, pipeline):
    param_sets = [(None, None)] + [(10000, snr) for snr in [6, 8, 10, 15, 20, 30]]

    # --- Group-based removal ---
    grouped = airmass.group_observations(observations)
    total_groups = len(grouped)

    for i in tqdm.tqdm(range(total_groups)):
        subset = [obs for group in grouped[i:] for obs in group]

        for R, snr in param_sets:
            # print(f"Group removal: {len(subset)} observations | R={R}, snr_desired={snr}")
            pipeline(subset, R=R, snr_desired=snr,sm='Group',rm='')

    # --- Percentage-based random removal ---
    total = len(observations)
    step_size = math.ceil(total * 0.10)
    remaining = observations[:]
    removed = set()

    num_steps = total // step_size

    for step in tqdm.tqdm(range(1, num_steps + 1)):  # skip step 0 (full dataset)
        # Determine how many to remove this step
        to_remove = random.sample([obs for obs in observations if obs not in removed], step_size)
        removed.update(to_remove)

        subset = [obs for obs in observations if obs not in removed]

        for R, snr in param_sets:
            # print(f"Random removal: {len(subset)} observations | R={R}, snr_desired={snr}")
            pipeline(subset,R=R, snr_desired=snr,sm='Random',rm=f'{step*10:.0f}')

def run_pipeline_variations_test(observations, pipeline):
    param_sets = [(10000, 20),(None,None)]

    # ------------------------
    # ✅ GROUP-BASED TESTS (simulate only 2 groups)
    # ------------------------

    grouped_all = airmass.group_observations(observations)
    test_groups = grouped_all[-1:]  # Keep only the last 2 groups

    for i in range(len(test_groups), 0, -1):  # 2 groups, then 1
        subset = [obs for group in test_groups[-i:] for obs in group]

        for R, snr in param_sets:
            print(f"Group-based test: {i} groups | {len(subset)} observations | R={R}, snr_desired={snr}")
            pipeline(subset,R=R, snr_desired=snr,sm='Group',rm='')

    # ------------------------
    # ✅ RANDOM REMOVAL TESTS (simulate 80% and 90% removal)
    # ------------------------

    total = len(observations)
    test_fractions = [ 0.10]  # Keep 20% and 10%

    for frac in test_fractions:
        subset_size = max(1, int(total * frac))
        subset = random.sample(observations, subset_size)

        for R, snr in param_sets:
            print(f"Random removal test: keeping {frac*100:.0f}% | {len(subset)} observations | R={R}, snr_desired={snr}")
            pipeline(subset, R=R, snr_desired=snr,sm='Random',rm=f'{frac*100:.0f}')

def run_pipeline_firstgroup(observations, pipeline):
    grouped_all = airmass.group_observations(observations)
    subset = grouped_all[0]
    param_sets = [(None, None)] + [(10000, snr) for snr in [6, 8, 10, 15, 20, 30]]
    for R, snr in tqdm.tqdm(param_sets):
        print(f"Group-based test: {len(subset)} observations | R={R}, snr_desired={snr}")
        pipeline(subset, R=R, snr_desired=snr, sm='Group', rm='')

# run_pipeline_variations(open_masterfiles.mercator(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\'),run_selection)
run_pipeline_firstgroup(open_masterfiles.mercator(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\dataset_omar\\'),run_selection)