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
        selectionpath = r'\no_selection'
    else:
        selectionpath = '\selection'
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
        print(i+1)
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
def feed_selection_into_pipeline(filelist, LS_brick_folder,LS_plot_folder,datagrid_folder,Sumplot_folder,selection_method='Group',randomremoval_percent_removed = '0'):
    linelist = open_masterfiles.open_linelist(str(converted_Data_folder) + r'\linelists\final_lls\linelist_no_hy.txt')
    # filelistpaths = glob.glob(r'D:\peter\Master_Thesis\Datareduction\Converted_Data\test\pipelinetest\masterfiles\\'+'*.txt')
    # filelist_nested = airmass.group_observations(open_masterfiles.mercator(manual_filelist=filelistpaths))
    # Wat hierboven staat als parameters meegeven
    filelist_flattened = flatten_list(filelist)
    filelist_sorted = sorted(filelist_flattened, key=lambda obj: obj.HJD)
    if selection_method=='Group':
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

    run_full_pipeline(filelist_sorted,linelist,datagrid_folder,LS_brick_folder,LS_plot_folder,Sumplot_folder,R=10000,SNR_desired=150,selectionstring=selectionstring)
