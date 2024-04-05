import pathlib
import os
folder_of_this_file = os.path.dirname(os.path.abspath(__file__))


def up(path_string):
    path = pathlib.Path(path_string)
    return str(path.parent.absolute())


def dir_check(path_string=folder_of_this_file):
    dirlist = os.listdir(up(path_string))
    check_folder = set(['Converted_Data', 'Data', 'Plots', 'Scripts']).issubset(set(dirlist))
    if not check_folder:
        raise Exception("unexpected directory structure.\n current directory"+up(path_string)+"\n current sub folders: "
                        + str(dirlist) + "\n expected sub folders: ['Converted_Data', 'Data', 'Plots', 'Scripts']")

def dir_paths(path_string=folder_of_this_file):
    dir_check(path_string)
    a = list(pathlib.Path(up(path_string)).glob('*'))
    return a


dir_paths()


# copy code below into file to execute path check of that file
# import os
# import Path_check
# folder_of_this_file = os.path.dirname(os.path.abspath(__file__))
# Path_check.dir_check(folder_of_this_file)
