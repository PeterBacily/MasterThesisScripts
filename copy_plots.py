import os
import shutil


def flatten_folder_structure(basepath, dest_folder):
    """
    Copies all files from basepath (recursively) into dest_folder
    with names based on their relative paths, replacing folder
    separators with underscores.
    """
    if not os.path.exists(dest_folder):
        os.makedirs(dest_folder)

    for root, _, files in os.walk(basepath):
        for file in files:
            # Get relative path after basepath
            rel_path = os.path.relpath(root, basepath)

            if rel_path == ".":
                new_filename = file
            else:
                # Replace folder separators with underscores
                new_filename = rel_path.replace(os.sep, "_") + "_" + file

            src_file = os.path.join(root, file)
            dest_file = os.path.join(dest_folder, new_filename)

            # Handle possible filename collisions
            counter = 1
            base_name, ext = os.path.splitext(new_filename)
            while os.path.exists(dest_file):
                dest_file = os.path.join(dest_folder, f"{base_name}_{counter}{ext}")
                counter += 1

            shutil.copy2(src_file, dest_file)

    print(f"All files have been copied to {dest_folder}.")




def move_non_smoothed_files(source_folder, dest_folder):
    """
    Moves every file in source_folder whose filename does NOT contain 'smoothed'
    into dest_folder. Only affects the top-level files in the folder, not subfolders.
    """
    # Make sure the destination exists
    os.makedirs(dest_folder, exist_ok=True)

    for filename in os.listdir(source_folder):
        file_path = os.path.join(source_folder, filename)

        # Only act on files, not subfolders
        if os.path.isfile(file_path):
            if "smoothed" not in filename.lower():  # case-insensitive check
                new_path = os.path.join(dest_folder, filename)

                # Handle filename collisions
                counter = 1
                base, ext = os.path.splitext(filename)
                while os.path.exists(new_path):
                    new_path = os.path.join(dest_folder, f"{base}_{counter}{ext}")
                    counter += 1

                shutil.move(file_path, new_path)
                print(f"Moved: {file_path} → {new_path}")

    print("Move complete.")

# Example usage:
source_folder = r"D:\peter\Master_Thesis\Datareduction\Plots\LS_periodogram\All_specremoved_summed_smoothed"
dest_folder = r"D:\peter\Master_Thesis\Datareduction\Plots\LS_periodogram\All_specremoved_summed"
move_non_smoothed_files(source_folder, dest_folder)

# Example usage:
# basepath = r"D:\peter\Master_Thesis\Datareduction\Plots\LS_periodogram\summed\mercator_apolines_rebin05"
# dest_folder = r"D:\peter\Master_Thesis\Datareduction\Plots\LS_periodogram\All_specremoved_summed_smoothed"
# flatten_folder_structure(basepath, dest_folder)
