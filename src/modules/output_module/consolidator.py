import os


def copy_files(list_of_files, destination):
    for file in list_of_files:
        os.system(f'cp {file} {destination}')

def move_files(list_of_files, destination):
    for file in list_of_files:
        os.system(f'mv {file} {destination}')

def rename_file(file_name_in, file_name_out):
    '''
    INPUT:
        File name including full path
    OUTPUT:
        New file to create including full path
    '''
    ### TODO: file not found
    os.system(f'mv {file_name_in} {file_name_out}')

