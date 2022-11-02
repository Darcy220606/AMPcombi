#!/bin/python3

import os
import sys
import pathlib
from amp_database import download_DRAMP

def check_samplelist(samplelist, tools, path):
    if(samplelist==[]):
        print('<--sample-list> was not given, sample names will be inferred from directory names')
        for dirpath, subdirs, files in os.walk(path):
            for dir in subdirs:
                if (dir not in tools):
                    samplelist.append(dir)
        return list(set(samplelist))
    else:
        return samplelist

def check_pathlist(filepaths, samplelist, fileending, path):
    if(filepaths==[]):
        print('<--path-list> was not given, paths to AMP-results-files will be inferred')
        for sample in samplelist:
            pathlist = []
            for dirpath, subdirs, files in os.walk(path):
                for file in files:
                    if ((sample in dirpath) and ((list(filter(file.endswith, fileending))!=[]))):
                        pathlist.append(dirpath+'/'+file)
            filepaths.append(pathlist)
        return filepaths
    else:
        return filepaths

def check_faa_path(faa_path, samplename):
    path_list = list(pathlib.Path(faa_path).rglob(f"*{samplename}*.faa"))
    if (len(path_list)>1):
        sys.exit(f'AMPcombi interrupted: There is more than one .faa file for {samplename} in the folder given with --faa_path')
    elif(not path_list):
        sys.exit(f'AMPcombi interrupted: There is no .faa file containing {samplename} in the folder given with --faa_path')
    return path_list[0]

def check_ref_database(database):
    if((database==None) and (not os.path.exists('amp_ref_database'))):
        print('<--AMP_database> was not given, the current DRAMP general-AMP database will be downloaded and used')
        database = 'amp_ref_database'
        os.makedirs(database, exist_ok=True)
        db = database
        download_DRAMP(db)
        return db    
    elif ((not database==None)):
        if (os.path.exists(database)):
            db = database
            print(f'<--AMP_database> = ${db} is found and will be used')
            return db
        if (not os.path.exists(database)):
            sys.exit(f'Reference amp database path {database} does not exist, please check the path.')
    elif((database==None) and (os.path.exists('amp_ref_database'))):
        print('<--AMP_database> = DRAMP is already downloaded and will be reused')
        database = 'amp_ref_database'
        db = database
        return db

def check_path(path):
    return os.path.exists(path) #returns True or False

def check_directory_tree(path, tools, samplelist):
    print(f'Checking directory tree {path} for sub-directories \n ')
    # get first level of sub-directories, check if at least one is named by a tool-name
    subdirs_1 = [x for x in os.listdir(path) if x in tools]
    if (not subdirs_1):
        sys.exit(f'AMPcombi interrupted: First level sub-directories in {path} are not named by tool-names. Please check the directories names and the keys given in "--tooldict". \n ')
    else:
        print('First level sub-directories passed check.')
    # get second level of sub-directories, check if at least one is named by a sample-name
    subdirs_2 = []
    for dir in subdirs_1:
        subdirs = [x for x in os.listdir(path+dir) if x in samplelist]
        if (subdirs):
            subdirs_2.append(subdirs)
    if (not subdirs_2):
        sys.exit(f'AMPcombi interrupted: Second level sub-directories in {path} are not named by sample-names. Please check the directories names and the names given as "--sample_list" \n ')
    else:
        print('Second level sub-directories passed check')
    print('Finished directory check')

def check_input_complete(path, samplelist, filepaths, tools):
    # 1. Head folder does not exist and filepaths-list was not given
    if((not check_path(path)) and (not filepaths)):
        sys.exit('AMPcombi interrupted: Please provide the correct path to either the folder containing all amp files to be summarized (--amp_results) or the list of paths to the files (--path_list)')
    # 2. Head folder does not exist, filepaths-list was given but no samplelist
    elif((not check_path(path)) and (filepaths) and (not samplelist)):
        sys.exit('AMPcombi interrupted: Please provide a list of sample-names (--sample_list) in addition to --path_list')
    # 3. Head folder does not exist, filepaths- and samplelist are given:
    elif((not check_path(path)) and (not filepaths) and (not samplelist)):
        for file in filepaths:
            print(f'in check_input_complete the file in filepath is:')
            # 3.1. check if paths in filepath-list exist
            if(not check_path(file)):
                sys.exit(f'AMPcombi interrupted: The path {file} does not exist. Please check the --path_list input.')
            # 3.2. check if paths contain sample-names from samplelist
            if(not any(n in file for n in samplelist)):
                sys.exit(f'AMPcombi interrupted: The path {file} does not contain any of the sample-names given in --sample_list')
    # 4. Head folder and sample-list are given
    elif((check_path(path)) and (not samplelist)):
        check_directory_tree(path, tools, samplelist)
