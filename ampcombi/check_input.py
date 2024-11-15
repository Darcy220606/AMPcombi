#!/bin/python3

# This script verifies the inputs for ampcombi, optional and required.

import os
import sys
import pathlib
from amp_database import download_ref_db

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
    if(os.path.isdir(faa_path)):
        path_list = list(pathlib.Path(faa_path).rglob(f"*{samplename}*.faa"))
        if (len(path_list)>1):
            sys.exit(f'AMPcombi interrupted: There is more than one .faa file for {samplename} in the folder given with --faa')
        elif(not path_list):
            sys.exit(f'AMPcombi interrupted: There is no .faa file containing {samplename} in the folder given with --faa')
        return path_list[0]
    elif(os.path.isfile(faa_path)):
        return faa_path
    else:
        sys.exit(f'AMPcombi interrupted: The input given with --faa does not seem to be a valid directory or file. Please check.')

def check_gbk_path(gbk_path, samplename):
    if(os.path.isdir(gbk_path)):
        path_list = list(pathlib.Path(gbk_path).rglob(f"*{samplename}*.gbff")) #gbk or gbff[gG][bB][fFkK]
        if (len(path_list)>1):
            sys.exit(f'AMPcombi interrupted: There is more than one .gbk/.gbff file for {samplename} in the folder given with --gbk')
        elif(not path_list):
            sys.exit(f'AMPcombi interrupted: There is no .gbk/.gbff file containing {samplename} in the folder given with --gbk')
        return path_list[0]
    elif(os.path.isfile(gbk_path)):
        return gbk_path
    else:
        sys.exit(f'AMPcombi interrupted: The input given with --gbk does not seem to be a valid directory or file. Please check.')

def check_interpro_path(interpro_path, samplename):
    if interpro_path != None:
        if os.path.isdir(interpro_path):
            path_list = list(pathlib.Path(interpro_path).rglob(f"*{samplename}*.tsv"))
            if len(path_list) > 1:
                sys.exit(f'AMPcombi interrupted: There is more than one .tsv file for {samplename} in the folder given with --interproscan_output')
            elif not path_list:
                sys.exit(f'AMPcombi interrupted: There is no .tsv file containing {samplename} in the folder given with --interproscan_output')
            return path_list[0]
        elif os.path.isfile(interpro_path):
            return interpro_path
    else:
        print("No InterproScan files were provided. Workflow continuing ....")
        return None
  
def check_ref_database(database, database_dir, threads):
    valid_databases = ['DRAMP', 'APD', 'UniRef100']
    local_db_path = f'amp_{database}_database'
    #if database not in valid_databases:
    #    sys.exit(f"AMPcombi interrupted: {database} is not a valid AMP database. Choose from {valid_databases}.")    
    if ((database in valid_databases) and (database_dir != None)):
        if (os.path.exists(database_dir)):
            db = database_dir
            print(f'<--AMP_database> = {db} is found and will be used')
            return db
        else: 
            sys.exit(f'AMPcombi interrupted: Please check the path again to the database folder provided in --amp_database_dir {db}')
    elif((database in valid_databases) and (database_dir is None)):
        os.makedirs(local_db_path, exist_ok=True)
        db = local_db_path
        #make sure the directory is not empty 
        if not os.listdir(db):
            print(f'<--AMP_database> = ${database} will be downloaded and used')
            download_ref_db(db, database, threads)
            return db
        else:
            print(f'The {db} contains files and so database will not be redownloaded.')
            return db
    elif((database in valid_databases) and (database_dir is None) and (os.path.exists(local_db_path))):
        print(f'<--AMP_database> = {database} is already downloaded and will be reused')
        db = local_db_path
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
