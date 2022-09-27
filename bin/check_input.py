#!/bin/python3

import os
from amp_database import download_DRAMP

# TODO: Check input: either --amp-results directory OR --path-list has to be given
# TODO: check function should print INFO to screen

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
                    if ((sample in dirpath)&((list(filter(file.endswith, fileending))!=[]))):
                        pathlist.append(dirpath+'/'+file)
            filepaths.append(pathlist)
        return filepaths
    else:
        return filepaths

def check_dfshape(df1, df2):
    if (df1.shape[0] != df2.shape[0]):
        print(f'ERROR: different row number in tool output and faa file. Ensembleamppred output could not be included in the summary')
        return False
    else: 
        return True

def check_ref_database(database, outdir):
    if(database==[]):
        print('<--AMP_database> was not given, the current DRAMP general-AMP database will be downloaded and used')
        database = os.path.join(outdir, r'amp_ref_database')
        os.makedirs(database, exist_ok=True)
        db = database
        download_DRAMP(db)
        return db
    else:
        db = database
        return db
    #else:
        #if os.path.exists(database):
            #db = database
            #return db
            #else:
            # #if not os.path.exists(database):
            #raise ValueError('Reference amp database does not exist, please check the path.')
