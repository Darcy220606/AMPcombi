#!/bin/python3

import os

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