#!/bin/python3

import os
import argparse
# import functions from sub-scripts to main:
from reformat_tables import *
from amp_fasta import *
from check_input import *


# Define input arguments:
parser = argparse.ArgumentParser()

#TODO: add -h as argument and only print helpmessages if --help/-h is called
#TODO: print a nice AMPcombi header
#TODO: find out if parser and additional argument should be before or part of main()

parser.add_argument("--amp_results", dest="amp", nargs='?', help="enter the path to the folder that contains the different tool's output files in sub-folders named by sample name. \n If paths are to be inferred, sub-folders in this results-directory have to be organized like '/amp_results/toolsubdir/samplesubdir/tool.sample.filetype'",
                    type=str, default="../amp_results/")
parser.add_argument("--sample_list", dest="samples", nargs='?', help="enter a list of sample-names, e.g. ['sample_1', 'sample_2', 'sample_n']. \n If not given, the sample-names will be inferred from the folder structure",
                    type=list, default=[])
parser.add_argument("--path_list", dest="files", nargs='?', help="enter the list of paths to the files to be summarized as a list of lists, e.g. [['path/to/my/sample1.ampir.tsv', 'path/to/my/sample1.amplify.tsv'], ['path/to/my/sample2.ampir.tsv', 'path/to/my/sample2.amplify.tsv']]. \n If not given, the file-paths will be inferred from the folder structure",
                    type=list, default=[])
parser.add_argument("--outdir", dest="out", help="enter the name of the output directory",
                    type=str, default="../amp_summary/")
parser.add_argument("--cutoff", dest="p", help="enter the probability cutoff for AMPs",
                    type=int, default=0)
parser.add_argument("--faa_folder", dest="faa", help="enter the path to the folder containing the reference .faa files. Filenames have to contain the corresponding sample-name, i.e. sample_1.faa",
                    type=str, default='../test_faa/')
parser.add_argument("--tooldict", dest="tools", help="enter a dictionary of the AMP-tools used with their output file endings (as they appear in the directory tree), \n Tool-names have to be written as in default:\n default={'ampir':'ampir.tsv', 'amplify':'amplify.tsv', 'macrel':'macrel.tsv', 'hmmer_hmmsearch':'hmmsearch.txt', 'ensembleamppred':'ensembleamppred.txt'}",
                    type=dict, default={'ampir':'ampir.tsv', 'amplify':'amplify.tsv', 'macrel':'macrel.tsv', 'neubi':'neubi.fasta' ,'hmmer_hmmsearch':'hmmsearch.txt', 'ensembleamppred':'ensembleamppred.txt'})

# print help message for user
parser.print_help()

# get command line arguments
args = parser.parse_args()

# assign input arguments to variables
path = args.amp
samplelist_in = args.samples
filepaths_in = args.files
outdir = args.out
p = args.p
faa_path = args.faa
tooldict = args.tools

# additional variables
# extract list of tools from input dictionary. If not given, default dict contains all possible tools
tools = [key for key in tooldict]
# extract list of tool-output file-endings. If not given, default dict contains default endings.
fileending = [val for val in tooldict.values()]
# check input sample-list and create sample-list if input empty
samplelist = check_samplelist(samplelist_in, tools, path)
# check input filepaths and create list of list of filepaths per sample if input empty
filepaths = check_pathlist(filepaths_in, samplelist, fileending, path)

# create output directory
os.makedirs(outdir, exist_ok=True)

#########################################
# MAIN FUNCTION
#########################################
if __name__ == "__main__":
    #print_header()
    amp_faa_paths = []
    for i in range(0, len(samplelist)):
        main_list = []
        print(f'Processing AMP-files from sample: {samplelist[i]}')
        # fill main_list with tool-output filepaths for sample i
        read_path(main_list, filepaths[i], p, tooldict, faa_path, samplelist[i])
        # use main_list to create the summary file for sample i
        summary_df = summary(main_list, samplelist[i], faa_path, outdir)
        print(f'The summary file for {samplelist[i]} was saved to {outdir}')
        # Generate the AMP-faa.fasta for sample i
        out_path = outdir+samplelist[i]+'_amp.faa'
        faa_name = faa_path+samplelist[i]+'.faa'
        amp_fasta(summary_df, faa_name, out_path)
        amp_faa_paths.append(out_path)
        #call: check download
        #call: function that runs Diamond.bash
        #call: read Diamond output and add to summary
        print(f'The fasta containing AMP sequences for {samplelist[i]} was saved to {outdir} \n')
    print('Your AMPcombi summaries are now available in the output folder!')