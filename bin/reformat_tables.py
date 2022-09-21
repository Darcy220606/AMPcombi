#!/bin/python3

# TITLE: Reformat the AMP output tables

import os
import pandas as pd
import argparse
from Bio import SeqIO

# Define input arguments:
parser = argparse.ArgumentParser()

#TODO: add -h as argument and only print helpmessages if --help/-h is called
#TODO: print a nice AMPcombi header
parser.add_argument("--amp_results", dest="amp", nargs='?', help="enter the path to the folder that contains the different tool's output files in sub-folders named by sample name. \n If paths are to be inferred, sub-folders in this results-directory have to be organized like '/amp_results/toolsubdir/samplesubdir/tool.sample.tsv'",
                    type=str, default="../amp_results/")
parser.add_argument("--sample_list", dest="samples", nargs='?', help="enter a list of sample-names, e.g. ['sample_1', 'sample_2', 'sample_n']. \n If not given, the sample-names will be inferred from the folder structure",
                    type=list, default=[])
parser.add_argument("--path_list", dest="files", nargs='?', help="enter the list of paths to the files to be summarized as a list of lists, e.g. [['path/to/my/sample1.ampir.tsv', 'path/to/my/sample1.amplify.tsv'], ['path/to/my/sample2.ampir.tsv', 'path/to/my/sample2.amplify.tsv']]. \n If not given, the file-paths will be inferred from the folder structure",
                    type=list, default=[])
parser.add_argument("--outdir", dest="out", help="enter the name of the output directory",
                    type=str, default="../amp_summary/")
parser.add_argument("--cutoff", dest="p", help="enter the probability cutoff for AMPs",
                    type=int, default=0.5)
parser.add_argument("--faa", dest="faa", help="enter the path to the folder containing the reference .faa files. Filenames have to contain the corresponding sample-name, i.e. sample_1.faa",
                    type=str, default='../test_faa/')

# print help message for user
parser.print_help()

# get command line arguments
args = parser.parse_args()

# assign input arguments to variables
path = args.amp
samplelist = args.samples
filepaths = args.files
outdir = args.out
p = args.p

# additional variables
# TODO: flexibilize this input (also see below): add this to input args, user can provide a dict of 'tool':'tool-fileending'
tools = ['ampir', 'amplify', 'hmmer_hmmsearch', 'macrel'] # to include more tools, add names here
fileending = ['ampir.tsv', 'amplify.tsv', 'macrel.tsv', 'hmmsearch.txt'] # add endings of new tools, FLEXIBILIZE THIS

# create output directory
os.makedirs(outdir, exist_ok=True)

# TODO: Check input: either --amp-results directory OR --path-list has to be given
# TODO: check function should print INFO to screen

#########################################
# GENERATE LIST OF AMP-FILE-LISTS
#########################################
# list has to be processed per sample (to create summary per sample)
# list paths to target files per sample: list = [[pathlist_sample_1], ..., [pathlist_sample_n]]
toollist = []
pathlist = []

if(samplelist==[]):
    print('<--sample-list> was not given, sample names will be inferred from directory names')
    for dirpath, subdirs, files in os.walk(path):
        for dir in subdirs:
            if (dir in tools):
                toollist.append(dir)
            else: 
                samplelist.append(dir)
    samplelist = list(set(samplelist))


if(filepaths==[]):
    print('<--path-list> was not given, paths to AMP-results-files will be inferred')
    for sample in samplelist:
        for dirpath, subdirs, files in os.walk(path):
            for file in files:
                if ((sample in dirpath)&((list(filter(file.endswith, fileending))!=[]))):
                    pathlist.append(dirpath+'/'+file)
        filepaths.append(pathlist)
        #reset the pathlist for next sample
        pathlist = []

#########################################
# FUNCTIONS: READ TOOL OUTPUT TO DF
#########################################

#########################################
    #  AMP_ampir
#########################################
def ampir(path): 
    # Dictionary to rename columns
    ampir_dict = {'seq_name':'contig_id', 'seq_aa':'seq_aa', 'prob_AMP':'prob_ampir'}
    # read file as df and rename columns
    ampir_df = pd.read_csv(path, sep='\t').rename(columns=ampir_dict) 
    # cut contig_id to remove extra information added by tool
    ampir_df['contig_id']=ampir_df['contig_id'].apply(lambda x: x.split()[0])
    # apply probability cutoff
    ampir_df = ampir_df[(ampir_df['prob_ampir']>=p)]
    return ampir_df[['contig_id', 'prob_ampir']]

#########################################
    #  AMP_amplify
#########################################
def amplify(path):
    amplify_dict = {'Sequence_ID':'contig_id', 'Sequence':'seq_aa', 'Length':'length', 'Charge':'charge', 'Probability_score':'prob_amplify', 'AMPlify_log_scaled_score':'log_score', 'Prediction':'prediction'}
    amplify_df = pd.read_csv(path, sep='\t').rename(columns=amplify_dict).dropna()
    # apply probability cutoff
    amplify_df = amplify_df[(amplify_df['prob_amplify']>=p)]
    return amplify_df[['contig_id', 'prob_amplify']]

#########################################
    #  AMP_macrel
#########################################
def macrel(path):
    macrel_dict = {'Access':'contig_id', 'Sequence':'seq_aa', 'AMP_family':'amp_family', 'AMP_probability':'prob_macrel', 'Hemolytic':'hemolytic', 'Hemolytic_probability':'prob_hemo'}
    #set header to second row to skip first line starting with #
    macrel_df = pd.read_csv(path, sep='\t', header=[1]).rename(columns=macrel_dict)
    # apply probability cutoff
    macrel_df = macrel_df[(macrel_df['prob_macrel']>=p)]
    return macrel_df[['contig_id', 'prob_macrel']]

#########################################
    #  AMP_hmmsearch
#########################################
def hmmsearch(path):
    hmmer_dict = {'level_0':'evalue_hmmer', 'level_1':'score_hmmer', 'level_2':'bias', 'level_3':'eval_domain', 'level_4':'score_domain', 'level_5':'bias_domain', 'level_6':'exp_dom', '-------':'N_dom', '------':'contig_id'}
    hmmer_df = pd.read_table(path, delim_whitespace=True, header=[15]).reset_index().rename(columns=hmmer_dict)
    hmmer_df = hmmer_df.drop(hmmer_df.iloc[:,9:17], axis=1) #drop unnecessary columns
    for index, row in hmmer_df.iterrows():
        if (row.str.contains('Domain').any()):              #identify index of first row with 'Domain'
            i = index
            break
    hmmer_df = hmmer_df[hmmer_df.index<i]                   #only keep rows previous to index i
    return hmmer_df[['contig_id', 'evalue_hmmer']]

#########################################
# FUNCTION: READ DFs PER SAMPLE 
#########################################
# For one sample: parse filepaths and read files to dataframes, create list of dataframes
# TODO: for more flexible input: toollist and fileendings could be given in a dict and then create the lists (for now hardcoded)
def read_path(df_list, file_list):
    for path in file_list:
        if(path.endswith(fileending[0])):
            print('found ampir file')
            df_list.append(ampir(path))
        elif(path.endswith(fileending[1])):
            print('found amplify file')
            df_list.append(amplify(path))
        elif(path.endswith(fileending[2])):
            print('found macrel file')
            df_list.append(macrel(path))
        elif(path.endswith(fileending[3])):
            print('found hmmersearch file')
            df_list.append(hmmsearch(path))
        else:
            print('No AMP-output-files could be found with the given path. \n Please check your file paths and file endings or use the <--path-list> command')
            break

#########################################
# FUNCTION: MERGE DATAFRAMES
#########################################
# merge dataframes from list to summary output per sample
def summary(df_list, samplename):
    merge_df = pd.DataFrame(columns=['contig_id'])
    for df in df_list:
        merge_df = pd.merge(merge_df, pd.DataFrame(df) , how='outer', on='contig_id')
    merge_df = merge_df.fillna(0)
    # TODO: sort by highest p-values AND/OR put results found by most tools first before output
    merge_df.to_csv(outdir+'/'+samplename+'_AMPsummary.csv', sep=',')
    #return merge_df

#########################################
# FUNCTION: ADD AA-SEQUENCE
#########################################
# TODO: function to add the amino-acid sequence extracted from the faa


#########################################
# MAIN FUNCTION
#########################################
if __name__ == "__main__":
    #print_header()
    main_list = []
    for i in range(0, len(samplelist)):
        print(f'Processing AMP-files from sample: {samplelist[0]}')
        read_path(main_list, filepaths[i])
        summary(main_list, samplelist[i])
        print(f'The summary file for {samplelist[i]} was saved to {outdir}')
        # reset main_list for next sample summary
        main_list=[]
