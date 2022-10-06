#!/usr/bin/env python3

import os
import argparse
import warnings
# import functions from sub-scripts to main:
from reformat_tables import *
from amp_fasta import *
from check_input import *
from amp_database import *
from print_header import *
from contextlib import redirect_stdout

# Define input arguments:
parser = argparse.ArgumentParser(prog = 'ampcombi', formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=('''\
    .............................................................................
                                    *AMP-combi*
    .............................................................................
                This tool parses the results of amp prediction tools 
    and aligns the hits against reference databases for functional classification.
            For detailed usage documentation please refer to <github_repo>
    .............................................................................'''),
                                epilog='''Thank you for running AMP-combi!''',
                                add_help=True)

#TODO: Find out if parser and additional argument should be before or part of main()
#TODO: Consider adding a parameters to keep the contigs with all non-amps and not to filter them out 
# For stand a lone tool JFY suggests having a similar tool function as HAMRONIZATION, i.e. seperate the tools into different sub-modules
# -- subtools =2: parse and summaruze (classify) for parsing --amppir --amplify --n

parser.add_argument("--amp_results", dest="amp", nargs='?', help="Enter the path to the folder that contains the different tool's output files in sub-folders named by sample name. \n If paths are to be inferred, sub-folders in this results-directory have to be organized like '/amp_results/toolsubdir/samplesubdir/tool.sample.filetype' \n (default: %(default)s)",
                    type=str, default="./test_files/")
parser.add_argument("--sample_list", dest="samples", nargs='?', help="Enter a list of sample-names, e.g. ['sample_1', 'sample_2', 'sample_n']. \n If not given, the sample-names will be inferred from the folder structure",
                    type=list, default=[])
parser.add_argument("--path_list", dest="files", nargs='?', help="Enter the list of paths to the files to be summarized as a list of lists, e.g. [['path/to/my/sample1.ampir.tsv', 'path/to/my/sample1.amplify.tsv'], ['path/to/my/sample2.ampir.tsv', 'path/to/my/sample2.amplify.tsv']]. \n If not given, the file-paths will be inferred from the folder structure",
                    type=list, default=[])
parser.add_argument("--outdir", dest="out", help="Enter the name of the output directory \n (default: %(default)s)",
                    type=str, default="./ampcombi_results/")
parser.add_argument("--cutoff", dest="p", help="Enter the probability cutoff for AMPs \n (default: %(default)s)",
                    type=int, default=0)
parser.add_argument("--faa_folder", dest="faa", help="Enter the path to the folder containing the reference .faa files. Filenames have to contain the corresponding sample-name, i.e. sample_1.faa \n (default: %(default)s)",
                    type=str, default='./test_faa/')
parser.add_argument("--tooldict", dest="tools", help="Enter a dictionary of the AMP-tools used with their output file endings (as they appear in the directory tree), \n Tool-names have to be written as in default:\n default={'ampir':'ampir.tsv', 'amplify':'amplify.tsv', 'macrel':'macrel.tsv', 'hmmer_hmmsearch':'hmmsearch.txt', 'ensembleamppred':'ensembleamppred.txt'}",
                    type=dict, default={'ampir':'ampir.tsv', 'amplify':'amplify.tsv', 'macrel':'macrel.tsv', 'neubi':'neubi.fasta', 'hmmer_hmmsearch':'hmmsearch.txt', 'ensembleamppred':'ensembleamppred.txt'})
parser.add_argument("--amp_database", dest="ref_db", nargs='?', help="Enter the path to the folder containing the reference database files (.fa and .tsv); a fasta file and the corresponding table with functional and taxonomic classifications. \n (default: DRAMP database)",
                    type=str, default=None)
parser.add_argument("--log", dest="log_file", nargs='?', help="Silences the standardoutput and captures it in a log file)",
                    type=bool, default=False)

#parser.add_argument('--help', default=argparse.SUPPRESS, nargs='?', help='Show this help message and exit.')

# print help message for user
#parser.print_help()

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
database = args.ref_db

# additional variables
# extract list of tools from input dictionary. If not given, default dict contains all possible tools
tools = [key for key in tooldict]
# extract list of tool-output file-endings. If not given, default dict contains default endings.
fileending = [val for val in tooldict.values()]

# create output directory
os.makedirs(outdir, exist_ok=True)

# supress panda warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#########################################
# MAIN FUNCTION
#########################################
def main_workflow():
    # print_header()
    print_header()
    # check input parameters
    check_input_complete(path, samplelist_in, filepaths_in, tools)
    # check input sample-list and create sample-list if input empty
    samplelist = check_samplelist(samplelist_in, tools, path)
    # check input filepaths and create list of list of filepaths per sample if input empty
    filepaths = check_pathlist(filepaths_in, samplelist, fileending, path)
    # check amp_ref_database filepaths and create a directory if input empty
    db = check_ref_database(database, outdir)

    # generate summary for each sample
    amp_faa_paths = []
    create_diamond_ref_db(db)
    for i in range(0, len(samplelist)):
        main_list = []
        print('\n ########################################################## ')
        print(f'Processing AMP-files from sample: {samplelist[i]}')
        os.makedirs(outdir + '/'+ samplelist[i], exist_ok=True)
        # fill main_list with tool-output filepaths for sample i
        read_path(main_list, filepaths[i], p, tooldict, faa_path, samplelist[i])
        # use main_list to create the summary file for sample i
        summary_df = summary(main_list, samplelist[i], faa_path, outdir)
        # Generate the AMP-faa.fasta for sample i
        out_path = outdir+ '/'+samplelist[i] +'/'+samplelist[i]+'_amp.faa'
        faa_name = faa_path+samplelist[i]+'.faa'
        amp_fasta(summary_df, faa_name, out_path)
        amp_faa_paths.append(out_path)
        print(f'The fasta containing AMP sequences for {samplelist[i]} was saved to {outdir}/{samplelist[i]}/ \n')
        amp_matches = outdir + '/'+samplelist[i] +'/'+samplelist[i]+'_diamond_matches.txt'
        print(f'The diamond alignment for {samplelist[i]} in process....')
        diamond_df = diamond_alignment(db, amp_faa_paths, amp_matches)
        print(f'The diamond alignment for {samplelist[i]} was saved to {outdir}/{samplelist[i]}/.')
        # Merge summary_df and diamond_df
        complete_summary_df = pd.merge(summary_df, diamond_df, on = 'contig_id', how='left')
        complete_summary_df.to_csv(outdir +'/'+samplelist[i] +'/'+samplelist[i]+'_ampcombi.csv', sep=',')
        print(f'The summary file for {samplelist[i]} was saved to {outdir}/{samplelist[i]}/.')
        
def main():
    if args.log_file == True:
        with open(f'{outdir}/ampcombi.log', 'w') as f:
            with redirect_stdout(f):
                main_workflow()
    else: main_workflow()
    

if __name__ == "__main__":
    main()