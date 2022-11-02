#!/usr/bin/env python3

import os
import argparse
import warnings
from contextlib import redirect_stdout
from version import __version__
import json
import os.path
# import functions from sub-scripts to main:
from reformat_tables import *
from amp_fasta import *
from check_input import *
from amp_database import *
from print_header import *
from visualise_complete_summary import *

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

parser.add_argument("--amp_results", dest="amp", nargs='?', help="Enter the path to the folder that contains the different tool's output files in sub-folders named by sample name. \n If paths are to be inferred, sub-folders in this results-directory have to be organized like '/amp_results/toolsubdir/samplesubdir/tool.sample.filetype' \n (default: %(default)s)",
                    type=str, default='./test_files/')
parser.add_argument("--sample_list", dest="samples", nargs='*', help="Enter a list of sample-names, e.g. sample_1 sample_2 sample_n. \n If not given, the sample-names will be inferred from the folder structure",
                    default=[])
parser.add_argument("--path_list", dest="files", nargs='*', action='append', help="Enter the list of paths to the files to be summarized as a list of lists, e.g. --path_list path/to/my/sample1.ampir.tsv path/to/my/sample1.amplify.tsv --path_list path/to/my/sample2.ampir.ts path/to/my/sample2.amplify.tsv. \n If not given, the file-paths will be inferred from the folder structure",
                    default=[])
parser.add_argument("--cutoff", dest="p", help="Enter the probability cutoff for AMPs \n (default: %(default)s)",
                    type=int, default=0)
parser.add_argument("--faa_folder", dest="faa", help="Enter the path to the folder containing the reference .faa files. Filenames have to contain the corresponding sample-name, i.e. sample_1.faa \n (default: %(default)s)",
                    type=str, default='./test_faa/')
parser.add_argument("--tooldict", dest="tools", help="Enter a dictionary of the AMP-tools used with their output file endings (as they appear in the directory tree), \n Tool-names have to be written as in default:\n default={'ampir':'ampir.tsv', 'amplify':'amplify.tsv', 'macrel':'macrel.tsv', 'hmmer_hmmsearch':'hmmsearch.txt', 'ensembleamppred':'ensembleamppred.txt'}",
                    type=str, default='{"ampir":"ampir.tsv", "amplify":"amplify.tsv", "macrel":"macrel.tsv", "neubi":"neubi.fasta", "hmmer_hmmsearch":"hmmsearch.txt", "ensembleamppred":"ensembleamppred.txt"}')
parser.add_argument("--amp_database", dest="ref_db", nargs='?', help="Enter the path to the folder containing the reference database files (.fa and .tsv); a fasta file and the corresponding table with functional and taxonomic classifications. \n (default: DRAMP database)",
                    type=str, default=None)
parser.add_argument("--complete_summary", dest="complete", nargs='?', help="Concatenates all sample summaries to one final summary and outputs both csv and interactive html files",
                    type=bool, default=False)
parser.add_argument("--log", dest="log_file", nargs='?', help="Silences the standard output and captures it in a log file)",
                    type=bool, default=False)
parser.add_argument("--threads", dest="cores", nargs='?', help="Changes the threads used for DIAMOND alignment (default: %(default)s)",
                    type=bool, default='4')
parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

# get command line arguments
args = parser.parse_args()

# assign input arguments to variables
path = args.amp
samplelist_in = args.samples
filepaths_in = args.files
p = args.p
faa_path = args.faa
tooldict = json.loads(args.tools)
database = args.ref_db
complete_summary = args.complete
threads = args.cores

# additional variables
# extract list of tools from input dictionary. If not given, default dict contains all possible tools
tools = [key for key in tooldict]
# extract list of tool-output file-endings. If not given, default dict contains default endings.
fileending = [val for val in tooldict.values()]

# supress panda warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#########################################
# MAIN FUNCTION
#########################################
def main_workflow():
    # print AMPcombi header
    print_header()
    # check input sample-list and create sample-list if input empty
    samplelist = check_samplelist(samplelist_in, tools, path)
    # check input parameters
    check_input_complete(path, samplelist, filepaths_in, tools)
    # check input filepaths and create list of list of filepaths per sample if input empty
    filepaths = check_pathlist(filepaths_in, samplelist, fileending, path)
    # check amp_ref_database filepaths and create a directory if input empty
    db = check_ref_database(database)
    # initiate a final_summary dataframe to concatenate each new sample-summary
    if (complete_summary):
        complete_summary_df = pd.DataFrame([])

    # generate summary for each sample
    amp_faa_paths = []
    create_diamond_ref_db(db,threads)
    for i in range(0, len(samplelist)):
        main_list = []
        print('\n ########################################################## ')
        print(f'Processing AMP-files from sample: {samplelist[i]}')
        os.makedirs(samplelist[i], exist_ok=True)
        # fill main_list with tool-output filepaths for sample i
        read_path(main_list, filepaths[i], p, tooldict, faa_path, samplelist[i])
        # get the path to the samples' corresponding faa file
        faa_name = check_faa_path(faa_path, samplelist[i])
        # use main_list to create the summary file for sample i
        summary_df = summary(main_list, samplelist[i], faa_name)
        # Generate the AMP-faa.fasta for sample i
        out_path = samplelist[i] +'/'+samplelist[i]+'_amp.faa'
        amp_fasta(summary_df, faa_name, out_path)
        amp_faa_paths.append(out_path)
        print(f'The fasta containing AMP sequences for {samplelist[i]} was saved to {samplelist[i]}/ \n')
        amp_matches = samplelist[i] +'/'+samplelist[i]+'_diamond_matches.txt'
        print(f'The diamond alignment for {samplelist[i]} in progress ....')
        diamond_df = diamond_alignment(db, amp_faa_paths, amp_matches, threads)
        print(f'The diamond alignment for {samplelist[i]} was saved to {samplelist[i]}/.')
        # Merge summary_df and diamond_df
        sample_summary_df = pd.merge(summary_df, diamond_df, on = 'contig_id', how='left')
        # Insert column with sample name on position 0
        sample_summary_df.insert(0, 'name', samplelist[i])
        # Write sample summary into sample output folder
        sample_summary_df.to_csv(samplelist[i] +'/'+samplelist[i]+'_ampcombi.csv', sep=',', index=False)
        print(f'The summary file for {samplelist[i]} was saved to {samplelist[i]}/.')
        if (complete_summary):
        # concatenate the sample summary to the complete summary and overwrite it
            complete_summary_df = pd.concat([complete_summary_df, sample_summary_df])
            complete_summary_df.to_csv('AMPcombi_summary.csv', sep=',', index=False)
            html_generator() 
        else: 
            continue
    if (complete_summary):
        print(f'\n FINISHED: The AMPcombi_summary.csv and AMPcombi_summary.html files were saved to your current working directory.')
    else: 
        print(f'\n FINISHED: AMPcombi created summaries for all input samples.')

def main():
    if (args.log_file == True and not os.path.exists('ampcombi.log')):
        with open(f'ampcombi.log', 'w') as f:
            with redirect_stdout(f):
                main_workflow()
    elif(args.log_file == True and os.path.exists('ampcombi.log')):
        with open(f'ampcombi.log', 'a') as f:
            with redirect_stdout(f):
                main_workflow()
    else: main_workflow()
    

if __name__ == "__main__":
    main()