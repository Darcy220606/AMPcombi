#!/usr/bin/env python3

import os
import argparse
import warnings
import pandas as pd
from contextlib import redirect_stdout

# Define input arguments:
parser = argparse.ArgumentParser(prog = 'summarize', formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=('''\
    .............................................................................
                                    *ampcombi-summarize*
    .............................................................................
                This sub-tool summarizes the results of ampcombi-parse module. 
                For detailed usage documentation please refer to <github_repo>
    .............................................................................'''),
                                epilog='''Thank you for running AMPcombi!''',
                                add_help=True)

parser.add_argument("--summary_list", dest="summary", nargs='*', action='append', help="Enter a list that contains the path to the ampcombi summary files. \n  E.g. '/Sample_1/sample_1_ampcombi.csv' \n (default: %(default)s)",
                    default=[])

# get command line arguments
args = parser.parse_args()

# assign input arguments to variables
summaries = args.summary

# supress panda warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

#########################################
# MAIN FUNCTION
#########################################
def summarize(summaries):
    print(f'\n Starting AMPcombi-summarize')
    # initiate a final_summary dataframe to concatenate each new sample-summary
    complete_summary_df = pd.DataFrame([])
    # concatenate the sample summary to the complete summary and overwrite it
    for s in summaries:
        complete_summary_df = pd.concat([complete_summary_df, pd.read_csv(s)])
        #print(sum.shape)
    complete_summary_df.to_csv('AMPcombi_summary.csv', sep=',', index=False)
    print(f'\n FINISHED: The AMPcombi_summary.csv file was saved to your current working directory.')
    
def main():
    if (os.path.exists(f'ampcombi.log') and args.log_file == True):
        with open(f'ampcombi.log', 'a') as f:
            with redirect_stdout(f):
                summarize(summaries)
    else: summarize(summaries)
    
if __name__ == "__main__":
    main()