#!/bin/python3

# TITLE: Concatenate AMPcombi sample summaries

import pandas as pd
import os
import sys

#########################################
# FUNCTION: CONCATENATE AMPCOMBI SUMMARIES
#########################################
# concatenates AMPcombi files into the final complete summary
def concatenate_summaries(file_path):
    if not os.path.exists(file_path):
        sys.exit(f'AMPcombi interrupted: File path {file_path} does not exist. Please correct the file path')
    else:
        # all paths exists
        dfs = []
        df = pd.read_csv(file_path, delimiter='\t')
        # append the dataframe to the list
        dfs.append(df)                                        
        # concatenate the dataframes in the list
        concatenated_dfs = pd.concat(dfs, ignore_index=True, sort=False)
        # get all unique column names from all dataframes
        unique_columns = set()
        for df in dfs:
            unique_columns.update(df.columns)
        # create a new dataframe with all possible columns
        new_df = pd.DataFrame(columns=list(unique_columns))
        # concatenate the original dataframe with the new dataframe
        concatenated_df = pd.concat([new_df, concatenated_dfs], ignore_index=True, sort=False)
        # remove duplicates rows
        concatenated_df = concatenated_df.drop_duplicates()
        # Save the concatenated dataframe to a file
        concatenated_df.to_csv('Ampcombi_summary.tsv', index=False, sep='\t')        
