#!/bin/python3

# TITLE: Concatenate AMPcombi sample summaries

import pandas as pd
import os
import sys

#########################################
# FUNCTION: CONCATENATE AMPCOMBI SUMMARIES
#########################################
# concatenates AMPcombi files into the final complete summary
def concatenate_summaries(file_list):
    max_columns = 0
    max_columns_df = None
    for file_path in file_list:
        if not os.path.exists(file_path):
            sys.exit(f'AMPcombi interrupted: File path {file_path} does not exist. Please correct the file path')
        else:
            # all paths exists
            df = pd.read_csv(file_path, delimiter='\t')
            # Check if this dataframe has more columns than the current max
            if df.shape[1] > max_columns:
                max_columns = df.shape[1]
                max_columns_df = df
    
    if max_columns_df is None:
        sys.exit('AMPcombi interrupted: No valid dataframes found in the file list')

    # Create a new dataframe with column headers from the dataframe with the most columns
    concatenated_df = pd.DataFrame(columns=max_columns_df.columns)

    # Concatenate the dataframes
    for file_path in file_list:
        df = pd.read_csv(file_path, delimiter='\t')
        concatenated_df = pd.concat([concatenated_df, df], ignore_index=True, sort=False)

    # Remove duplicate rows
    concatenated_df = concatenated_df.drop_duplicates()

    # Save the concatenated dataframe to a file
    concatenated_df.to_csv('Ampcombi_summary.tsv', index=False, sep='\t')