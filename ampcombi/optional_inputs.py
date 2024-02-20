#!/bin/python3

# TITLE: Optional inputs to expand sample details

import pandas as pd
from Bio import SeqIO
import os

#########################################
# FUNCTION: Add sample metadata
#########################################
def sample_metadata_addition(merged_df, smetadata):
    """
    Adds the sample metadata only when turned on.
    Important to note: the first column must have the sample name and the should be tab seperated.
    """
    if smetadata == None :
        return merged_df
    else:
        # add the samples metadata 
        metadata_df = pd.read_csv(smetadata, sep='\t')
        metadata_df.rename(columns={metadata_df.columns[0]: 'name'}, inplace=True)
        # merge it to the df using sample name 'name' as common
        df1 = merged_df.merge(metadata_df, on='name', how='left')
        return df1

#########################################
# FUNCTION: Add contig metadata
#########################################
def contig_metadata_addition(merged_df, cmetadata):
    """
    Adds the contig metadata only when turned on.
    Important to note: the first column must have the sample name and the second name must have contig_id
                       the table should be tab seperated.
    """
    if cmetadata == None :
        return merged_df
    else:
        # add the samples metadata 
        metadata_df = pd.read_csv(cmetadata, sep='\t')
        # rename first column : must contain the sample names
        metadata_df.rename(columns={metadata_df.columns[0]: 'name'}, inplace=True)
        # rename second column _ must contain the contig names
        metadata_df.rename(columns={metadata_df.columns[1]: 'contig_name'}, inplace=True)
        # ensure that the dataframe is not empty or else KeyError arises
        if not merged_df.empty:
            # remove duplicate entries
            merged_df.drop_duplicates(inplace=True)
            # merge it to the df using sample name 'name' as common
            df2 = pd.merge(merged_df, metadata_df, on=['contig_name', 'name'], how='left')
        else:
            df2 = merged_df
        return df2
