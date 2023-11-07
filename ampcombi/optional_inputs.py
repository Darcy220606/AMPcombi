#!/bin/python3

# TITLE: Optional inputs to expand sample details

import pandas as pd
from Bio import SeqIO
import os

#########################################
# FUNCTION: Add sample metadata
#########################################
def metadata_addition(merged_df, metadata):
    """
    Adds the sample metadata only when turned on.
    Important to note: the first column must have the sample name and the should be tab seperated.
    """
    if metadata == None :
        return merged_df
    else:
        # add the samples metadata 
        metadata_df = pd.read_csv(metadata, sep='\t')
        metadata_df.rename(columns={metadata_df.columns[0]: 'name'}, inplace=True)
        # merge it to the df using sample name 'name' as common
        df1 = merged_df.merge(metadata_df, on='name', how='left')
        return df1
    
