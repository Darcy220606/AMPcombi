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

##########################################
#  FUNCTION: Parse interproscan file 
##########################################
def parse_interproscan(summary_df, interpro_name, interpro_filter_values):
    """
    Merges the ouput from interproscan only when turned on.
    IMPORTANT: It parses ONLY the output from interproscan using:
    ' --applications PANTHER,ProSiteProfiles,ProSitePatterns,Pfam '
    If more applications are used please adapt the column headers here.
    The 'filter_values' are values defined by the user that should be in the description for the hit to be removed.
    """
    if interpro_name != None :
        print(f'Adding ingterproscan and removing hits with {interpro_filter_values} ...')
        # read in the tsv and add column headers to the output
        interpro_headers = ['contig_id','MD5sums','lengths','analysis','accession','description','start','stop','evalue','exists','date','interpro_accession','interpro_description','unknown_1','unknown_2']
        interpro_df = pd.read_csv(interpro_name, sep='\t', names=interpro_headers)
        # remove unnecessary interpro headers
        interpro_noncolumns = ['MD5sums','lengths','analysis','start','stop','evalue','exists','date','unknown_1','unknown_2']
        interpro_df.drop(interpro_noncolumns, axis=1, inplace=True)
        # group the columns by unique cds_id using aggregate
        concatenate_columns = ['accession','description','interpro_accession','interpro_description']
        interpro_df_grouped = interpro_df.groupby('contig_id')[concatenate_columns].agg(', '.join).reset_index() # contig_id means here cds_ID but that changes only with gbk file
        # merge interpro results
        ampcombi_df_interpro = pd.merge(summary_df, interpro_df_grouped, on='contig_id', how='left')
        # remove ribosomal proteins or any value given in a list (def non AMPs)
        values_to_remove = interpro_filter_values.split(",")
        ampcombi_df_interpro['description'] = ampcombi_df_interpro['description'].astype(str)
        ampcombi_df_interpro_filt = ampcombi_df_interpro[~ampcombi_df_interpro['description'].str.contains('|'.join(values_to_remove), case=False)]        
        return ampcombi_df_interpro_filt
    else:
        return summary_df
