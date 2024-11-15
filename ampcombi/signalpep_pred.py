#!/bin/python3

#TITLE: Detect signaling peptides using signalp

import os
import sys
import subprocess
import pandas as pd

#########################################
# FUNCTION: PARSE INPUT
#########################################
def table_to_fasta_sp(ampcombi_file):
    """
    This parses the ampcombi file and converts table to fasta
    """
    # mmseqsclust final table
    ampcombi = pd.read_csv(ampcombi_file, sep='\t')
    # add the corresponding seq_headers to ampcombi table for merging (DONT FORGET TO DROP IT AM ENDE)
    ampcombi['seq_headerss'] = ampcombi['sample_id'] + '!' + ampcombi['contig_id'] + '!' + ampcombi['CDS_id']
    selected_columns = ['seq_headerss', 'aa_sequence']
    ampcombi_df = ampcombi[selected_columns].copy()
    # table to fasta
    with open('representative_seq.fasta', 'w') as fasta_file:
        for _, row in ampcombi_df.iterrows():
            header = row['seq_headerss']
            sequence = row['aa_sequence']
            fasta_entry = f'>{header}\n{sequence}\n'
            fasta_file.write(fasta_entry)
    return ampcombi

#########################################
# FUNCTION: PREDICT SIGNAL PEPTIDE
#########################################
def signalp6(wdir, signalp6_model, ampcombi):
    """
    This detects the presence of signaling peptides if any.
    By running SignalP6slow-sequential
    """
    signalp_commands =f'signalp6 --fastafile representative_seq.fasta --model_dir {signalp6_model} --organism other --output_dir {wdir} --format all --mode slow-sequential'
    # execute the signalp6 command
    try:
        subprocess.run(signalp_commands, shell=True, text=True)
    except subprocess.CalledProcessError as e:
        print("Error running SignalP6 , please check the log file.", e)
            
    # remove unnecessary files from signalip dir:
    extensions = ['plot.txt', '.json', '.gff3', 'entries.fasta']
    file_list = os.listdir(wdir)
    for file_name in file_list:
        for extension in extensions:
            if file_name.endswith(extension):
                file_path = os.path.join(wdir, file_name)
                os.remove(file_path)
    
    # parse the results from signalp6
    # map the SP detected to ampcombi
    signalip = pd.read_csv(f'{wdir}/prediction_results.txt', sep='\t', header=[1])
    # rename header
    signalip = signalip.rename(columns={'# ID': 'seq_headerss'})
    # assign yes or no according to SP
    signalip['signal_peptide'] = signalip['CS Position'].str.contains('CS').map({True: 'yes', False: 'no'})
    selected_columns = ['seq_headerss', 'signal_peptide']
    signalip_df = signalip[selected_columns].copy()

    # grab the seqheaderss and index for alphafold later
    selected_columns = ['seq_headerss']
    pred_results = signalip_df[selected_columns].copy()
    pred_results['indices'] = range(len(pred_results))
    pred_results.to_csv(f'{wdir}/prediction_results_index.tsv',sep='\t', index=False)

    # add whether a signal peptide was found to ampcombi hits results
    rep_seqs_signal = pd.merge(ampcombi, signalip_df, on=['seq_headerss'], how='left')
    rep_seqs_signal = rep_seqs_signal.drop_duplicates()
    rep_seqs_signal.to_csv("Ampcombi_summary_cluster_SP.tsv",sep='\t', index=False)
    return rep_seqs_signal

#########################################
# FUNCTION: REMOVE CLUSTERS THAT HAVE NO SP
#########################################
def remove_clusters_no_sp(rep_seqs_signal):
    """
    This removes clusters that have no signaling peptide detected by signalp6
    """
    # table to dictionary (cluster_id and signal_peptide) by uniq cluster_id 
    selected_columns = ['cluster_id', 'signal_peptide']
    clusters = rep_seqs_signal[selected_columns].copy()

    # if cluster_id has no 'yes' then remove the dictionary with that cluster
    # iterate over unique values in the first column
    clusters_dict = {}
    listdict = []
    for index, value in enumerate(clusters.iloc[:, 0].unique(), start=0):
        clusters_dict = {}  # Initialize the dictionary for each iteration
        values_list = clusters.loc[clusters.iloc[:, 0] == value, clusters.columns[1]].tolist()
        clusters_dict['cluster_id'] = value
        clusters_dict['signal_peptide'] = values_list
        if 'yes' in clusters_dict['signal_peptide']:
            listdict.append(clusters_dict.copy())

    cluster_ids = [d['cluster_id'] for d in listdict]

    ampcombi_filtered = rep_seqs_signal[rep_seqs_signal['cluster_id'].isin(cluster_ids)]
    
    # remove duplicates 
    ampcombi_filtered = ampcombi_filtered.drop_duplicates(keep='first')
    ampcombi_filtered.to_csv("Ampcombi_summary_cluster_SP_onlyclusterswithSP.tsv",sep='\t', index=False)
