#!/bin/python3

# TITLE: Cluster the AMP hits using MMseqs2 after generating the complete summary

import tempfile
import shutil
import subprocess
import os
import pandas as pd

#########################################
# FUNCTION: PARSE INPUT FILE
#########################################
def parsing_input_for_cluster(merged_df):
    """
    This parses the ampcombi file and prepares it for downstream processing by converting sequences to fasta files. 
    """
    # make a column for aa_sequnces headers
    merged_df['seq_headers'] = merged_df['sample_id'] + '!' + merged_df['contig_id'] + '!' + merged_df['CDS_id']
    # convert table to sequnecs in temporary folder
    temp_folder = 'clusters'
    os.makedirs(temp_folder, exist_ok=True)
    #os.makedirs(temp_folder, exist_ok=True)
    
    fasta_entries = []
    for index, row in merged_df.iterrows():
        sequence_id = row['seq_headers']
        sequence = row['aa_sequence']
        fasta_entry = f'>{sequence_id}\n{sequence}\n'
        fasta_entries.append(fasta_entry)
    fasta_content = ''.join(fasta_entries)
    # write to file 
    output_file_path = os.path.join(temp_folder, 'merged_fastas.fasta')
    with open(output_file_path, 'w') as fasta_file:
        fasta_file.write(fasta_content)
    return merged_df

########################################
# FUNCTION: CLUSTER USING MMSEQS2 
#########################################
def mmseqs_cluster(cov_mod,cluster_mode,coverage,seq_id,sensitivity,threads):
    """
    This clusters the sequences using mmseqs2 (please refer to https://mmseqs.com/latest/userguide.pdf).
    - This mode of clustering : - combines prefiltering - alignment -clustering : so will take long to process. 
                                - can choose single step or cascascade clustering please refer to (https://github.com/soedinglab/MMseqs2/wiki#clustering-databases-using-mmseqs-cluster-or-mmseqs-linclust) for more info
                                - the clustering module can need large amounts of memory
    - The coverage is bidirectiional query and target, so should have a coverage of 100% : --cov-mode 0 -c 1
    - The sequence identity --min-seq-id 0.4
    - The sensitivity of alignment for clustering can be changed with -s 7.5 (very sensitive) and 4.0, which is slightly less sensitive but faster (default)
    - The cluster mode: --cluster-mode 0 (default). At each step, it forms a cluster containing the representative sequence with the most alignments above the special or default thresholds with other sequences of the database and these matched sequences. Then, the sequences contained in the cluster are removed and the next representative sequence is chosen.
    # NOTE: cluster-mode 1 uses CD-HIT algorithm which first sets the cluster according to length and then clusters
    # NOTE: reassign cluster: --cluster-reassign removes sequences from the cascaded cluster result that do not fulfill the cluster criteria and reassigns them (if possible) to a different cluster
    # NOTE: --cov-mode 3 : only sequences are clustered that have a sequence length overlap greater than X% of the query sequence
    """
    # define mmseqs db and cluster and align commands 
    os.makedirs('tmp', exist_ok=True)
    mmseqs_commands = [
        f'mmseqs createdb ./clusters/merged_fastas.fasta ./clusters/sequenceDB -v 1',
        # The combination chosen was tested on the simulated data and its values that make most sense 
        f'mmseqs cluster ./clusters/sequenceDB ./clusters/sequenceDBclu ./tmp --cov-mode {cov_mod} --cluster-mode {cluster_mode} -c {coverage} --min-seq-id {seq_id} -s {sensitivity} --threads {threads} -v 1',
        f'mmseqs createtsv ./clusters/sequenceDB ./clusters/sequenceDB ./clusters/sequenceDBclu ./clusters/sequenceDBclu.tsv --threads {threads} -v 1'
    ]
    # execute the MMseqs2 command
    for command in mmseqs_commands:
        try:
            # subprocess.run(mmseqs_create_db, shell=True, text=True)
            subprocess.run(command, shell=True, text=True)
        except subprocess.CalledProcessError as e:
            print("Error running MMseqs2 clustering, please check the log file.", e)
    
########################################
# FUNCTION: COMPILE THE CLUSTERS IN TABLE
#########################################
def compile_clusters(merged_df,retain_clusters_with, remove_singletons, min_cluster_members):
    # collect the clusters in a list and df
    cluster = "./clusters/sequenceDBclu.tsv"
    clusters = pd.read_csv(cluster, sep='\t', header=None)

    listdict = []

    # iterate over unique values in the first column
    for index, value in enumerate(clusters.iloc[:, 0].unique(), start=0):
        dictionary = {}
        # get the corresponding values from the second column as a list
        values_list = clusters.loc[clusters.iloc[:, 0] == value, 1].tolist()
        # remove the value if its = to the key (no clusters found)
        values_list = [v for v in values_list if v != value]
        # convert float to string in values_list
        values_list = [str(value) for value in values_list]
        # list to string
        values_string = '!!'.join(values_list)
        # add the list of values as an item for the key in the dictionary
        dictionary['seq_headers'] = value
        dictionary['cluster'] = value + '!!' + values_string + '*'
        # remove all keys that have no clusters ~singletons
        if '!!*' in dictionary['cluster']:
            continue
        else:
            # remove asterisk
            dictionary['cluster'] = dictionary['cluster'][:-1]
            dictionary['index'] = index
            # count the number of !! == values
            dictionary['total_cluster_members'] = int(dictionary['cluster'].count('!!') +1)
            # split the values in cluster
            dictionary['cluster'] = dictionary['cluster'].split('!!')
            # remove clusters that dont have 'specific value' :::: remove clusters that have only metaspades
            if retain_clusters_with is not None:
                if any(retain_clusters_with in value for value in dictionary['cluster']):
                    # remove clusters with sp. amount
                    if dictionary['total_cluster_members'] > min_cluster_members:
                        listdict.append(dictionary.copy())
            else:
                for value in dictionary['cluster']:
                    # remove clusters with  sp. amount
                    if dictionary['total_cluster_members'] > min_cluster_members:
                        listdict.append(dictionary.copy())
                
 
    # make a dataframe with sequence headers and index
    new_df = pd.DataFrame(merged_df['seq_headers'])
    new_df['cluster_id'] = None 

    # iterate over the df and dictionary and add the corresponding values from dict to df, comparison
    for index, row in new_df.iterrows():
        seq_header = str(row['seq_headers'])
        for dicto in listdict:
            for value in dicto['cluster']:
                if seq_header in str(value):
                    # adds the dict key value to df column
                    new_df.at[index, 'cluster_id'] = dicto['index']
                    # break teh inner loop to avoid unencessary iterations
                    break  
                    
    # write the clusters represnetative (e.g. after filtering out the modern clusters and the clusters that have < 3 members) in a file (dict to df)
    representative_seq = pd.DataFrame(listdict, columns=['seq_headers', 'index', 'total_cluster_members'])
    representative_seq.to_csv(f'Ampcombi_summary_cluster_representative_seq.tsv', sep='\t', index=False) # REMOBER TO ACTIVATE IN FUNCTION

    # clear the dictionary to clear memory allocated
    listdict.clear()                      
    dictionary.clear()

    # merge to ampcombi_summary and remove those not found in clusters (remove singletons)
    ampcombi_cluster = pd.merge(merged_df, new_df, on=['seq_headers'], how='left')

    # remove the new_df to release memory
    del new_df
    # remove temp dir for the clustering step
    shutil.rmtree('./clusters')
    shutil.rmtree('./tmp')
    
    # remove singletons (= None)
    if remove_singletons == True:
        ampcombi_cluster = ampcombi_cluster.dropna(subset=['cluster_id'])

    # remove cds_ID to grab only the gbk files corresponding to the dereplicated hits 
    ampcombi_cluster['seq_headers'] = ampcombi_cluster['seq_headers'].str.split('!', n=2).str[:2].str.join('!')
    ampcombi_cluster['seq_headers'] = ampcombi_cluster['seq_headers'].str.replace('!', '_')
    ampcombi_cluster.to_csv(f'Ampcombi_summary_cluster.tsv', sep='\t', index=False)
