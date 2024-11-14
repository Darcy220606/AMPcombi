#!/bin/python3

# TITLE: Reformat the AMP output tables

import pandas as pd
from Bio import SeqIO
import os

import tempfile
import shutil
import subprocess

#########################################
# FUNCTION: KEEP ONLY LINES WITH KEYWORD
#########################################
def trim_text(filepath, key):
    lines = []
    # read file
    with open(filepath, 'r') as fp:
        # read an store all lines into list
        lines = fp.readlines()

    # Write file
    with open(filepath, 'w') as fp:
        # iterate each line
        for line in lines:
            if key in line:
                fp.write(line)

#########################################
# FUNCTION: KEEP ONLY LINES WITH KEYWORD
#########################################
def check_dfshape(df1, df2):
    if (df1.shape[0] != df2.shape[0]):
        print(f'ERROR: different row number in tool output and faa file. Ensembleamppred output could not be included in the summary')
        return False
    else: 
        return True


#########################################
# FUNCTIONS: READ TOOLS' OUTPUT TO DFs
#########################################

#########################################
    #  AMP_ampir
#########################################
def ampir(path, p): 
    # Dictionary to rename columns
    ampir_dict = {'seq_name':'contig_id', 'seq_aa':'seq_aa', 'prob_AMP':'prob_ampir'}
    # read file as df and rename columns
    ampir_df = pd.read_csv(path, sep='\t').rename(columns=ampir_dict) 
    # cut contig_id to remove extra information added by tool
    ampir_df['contig_id']=ampir_df['contig_id'].apply(lambda x: x.split()[0])
    # apply probability cutoff
    ampir_df = ampir_df[(ampir_df['prob_ampir']>=p)]
    ampir_df['prob_ampir'] = ampir_df['prob_ampir'].round(3)
    return ampir_df[['contig_id', 'prob_ampir']]

#########################################
    #  AMP_amplify
#########################################
def amplify(path, p):
    amplify_dict = {'Sequence_ID':'contig_id', 'Sequence':'seq_aa', 'Length':'length', 'Charge':'charge', 'Probability_score':'prob_amplify', 'AMPlify_log_scaled_score':'log_score', 'Prediction':'prediction'}
    amplify_df = pd.read_csv(path, sep='\t').rename(columns=amplify_dict)
    # apply probability cutoff
    amplify_df = amplify_df[(amplify_df['prob_amplify']>=p)]
    amplify_df['prob_amplify'] = amplify_df['prob_amplify'].round(3)
    return amplify_df[['contig_id', 'prob_amplify']]

#########################################
    #  AMP_macrel
#########################################
def macrel(path, p):
    macrel_dict = {'Access':'contig_id', 'Sequence':'seq_aa', 'AMP_family':'amp_family', 'AMP_probability':'prob_macrel', 'Hemolytic':'hemolytic', 'Hemolytic_probability':'prob_hemo'}
    #set header to second row to skip first line starting with #
    macrel_df = pd.read_csv(path, sep='\t', header=[1]).rename(columns=macrel_dict)
    # apply probability cutoff
    macrel_df = macrel_df[(macrel_df['prob_macrel']>=p)]
    macrel_df['prob_macrel'] = macrel_df['prob_macrel'].round(3)
    return macrel_df[['contig_id', 'prob_macrel']]

#########################################
    #  AMP_neubi
#########################################
def neubi(path, p):
    neubi_seq = SeqIO.parse(open(path), 'fasta')
    # initiate the dataframe containing contig ids, aa-sequences and probability in three columns
    neubi_df = pd.DataFrame(columns=['contig_id', 'aa_sequence', 'prob_neubi'])
    # append contig information to df (p is last value in header following '|' symbol)
    for contig in neubi_seq:
        contig_id, sequence, description = contig.id, str(contig.seq), float(contig.description.split("|",1)[1])
        neubi_df = neubi_df.append({'contig_id':contig_id, 'aa_sequence':sequence, 'prob_neubi':description}, ignore_index=True)
    neubi_df = neubi_df[(neubi_df['prob_neubi']>=p)]
    neubi_df['prob_neubi'] = neubi_df['prob_neubi'].round(3)
    return neubi_df[['contig_id', 'prob_neubi']]

#########################################
    #  AMP_ampgram
#########################################
def ampgram(path, p): 
    # dictionary to rename columns
    ampgram_dict = {'single_prot_pred':'prob_ampgram'}
    # read file as df and rename columns
    ampgram_df = pd.read_csv(path, sep='\t').rename(columns=ampgram_dict) 
    # apply probability cutoff
    ampgram_df = ampgram_df[(ampgram_df['prob_ampgram']>=p)]
    ampgram_df['prob_ampgram'] = ampgram_df['prob_ampgram'].round(3)
    return ampgram_df[['contig_id', 'prob_ampgram']]

##########################################
    #  AMP_hmmsearch (single and multi HMM models)
##########################################
def hmmsearch(path, hmmevalue):
    # create temp folder to convert hmm to csv
    temp = tempfile.mkdtemp()
    # define the command to execute
    command = f'hmm_to_csv_input_file.py -i {path} -o ./temp/hmm_stats.csv'
    subprocess.run(command, text=True, shell=True)
    # open the csv file and extract what you need from there
    hmm_df = pd.read_csv('./temp/hmm_stats.csv', sep=',', header = None, skiprows=1,usecols=[1, 3, 11])
    hmm_df.columns = ['HMM_model','evalue_hmmer','contig_id']
    hmm_df['evalue_hmmer'] = hmm_df['evalue_hmmer']
    hmm_df['evalue_hmmer'] = pd.to_numeric(hmm_df['evalue_hmmer'], errors='coerce')
    # remove any hits below evalue specified
    if hmmevalue is not None:
        # make sure the evalues are float type 
        # remove any hits below evalue
        hmm_df = hmm_df[hmm_df['evalue_hmmer'] < float(hmmevalue)]
    return hmm_df[['HMM_model', 'evalue_hmmer', 'contig_id']] 

#########################################
    #  AMP_transformer
#########################################
def amptransformer(path, p):
    # Dictionary to rename columns
    amptransformer_dict = {'peptides':'contig_id','sequence':'seq_aa','Antimicrobial_Peptide_Prediction':'prob_amptransformer'}
    # read file as df and rename columns
    amptransformer_df = pd.read_csv(path, sep='\t').rename(columns=amptransformer_dict)
    # apply probability cutoff
    amptransformer_df = amptransformer_df[(amptransformer_df['prob_amptransformer']>=p)]
    amptransformer_df['prob_amptransformer'] = amptransformer_df['prob_amptransformer'].round(3)
    return amptransformer_df[['contig_id', 'prob_amptransformer']]

#########################################
    #  AMP_ensembleamppred
#########################################
def amppred(path, p):
    trim_text(path, 'Sequence')
    amppred_dict = {4:'index', 14:'prob_amppred'} #{'level_0':1, 'level_1':2, 'level_2':'index', 'level_3':3, 'level_4':4, 'level_5':5, '############':6, 'Prediction':7, 'results':8, 'by':9, 'EnsembleAMPPred>
    amppred_df = pd.read_csv(path, sep=' ', header=None).rename(columns=amppred_dict)
    amppred_df = amppred_df[(amppred_df['prob_amppred']>=p)]
    amppred_df['prob_amppred'] = amppred_df['prob_amppred'].round(3)
    return amppred_df[['index', 'prob_amppred']]

#########################################
# FUNCTION: READ DFs PER SAMPLE 
#########################################
# For one sample: parse filepaths and read files to dataframes, create list of dataframes
def read_path(df_list, file_list, p, hmmevalue, dict, faa_path, samplename):
    for path in file_list:
        try:
            if 'ampir' in path and path.endswith(dict.get('ampir')):
                print('found ampir file')
                df_list.append(ampir(path, p))
            if 'amplify' in path and path.endswith(dict.get('amplify')):
                print('found amplify file')
                df_list.append(amplify(path, p))
            if 'macrel' in path and path.endswith(dict.get('macrel')):
                print('found macrel file')
                df_list.append(macrel(path, p))
            if 'ampgram' in path and path.endswith(dict.get('ampgram')):
                print('found ampgram file')
                df_list.append(ampgram(path, p))
            if 'amptransformer' in path and path.endswith(dict.get('amptransformer')):
                print('found amptransformer file')
                df_list.append(amptransformer(path, p))
            if 'neubi' in path and path.endswith(dict.get('neubi')):
                print('found neubi file')
                df_list.append(neubi(path, p))
            if 'hmmer_hmmsearch' in path and path.endswith(dict.get('hmmer_hmmsearch')):
                print('found hmmersearch file')
                df_list.append(hmmsearch(path, hmmevalue))
            if 'ensembleamppred' in path and path.endswith(dict.get('ensembleamppred')):
                print('found ensemblamppred file')
                faa_filepath = faa_path+'/'+samplename+'.faa'
                faa_df = faa2table(faa_filepath)
                amppred_df = amppred(path, p)
                if(check_dfshape(amppred_df, faa_df)):
                    # add contig_ids via index numbers, because ensembleamppred only gives numbered sequences without ids, in the order of sequences in faa
                    amppred_df = pd.merge(amppred_df, faa_df.reset_index(), on='index')
                    amppred_df.drop(['index', 'aa_sequence'], axis=1)
                    df_list.append(amppred_df)
        except:
            print(f'No AMP-output-files could be found with the given path ({path}). \n Please check your file paths and file endings or use the <--path-list> command')
            #break
            pass

#########################################
# FUNCTION: MERGE DATAFRAMES
#########################################
# merge dataframes from list to summary output per sample
def summary(df_list, samplename, faa_path, aa_len):
    # initiate merge_df
    merge_df = pd.DataFrame(columns=['contig_id'])
    # merge all dfs in the df-list on contig_id
    for df in df_list:
        merge_df = pd.merge(merge_df, pd.DataFrame(df) , how='outer', on='contig_id')
    # replace all NAs (where a tool did not identify the contig as AMP) with 0
    merge_df = merge_df.fillna(0)
    # add amino-acid sequences
    faa_df = faa2table(faa_path)     
    merge_df = merge_df.merge(faa_df, how='inner', on='contig_id')
    # remove hits that have a hit lengths <100aa: TODO!!!
    # count the number of letters across the aa_sequence
    merge_df['aa_lengths'] =  merge_df['aa_sequence'].apply(lambda x: len(x))
    # retain hits below the aa lengths
    merge_df = merge_df[merge_df['aa_lengths'] <= aa_len]
    # sort by sum of p-values over rows
    merge_df = merge_df.set_index('contig_id')
    merge_df['p_sum']= merge_df.sum(axis=1)#.sort_values(ascending=False)
    merge_df = merge_df.sort_values('p_sum', ascending=False).drop(['p_sum', 'aa_lengths'], axis=1).reset_index()
    # cleanup remove temp dir
    # shutil.rmtree('./temp')
    return merge_df

#########################################
# FUNCTION: READ FAA TO TABLE
#########################################
# transform faa to dataframe with two columns
def faa2table(faa_path):
    # read the amino-acid fasta with SeqIO
    faa_seq = SeqIO.parse(open(faa_path), 'fasta')
    # initiate the dataframe containing contig ids and aa-sequences in two columns
    fasta_df = pd.DataFrame(columns=['contig_id', 'aa_sequence'])
    # append contig information to df
    for contig in faa_seq:
        contig_id, sequence = contig.id, str(contig.seq)
        fasta_df = fasta_df.append({'contig_id':contig_id, 'aa_sequence':sequence}, ignore_index=True)
    return fasta_df
