#!/bin/python3

# TITLE: Reformat the AMP output tables

import pandas as pd
from Bio import SeqIO
from check_input import check_dfshape

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
    return ampir_df[['contig_id', 'prob_ampir']]

#########################################
    #  AMP_amplify
#########################################
def amplify(path, p):
    amplify_dict = {'Sequence_ID':'contig_id', 'Sequence':'seq_aa', 'Length':'length', 'Charge':'charge', 'Probability_score':'prob_amplify', 'AMPlify_log_scaled_score':'log_score', 'Prediction':'prediction'}
    amplify_df = pd.read_csv(path, sep='\t').rename(columns=amplify_dict).dropna()
    # apply probability cutoff
    amplify_df = amplify_df[(amplify_df['prob_amplify']>=p)]
    return amplify_df[['contig_id', 'prob_amplify']]

#########################################
    #  AMP_ensembleamppred
#########################################
def amppred(path, p):
    trim_text(path, 'Sequence')
    amppred_dict = {4:'index', 14:'prob_amppred'} #{'level_0':1, 'level_1':2, 'level_2':'index', 'level_3':3, 'level_4':4, 'level_5':5, '############':6, 'Prediction':7, 'results':8, 'by':9, 'EnsembleAMPPred':10, '#############':'prob_amppred'}
    amppred_df = pd.read_csv(path, sep=' ', header=None).rename(columns=amppred_dict)
    amppred_df = amppred_df[(amppred_df['prob_amppred']>=p)]
    return amppred_df[['index', 'prob_amppred']]

#########################################
    #  AMP_macrel
#########################################
def macrel(path, p):
    macrel_dict = {'Access':'contig_id', 'Sequence':'seq_aa', 'AMP_family':'amp_family', 'AMP_probability':'prob_macrel', 'Hemolytic':'hemolytic', 'Hemolytic_probability':'prob_hemo'}
    #set header to second row to skip first line starting with #
    macrel_df = pd.read_csv(path, sep='\t', header=[1]).rename(columns=macrel_dict)
    # apply probability cutoff
    macrel_df = macrel_df[(macrel_df['prob_macrel']>=p)]
    return macrel_df[['contig_id', 'prob_macrel']]

#########################################
    #  AMP_neubi
#########################################
def neubi(path, p):
    neubi_seq = SeqIO.parse(open(path), 'fasta')
    #initiate the dataframe containing contig ids, aa-sequences and probability in three columns
    neubi_df = pd.DataFrame(columns=['contig_id', 'aa_sequence', 'prob_neubi'])
    #append contig information to df (p is last value in header following '|' symbol)
    for contig in neubi_seq:
        contig_id, sequence, description = contig.id, str(contig.seq), float(contig.description.split("|",1)[1])
        neubi_df = neubi_df.append({'contig_id':contig_id, 'aa_sequence':sequence, 'prob_neubi':description}, ignore_index=True)
    neubi_df = neubi_df[(neubi_df['prob_neubi']>=p)]
    return neubi_df[['contig_id', 'prob_neubi']]

#########################################
    #  AMP_hmmsearch
#########################################
def hmmsearch(path):
    hmmer_dict = {'level_0':'evalue_hmmer', 'level_1':'score_hmmer', 'level_2':'bias', 'level_3':'eval_domain', 'level_4':'score_domain', 'level_5':'bias_domain', 'level_6':'exp_dom', '-------':'N_dom', '------':'contig_id'}
    hmmer_df = pd.read_table(path, delim_whitespace=True, header=[15]).reset_index().rename(columns=hmmer_dict)
    hmmer_df = hmmer_df.drop(hmmer_df.iloc[:,9:17], axis=1) #drop unnecessary columns
    for index, row in hmmer_df.iterrows():
        if (row.str.contains('Domain').any()):              #identify index of first row with 'Domain'
            i = index
            break
    hmmer_df = hmmer_df[hmmer_df.index<i]                   #only keep rows previous to index i
    return hmmer_df[['contig_id', 'evalue_hmmer']]

#########################################
# FUNCTION: READ DFs PER SAMPLE 
#########################################
# For one sample: parse filepaths and read files to dataframes, create list of dataframes
def read_path(df_list, file_list, p, dict, faa_path, samplename):
    for path in file_list:
        if(path.endswith(dict['ampir'])):
            print('found ampir file')
            df_list.append(ampir(path, p))
        elif(path.endswith(dict['amplify'])):
            print('found amplify file')
            df_list.append(amplify(path, p))
        elif(path.endswith(dict['macrel'])):
            print('found macrel file')
            df_list.append(macrel(path, p))
        elif(path.endswith(dict['neubi'])):
            print('found neubi file')
            df_list.append(neubi(path, p))
        elif(path.endswith(dict['hmmer_hmmsearch'])):
            print('found hmmersearch file')
            df_list.append(hmmsearch(path))
        elif(path.endswith(dict['ensembleamppred'])):
            print('found ensemblamppred file')
            faa_filepath = faa_path+samplename+'.faa'
            faa_df = faa2table(faa_filepath)
            amppred_df = amppred(path, p)
            if(check_dfshape(amppred_df, faa_df)):
                # add contig_ids via index numbers, because ensembleamppred only gives numbered sequences without ids, in the order of sequences in faa
                amppred_df = pd.merge(amppred_df, faa_df.reset_index(), on='index')
                amppred_df.drop(['index', 'aa_sequence'], axis=1)
                df_list.append(amppred_df)
        else:
            print(f'No AMP-output-files could be found with the given path ({path}). \n Please check your file paths and file endings or use the <--path-list> command')
            break

#########################################
# FUNCTION: MERGE DATAFRAMES
#########################################
# merge dataframes from list to summary output per sample
def summary(df_list, samplename, faa_path, outdir):
    #initiate merge_df
    merge_df = pd.DataFrame(columns=['contig_id'])
    #merge all dfs in the df-list on contig_id
    for df in df_list:
        merge_df = pd.merge(merge_df, pd.DataFrame(df) , how='outer', on='contig_id')
    #replace all NAs (where a tool did not identify the contig as AMP) with 0
    merge_df = merge_df.fillna(0)
    #add amino-acid sequences
    faa_df = faa2table(faa_path+samplename+'.faa')
    merge_df = merge_df.merge(faa_df, how='inner', on='contig_id')
    # sort by sum of p-values over rows
    merge_df = merge_df.set_index('contig_id')
    merge_df['p_sum']= merge_df.sum(axis=1)#.sort_values(ascending=False)
    merge_df = merge_df.sort_values('p_sum', ascending=False).drop('p_sum', axis=1).reset_index()
    # write summary to outdir
    #merge_df.to_csv(outdir+'/'+samplename+'_AMPsummary.csv', sep=',')
    return merge_df

#########################################
# FUNCTION: READ FAA TO TABLE
#########################################
# transform faa to dataframe with two columns
def faa2table(faa_path):
    #read the amino-acid fasta with SeqIO
    faa_seq = SeqIO.parse(open(faa_path), 'fasta')
    #initiate the dataframe containing contig ids and aa-sequences in two columns
    fasta_df = pd.DataFrame(columns=['contig_id', 'aa_sequence'])
    #append contig information to df
    for contig in faa_seq:
        contig_id, sequence = contig.id, str(contig.seq)
        fasta_df = fasta_df.append({'contig_id':contig_id, 'aa_sequence':sequence}, ignore_index=True)
    return fasta_df