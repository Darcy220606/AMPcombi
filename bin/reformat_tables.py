#!/bin/python3

# TITLE: Reformat the AMP output tables

import os
import pandas as pd

path ='/media/aibrahim/INTENSO/Rosa-table-parsing/AMP_pipeline/amp/'

# Main input format should be in a /maindir and /subdirectories named by tools and within that /samplesubdir
# E.g. /maindir/toolsubdir/samplesubdir/*.tsv

# Done with the funcscan result in mind 

#########################################
#  AMP_ampir
#########################################
#def ampir(path):    
#   """
#   This is to reformat the ampir output result table .tsv. 
#   Input includes only the path to the directory.
#   """
filelist = []
for dirpath, subdirs, files in os.walk(path):
    for file in files:
        if(file.endswith(".ampir.tsv")):
            filelist.append(os.path.join(dirpath, file)) #grab the full directory of the files
            #print(filelist)
            for file in filelist:
                #print(file)
                fields = ['seq_name', 'prob_AMP']
                df = pd.read_csv(file, sep='\t', usecols=fields) #Retain the 1st and 3rd column by name
                df['seq_name']=df['seq_name'].str.split(" ", expand=True)[0] #Split strings in the first column and remove everythng after the first space
                df.rename({'seq_name': 'contig_id','prob_AMP': 'prob'}, axis=1, inplace=True) #Replace the column headers
                #print(df)
        else:
            None
                
#########################################
#  AMP_macrel
#########################################
#def macrel(path):    
#   """
#   This is to reformat the macrel output result table .tsv. 
#   Input includes only the path to the directory.
#   """
filelist = []
for dirpath, subdirs, files in os.walk(path):
    for file in files:
        if(file.endswith(".macrel.tsv")): ## CHECK WITH LOUISA IF SHE ALSO GETS A TEMP FILE AFTER EXTRACTION and should we rename the files in FUNCSCAN ???????
            filelist.append(os.path.join(dirpath, file)) #grab the full directory of the files
            #print(filelist)
            for file in filelist:
                #print(file)
                fields = ['Access', 'AMP_probability']
                df = pd.read_csv(file, sep='\t',skiprows=[0], usecols=fields) #Remove the first commented row and retain the the seqID and prob columns
                df.rename({'Access': 'contig_id','AMP_probability': 'prob'}, axis=1, inplace=True) #Replace the column headers
                #print(df)
        else:
            None

#########################################
#  AMP_amplify
#########################################
#def amplify(path):    
#   """
#   This is to reformat the amplify output result table .tsv. 
#   Input includes only the path to the directory.
#   """
filelist = []
for dirpath, subdirs, files in os.walk(path):
    for file in files:
        if(file.endswith(".amplify.tsv")): ## CHECK WITH LOUISA IF SHE ALSO GETS A TEMP FILE AFTER EXTRACTION and should we rename the files in FUNCSCAN ???????
            filelist.append(os.path.join(dirpath, file)) #grab the full directory of the files
            #print(filelist)
            for file in filelist:
                #print(file)
                fields = ['Sequence_ID', 'Probability_score']
                df = pd.read_csv(file, sep='\t', usecols=fields) #Retain the the seqID and prob columns
                df['Probability_score'] = df['Probability_score'].fillna(0) #Replace the NaNs with zeros
                df.rename({'Sequence_ID': 'contig_id','Probability_score': 'prob'}, axis=1, inplace=True) #Replace the column headers
                print(df)
        else:None

#########################################
#  AMP_hmmsearch
#########################################
#def hmmsearch(path):    
#   """
#   This is to reformat the hmmsearch output result table .tsv. 
#   Input includes only the path to the directory.
#   """
filelist = []
for dirpath, subdirs, files in os.walk(path):
    for file in files:
        if(file.endswith(".hmmsearch.txt")): ## CHECK WITH LOUISA IF SHE ALSO GETS A TEMP FILE AFTER EXTRACTION and should we rename the files in FUNCSCAN ???????
            filelist.append(os.path.join(dirpath, file)) #grab the full directory of the files
            #print(filelist)
            for file in filelist:
                #print(file)
                #dictionary to rename columns
                hmmer_dict = {'level_0':'evalue_hmmer', 'level_1':'score_hmmer', 'level_2':'bias', 'level_3':'eval_domain', 'level_4':'score_domain', 'level_5':'bias_domain', 'level_6':'exp_dom', '-------':'N_dom', '------':'contig_id'}
                #read the txt-file with read_table, rename the first columns, delete unnecessary columns and rows
                df = pd.read_table(file, delim_whitespace=True, header=[15]).reset_index().rename(columns=hmmer_dict).drop(df.iloc[:,9:17], axis=1).dropna()
                print(df)
        else:None



#ampir()
#macrel()
#amplify()
#hmmsearch()