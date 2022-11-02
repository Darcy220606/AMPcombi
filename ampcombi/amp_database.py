#!/bin/python3

# TITLE: Download the DRAMP database if input db empty AND and make database compatible for diamond

import pandas as pd
import requests
import os
from datetime import datetime
import subprocess
from Bio import SeqIO
import tempfile
import shutil

########################################
#  FUNCTION: DOWNLOAD DRAMP DATABASE AND CLEAN IT
#########################################
def download_DRAMP(db):
    ##Download the (table) file and store it in a results directory 
    url = 'http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.xlsx'
    r = requests.get(url, allow_redirects=True)
    with open(db +'/'+ 'general_amps.xlsx', 'wb') as f:
        f.write(r.content)
    ##Convert excel to tab sep file and write it to a file in the DRAMP_db directly with the date its downloaded
    date = datetime.now().strftime("%Y_%m_%d")
    ref_amps=pd.read_excel (db +'/'+ r'general_amps.xlsx')
    ref_amps.to_csv (db +'/' + f'general_amps_{date}.tsv', index = None, header=True,sep='\t')
    ##Download the (fasta) file and store it in a results directory 
    urlfasta = 'http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.fasta'
    z = requests.get(urlfasta)
    fasta_path = os.path.join(db + '/' + f'general_amps_{date}.fasta')
    with open(fasta_path, 'wb') as f:     
        f.write(z.content)
    ##Cleaning step to remove ambigous aminoacids from sequences in the database (e.g. zeros and brackets)
    new_fasta =  db + '/' + f'general_amps_{date}_clean.fasta'      
    seq_record = SeqIO.parse(open(fasta_path), "fasta")
    with open(new_fasta, 'w') as f:
        for record in seq_record:
            id, sequence = record.id, str(record.seq)
            letters = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W',',Y']
            new = ''.join(i for i in sequence if i in letters)
            f.write('>' + id + '\n' + new + '\n')
    return os.remove(fasta_path), os.remove(db +'/'+ r'general_amps.xlsx')

########################################
#  FUNCTION: CREATE DIAMOND COMPATIBLE DATBASE FORMATS
#########################################
def create_diamond_ref_db(db,threads):
    cwd = os.getcwd()
    for file in os.listdir(db):
        if file.endswith('.fasta'):
            path = os.path.join(os.path.abspath(db) + '/' + file)
            os.chdir(db)
            #process = subprocess.Popen([f'{scripts_path}/diamond_makedb.sh', path])
            subprocess.run('diamond_makedb.sh', text=True, input=f'{path}\n{threads}')
            os.chdir(cwd)
            print
            return path
        
########################################
#  FUNCTION: DIAMOND ALIGNMENT
#########################################
def diamond_alignment(db, amp_faa_paths, amp_matches,threads):
    #create temp folder and delete at the end
    cwd = os.getcwd()
    for path in amp_faa_paths:
        # align the query with the database
        temp = tempfile.mkdtemp()
        subprocess.run('diamond_alignment.sh', text=True, input=f'{path}\n{temp}\n{db}\n{threads}')
        shutil.move(temp+'/diamond_matches.tsv', amp_matches)
        shutil.rmtree(temp)
        # mege the diamond_alignment with the ref_db table
        dd_align = pd.read_csv(amp_matches, delimiter='\t')
        dd_align = dd_align[['target_id','contig_id','pident','evalue']]
        for file in os.listdir(db):
            if file.endswith('.tsv'):
                path_2 = os.path.join(os.path.abspath(db) + '/' + file)
                ref_db = pd.read_csv(path_2, delimiter='\t')
                ref_db.columns.values[0] = "target_id"
                merged = pd.merge(dd_align, ref_db, on='target_id',how='inner')
                return merged