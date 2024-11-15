#!/bin/python3

# TITLE: Download the ref database if input db empty AND and make database compatible for mmseqs

import pandas as pd
import requests
import os
import subprocess
import re
import pandas as pd
import tempfile
import shutil

from datetime import datetime
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

########################################
#  FUNCTION: DOWNLOAD DATABASES AND CLEAN DRAMP and APD
#########################################
def download_ref_db(db, database, threads):
    """
    Downloads a specified AMP (antimicrobial peptide) reference database based on the 
    provided database name and saves it to the specified directory. Supports downloading 
    databases from DRAMP, APD, and UniRef100.
    Parameters:
    ----------
    db : str
        The directory path where the downloaded database should be saved.
    database : str
        The name of the database to download. Must be one of 'DRAMP', 'APD', or 'UniRef100'.
    threads : int
        Number of threads to use when downloading the UniRef100 database with `mmseqs`.
    """
    # check which database was given 
    if database == 'DRAMP':
        #Download the file
        try:
            url = 'http://dramp.cpu-bioinfor.org/downloads/download.php?filename=download_data/DRAMP3.0_new/general_amps.txt'
            response = requests.get(url, allow_redirects=True)
            response.raise_for_status()  # Check for any download errors
            date = datetime.now().strftime("%Y_%m_%d")
            with open(db + '/' + f'general_amps_{date}.txt', 'wb') as file:
                file.write(response.content)
            print(f"File downloaded successfully and saved to {db}/general_amps_{date}.txt")
            #Create fasta version and clean it
            db_df = pd.read_csv(f'{db}/general_amps_{date}.txt', sep='\t')
            records = []
            valid_sequence_pattern = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")
            for index, row in db_df.iterrows():
                sequence = row['Sequence']
                if valid_sequence_pattern.match(sequence):
                    record = SeqRecord(Seq(sequence), id=str(row['DRAMP_ID']), description="")
                    records.append(record)
            output_file = f'{db}/general_amps_{date}.fasta'
            SeqIO.write(records, output_file, "fasta")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download DRAMP AMP general database file: {e}")
            return
    
    if database == 'APD':
        #Download the file 
        try:
            url = 'https://aps.unmc.edu/assets/sequences/APD_sequence_release_09142020.fasta'
            response = requests.get(url, allow_redirects=True, verify=False)  # Disable SSL verification due to site certificate issue
            response.raise_for_status()
            content = response.text
            print("APD AMP database downloaded successfully.")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download content: {e}")
            return
        #Save the content line-by-line exactly as is
        try:
            with open(db + '/' + 'APD_orig.fasta', 'w') as file:
                file.write(content)
            with open(f'{db}/APD.fasta', 'w') as output_handle:
                valid_sequence_pattern = re.compile("^[ACDEFGHIKLMNPQRSTVWY]+$")
                for record in SeqIO.parse(f'{db}/APD_orig.fasta', "fasta"):
                    sequence = str(record.seq)
                    if valid_sequence_pattern.match(sequence):
                        SeqIO.write(record, output_handle, "fasta")
            os.remove(db + '/' + 'APD_orig.fasta')
            print(f"APD AMP database saved successfully to {db}/APD.fasta")
            #fasta to table
            headers = []
            sequences = []
            seq_ids = []
            for i, record in enumerate(SeqIO.parse(f'{db}/APD.fasta', "fasta")):                 
                sequence_id = record.description.split('|')[0]
                headers.append(record.description)
                sequences.append(str(record.seq))
                seq_ids.append(sequence_id)
            db_df = pd.DataFrame({
                "APD_ID": seq_ids, 
                "APD_Description": headers,
                "APD_Sequence": sequences})
            db_df.to_csv(f'{db}/APD.txt', sep='\t', index=False, header=True)
            os.remove(db + '/' + 'APD.fasta')
            #table to fasta
            records = []
            for index, row in db_df.iterrows():
                sequence = row['APD_Sequence']
                record = SeqRecord(Seq(sequence), id=str(row['APD_ID']), description="")
                records.append(record)
            output_file = f'{db}/APD.fasta'
            SeqIO.write(records, output_file, "fasta")
        except Exception as e:
            print(f"Failed to save APD AMP database: {e}")
    
    if database == 'UniRef100':
        #Download the file
        try:
            os.makedirs(f'{db}/mmseqs2', exist_ok=True)
            command = f"mmseqs databases UniRef100 {db}/mmseqs2/ref_DB {db}/mmseqs2/tmp --remove-tmp-files true --threads {threads} -v 0"
            subprocess.run(command, shell=True, check=True)
            print(f"UniRef100 protein database downloaded successfully and saved to {db}/mmseqs2/UniRef100")
        except subprocess.CalledProcessError as e:
            print(f"Failed to download UniRef100 protein database: {e}")
        
########################################
#  FUNCTION: CREATE MMseqs2 COMPATIBLE DATABASE FORMATS
#########################################
def create_mmseqs_ref_db(db):
    """
    Creates an MMseqs2 reference database for a given database in the specified directory,
    either supplied by user or given.
    Parameters:
    ----------
    db : str
        The directory path where the database files should be.
    """
    try: 
        
        db_path = os.path.join(os.path.abspath(db))
        db_path_mmseqs = os.path.join(db_path, 'mmseqs2')
        
        if os.path.exists(db_path_mmseqs):
            print(f"MMseqs2 directory already exists in {db}.")
            return
        else:
            print(f"Creating MMseqs2 directory for fasta file in {db}...")
            os.makedirs(db_path_mmseqs)
            fasta_file = os.path.join(db_path, '*.fasta')
            output_db = os.path.join(db_path_mmseqs, f"ref_DB")
            command = f"mmseqs createdb {fasta_file} {output_db} -v 0"
            subprocess.run(command, shell=True, check=True)  
    
    except Exception as e:
        print(f"Failed to create MMseqs2 database with fasta files in {db}: {e}")    

########################################
#  FUNCTION: MMseqs2 ALIGNMENT
#########################################
def mmseqs_alignment(db, amp_faa_paths, amp_matches, threads):
    """
    Performs an MMseqs2 alignment for each FASTA file in the provided list of AMP (antimicrobial peptide) 
    file paths against a reference database located in the specified directory. 
    This function sets search parameters based on the specific database name within the directory 
    and creates temporary directories for intermediate files.

    Parameters:
    -----------
    db : str
        The path to the directory containing the reference database for alignment.
        
    amp_faa_paths : list of str
        A list of file paths for AMP FASTA files to be aligned against the reference database.
        
    amp_matches : str
        The output path for storing alignment results from `mmseqs convertalis` in the specified format.
        
    threads : int
        The number of CPU threads to use for MMseqs2 processes, enabling parallelization for faster processing.
    """
    for faa_path in amp_faa_paths:
        cwd = os.getcwd()
        db_path = os.path.join(os.path.abspath(db))
        db_path_mmseqs = os.path.join(db_path, 'mmseqs2')

        if any(database_item in os.path.basename(db_path) for database_item in ['DRAMP', 'APD', 'UniRef100']):
            if 'DRAMP' in os.path.basename(db_path) or 'APD' in os.path.basename(db_path):
                min_seq_id = 0.1
                s_value = 8.0
                max_seqs = 1000000
            elif 'UniRef100' in os.path.basename(db_path):
                min_seq_id = 0.6
                s_value = 5.0
                max_seqs = 500
        else: 
            min_seq_id = 0.1
            s_value = 8.0
            max_seqs = 1000000

        temp = tempfile.mkdtemp()
        command_create_db = f"mmseqs createdb {faa_path} {temp}/shortseq_DB -v 0"
        command_search = f"mmseqs search {temp}/shortseq_DB {db_path_mmseqs}/ref_DB {temp}/results_DB {temp}/tmp --threads {threads} --cov-mode 1 -c 0.6 --min-seq-id {min_seq_id} -s {s_value} -a --max-seqs {max_seqs} --strand 2 --max-accept 1 -k 5 --sort-results 1 -v 0"
        command_convertalis = f"mmseqs convertalis {temp}/shortseq_DB {db_path_mmseqs}/ref_DB {temp}/results_DB {amp_matches} --format-output query,target,evalue,pident,nident,tlen,tstart,tend,taln,theader,alnlen,qcov,tcov -v 0"
        try:
            subprocess.run(command_create_db, shell=True, check=True, text=True)
            subprocess.run(command_search, shell=True, check=True, text=True)
            subprocess.run(command_convertalis, shell=True, check=True, text=True)
            return amp_matches
        except subprocess.CalledProcessError as e:
            print("Failed to run mmseqs search:", e) 
        finally: 
            shutil.rmtree(temp)
        
########################################
#  FUNCTION: MMseqs2 Merge to summary
#########################################
def mmseq_alignment_merge(db, amp_faa_paths, amp_matches, threads, dbevalue, summary_df_filtered, sample):
    """
    Merges MMseqs alignment results with reference database information based on a specific e-value cutoff.
    Parameters:
    -----------
    db : str
        Path to the directory containing the reference database files.
    amp_faa_paths : list of str
        List of AMP FASTA file paths to align.
    amp_matches : str
        Path where the alignment results are saved.
    threads : int
        Number of threads for parallel processing.
    dbevalue : float
        E-value cutoff for filtering alignment results.
    """
    # check first if faa is empty or not 
    for faa_path in amp_faa_paths:
        if (not os.path.exists(faa_path)) or (os.stat(faa_path).st_size == 0):
            print(f'Alignment to ref. database is skipped for {sample}, because no AMPs were found after filtering with your thresholds.')
            return summary_df_filtered
        else:
            mmseqs_search_path = mmseqs_alignment(db, amp_faa_paths, amp_matches, threads)
            if (os.stat(mmseqs_search_path).st_size == 0):
                print(f'No AMP hits for {sample} were successfully aligned to the reference database. Consider changing the database.')
                return summary_df_filtered
            else:
                try:
                    dd_align = pd.read_csv(mmseqs_search_path, sep='\t', header= None)
                    #Rename columns according to dict
                    dd_align.columns = ["query", 
                                                  "target",
                                                  "evalue",
                                                  "pident",
                                                  "nident",
                                                  "tlen",
                                                  "tstart",
                                                  "tend",
                                                  "taln",
                                                  "theader",
                                                  "alnlen",
                                                  "qcov",
                                                  "tcov"]
                    #Reassign the contents of the file
                    dd_align.to_csv(amp_matches, sep='\t', index=False)
                    #Make sure the evalues are float type 
                    dd_align['evalue'] = dd_align['evalue'].astype(float)
                    #Remove only the classification below evalue
                    dd_align = dd_align[dd_align['evalue'] <= float(dbevalue)]

                    #Merge to the ref db to get more info according to the db
                    db_path = os.path.join(os.path.abspath(db))

                    if any(database_item in os.path.basename(db_path) for database_item in ['DRAMP', 'APD', 'UniRef100']):
                        if 'DRAMP' in os.path.basename(db_path):
                            for file in os.listdir(db_path):
                                if file.endswith('.txt'):
                                    dramp_file = os.path.join(db_path, file)
                                    dramp_df = pd.read_csv(dramp_file, sep='\t', usecols=['DRAMP_ID','Sequence','Name','Family','Source','Target_Organism','Gene','Swiss_Prot_Entry','PDB_ID'])
                                    merged = pd.merge(dd_align, dramp_df, left_on='target',right_on='DRAMP_ID', how='left')
                                    return merged
                        elif 'APD' in os.path.basename(db_path):
                            for file in os.listdir(db_path):
                                if file.endswith('.txt'):
                                    apd_file = os.path.join(db_path, file)
                                    apd_df = pd.read_csv(apd_file, sep='\t')
                                    merged = pd.merge(dd_align, apd_df, left_on='target',right_on='APD_ID', how='left')
                                    return merged
                        elif 'UniRef100' in os.path.basename(db_path):
                            # grab everything between the first space and the n=
                            dd_align['header'] = dd_align['theader'].str.extract(r'\s(.*?)\s*n=')
                            return dd_align
                except Exception as e:
                    print(f"Error reading mmseqs search/alignment results: {e}")
