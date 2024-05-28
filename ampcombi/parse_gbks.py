#!/usr/bin/env python3

import pandas as pd
import os
import re

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from collections import defaultdict

########################################
#  FUNCTION: Parse the AMPcombi file 
#########################################
def ampcombi_input(sample_ampcombi_file):
    """
    This parses the ampcombi summary files (- ampcombi_file) and 
    creates the main df and the filtered dict. to be used later.
    """
    # ampcombi_summary file input
    ampcombi_main = sample_ampcombi_file
    grab_col = ['name', 'contig_id'] # contig_id is CDS_id in this case
    ampcombi = ampcombi_main[grab_col]
    # create a dictionary wih the table contents
    ampcombi_dict = ampcombi.to_dict(orient='records')
    combined_data = defaultdict(list)
    # combine the ampcombi dict according to the key (reduces the size by reducing duplicates)
    for item in ampcombi_dict:
        name = item['name']
        combined_data[name].append(item['contig_id'])
    # rename the dict items
    ampcombi_dict_mod = [{'name': key, 'contig_id': items} for key, items in combined_data.items()]
    return ampcombi_main, ampcombi_dict_mod

########################################
#  FUNCTION: Parse the GBK/GBFF file 
#########################################
def gbk_parse(gbk_dir, stop_codon_window, ampcombi_dict_mod, transporter_window,filter_stop_codon, outgbk):
    """
    This parses the gbk files in two steps:
    
    Step1: Extracts whether a 'transporter' is present in the vicinity of the hit.
           This is around the hit -10 and +10 CDSs (type = gene is not considered!)
    
    Step2: Extracts whether there is a stop_codon in the vicinity of the CDS hit
           This is around the hit -20 and +20 codons.
           These windows can be changed accordingly by changing the respective parameters.
    """    
    listdict = []
    dict = {}
    
    # STOP codons(+1)
    stop_plus = ['TAG', 'TGA', 'TAA']
    # STOP codons(-1)
    stop_minus = ['ATC', 'ACT', 'ATT']
    # Window size for searching before the start and end positions
    window_size = stop_codon_window
    # Step size for sliding the window
    window_step = 3 
    
    # look for the file name with correct ampcombi sample
    for d in ampcombi_dict_mod:
        for record in SeqIO.parse(gbk_dir, "genbank"):
        
            #########
            #  Step1: Grab the vicinity CDSs and search for any 'transporter' gene
            #########
            # count the number of features in the record that contains the hit/CDS
            cds_count = sum(1 
                            for feature in record.features 
                            if feature.type == "CDS" and "locus_tag" in feature.qualifiers)
            # grab only the features with CDS and index from 0--
            cds_features = [feature 
                            for feature in record.features 
                            if feature.type == "CDS" and "locus_tag" in feature.qualifiers]
            # iterate over every feature in cds_features
            for i, feature in enumerate(cds_features):
                locus_tag = feature.qualifiers["locus_tag"][0]
                # check if locus_tag is present in the list of contig_ids in dict.
                if locus_tag in d['contig_id']:
                    # grab the index of the locus_tag in cds_features
                    index = cds_features.index(feature)
                    # calculate the start and end indices for the preceding and following CDSs = 10
                    start_index = max(0, index - transporter_window)
                    end_index = min(index + transporter_window, cds_count)
                    # grab the preceding and following CDS from the matching index
                    preceding_cds_features = cds_features[start_index:index]
                    following_cds_features = cds_features[index:end_index]
                    # check if the 'transporter' is a product in any of the preceding and following CDSs
                    for f in preceding_cds_features + following_cds_features:
                        if 'product' in f.qualifiers and 'transporter' in f.qualifiers['product'][0]:
                             dict['transporter_protein'] = f.qualifiers['product'][0]
                        else:
                            dict['transporter_protein'] = 'no'
                    
            #########
            #  Step2: Grab the contig_ids - CDS start and end - look for stop codon in vicinity
            #########
            for feature in record.features:
                if feature.type == "CDS" and "locus_tag" in feature.qualifiers:
                    locus_tag = feature.qualifiers["locus_tag"][0]
                    # for every CDS in the list of contig_ids check if its present in the features and add contig name and CDS location and direction
                    for item in d['contig_id']:
                        if locus_tag == item:
                            # add the contig_id and _name
                            dict['name'] = d['name']
                            dict['contig_name'] = record.id
                            dict['contig_id'] = locus_tag
                            # add the CDS location
                            start = feature.location.start.position 
                            end = feature.location.end.position
                            dict['CDS_start'] = start
                            dict['CDS_end'] = end
                            cds_dir = feature.location.strand
                            dict['CDS_dir'] = cds_dir
                            # add and search for any stop codons in the vicinity of the hit
                            sequence = record.seq
                            # define the search window
                            start_window = max(start - window_size, 0)  # Make sure not to go beyond the start of the sequence
                            end_window = min(end + window_size, len(sequence))  # Make sure not to go beyond the end of the sequence
                            # extract the window sequence before and after the CDS
                            start_window_seq = sequence[start_window:start]
                            end_window_seq = sequence[end:end_window]
                            # flag to check if any of the codons are found : keeps iterating until codon_found
                            codon_found = False
                            # check if the strand 5'-3'
                            if cds_dir == 1:
                                # iterate over the start window sequences first in a window of 3 bases from right to left
                                for i in range(len(start_window_seq), -1, -window_step):
                                    codon = start_window_seq[i:i+3]
                                    if codon in stop_plus:
                                        codon_found = True
                                        #print(f"Codon sequence '{codon}' was found in '{start}'")
                                        dict['CDS_stop_codon_found'] = str(codon)
                                        break # Exit the loop if any of the stop codon is found
                                # if the codon is not found in the start window, check the end window
                                if not codon_found:    
                                    for i in range(len(end_window_seq), -1, -window_step):
                                        codon = end_window_seq[i:i+3]
                                        if codon in stop_plus:
                                            codon_found = True
                                            #print(f"Codon sequence '{codon}' was found in '{start}'")
                                            dict['CDS_stop_codon_found'] = str(codon)
                                            break  # Exit the loop if any of the stop codon is found
                                # if the codon is not found in either the start or end window
                                if not codon_found:
                                    #print(f"Codon sequence '{stop_plus}' not found in the window up and down stream of the hit '{start}'") 
                                    dict['CDS_stop_codon_found'] = 'no' 
                                listdict.append(dict.copy())
                            # check if the strand 3'-5'
                            else:
                                # iterate over the start window sequences first in a window of 3 bases from right to left
                                for i in range(len(start_window_seq), -1, -window_step):
                                    codon = start_window_seq[i:i+3]
                                    if codon in stop_minus:
                                        codon_found = True
                                        #print(f"Codon sequence '{codon}' was found in '{start}'") 
                                        dict['CDS_stop_codon_found'] = str(codon)
                                        break # Exit the loop if any of the stop codon is found
                                # if the codon is not found in the start window, check the end window
                                if not codon_found:    
                                    for i in range(len(end_window_seq), -1, -window_step):
                                        codon = end_window_seq[i:i+3]
                                        if codon in stop_minus:
                                            codon_found = True
                                            #print(f"Codon sequence '{codon}' was found in '{start}'")
                                            dict['CDS_stop_codon_found'] = str(codon)
                                            break  # Exit the loop if any of the stop codon is found
                                # if the codon is not found in either the start or end window
                                if not codon_found:
                                    #print(f"Codon sequence '{stop_minus}' not found in the window up and down stream of the hit '{start}'") 
                                    dict['CDS_stop_codon_found'] = 'no'
                                listdict.append(dict.copy())
            
            # if flagged remove hits with no stop codons found
            if filter_stop_codon == True:
                listdict = [d for d in listdict if d.get('CDS_stop_codon_found') != 'no']
        
            #########
            #  Step3: Extract the new gbks that contain the hits
            #########
            for item in listdict:
                if item['name'] and item['contig_name'] == record.id:
                    name = item['name']
                    contig_name = item['contig_name']
                    new_seq_record = record
                    print(f'writing {name}_{contig_name}.gbk')
                    SeqIO.write(new_seq_record, f'{outgbk}/{name}_{contig_name}.gbk', "genbank")

    return listdict

########################################
#  FUNCTION: Extract contig names, and transporters if present and filter by stop codon presence and add 
#########################################
def gbkparsing(sample_ampcombi_file, gbk_dir, stop_codon_window, transporter_window, filter_stop_codon, outgbk):
    """
    This parses the ampcombi and gbk files, and extracts necessary information like:
    contig_name, CDS_location, stop codon and any transporter found in the vicinity,
    which increases the likelihood that the CDS/hit is an AMP.    
    """

    ampcombi_main, ampcombi_dict_mod = ampcombi_input(sample_ampcombi_file)
    listdict = gbk_parse(gbk_dir, stop_codon_window, ampcombi_dict_mod, transporter_window, filter_stop_codon, outgbk)
    
    # clear the dictionary to clear memory allocated
    ampcombi_dict_mod.clear()                      
    
    # merge ampcombi and the dict_df based on name and contig_id
    # convert the dictionary to a DataFrame
    dict_df = pd.DataFrame.from_dict(listdict)        
    try:
        merged_df = ampcombi_main.merge(dict_df, on=['name','contig_id'])
    except KeyError as e:
        #specify the error exactly : no sample name and contig id is found due to empty results
        if str(e) == "'name'":
            print("No hits were found in the sample using these parameters, skipping ...")
            merged_df = ampcombi_main
        else:
            import traceback
            traceback.print_exc()
            print(f"KeyError occurred: {e}")
            merged_df = ampcombi_main
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"An error occurred: {e}")
        merged_df = ampcombi_main
    
    # clear the dictionary to clear memory allocated
    listdict.clear()                      

    return merged_df
