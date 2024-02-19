#!/bin/python3

# TITLE: Estimate aminoacid sequence functionality: hydrophobicity, mwt, isoelectric point, etc.

import os
import re
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

########################################
#  FUNCTION: Molecular weight
#########################################
def mwt(row):
    sequence = row['aa_sequence']
    try:
        analysis = ProteinAnalysis(sequence)
        # round to 3 decimal places
        row['molecular_weight'] = round(analysis.molecular_weight(), 3)
    except ValueError:
        row['molecular_weight'] = None
    return row

########################################
#  FUNCTION: isoelectric point
#########################################
def ip(row):
    sequence = row['aa_sequence']
    try:
        analysis = ProteinAnalysis(sequence)
        row['isoelectric_point'] = round(analysis.isoelectric_point(), 3)
    except ValueError:
        row['isoelectric_point'] = None
    return row

########################################
#  FUNCTION: structure prediction
#########################################
def structure(row):
    sequence = row['aa_sequence']
    analysis = ProteinAnalysis(sequence)
    struct = analysis.secondary_structure_fraction()
    row['helix_fraction'] = round(struct[0], 3)
    row['turn_fraction'] = round(struct[1], 3)
    row['sheet_fraction'] = round(struct[2], 3)
    return row

########################################
#  FUNCTION: structure prediction
#########################################
def hydrophobicity(row):
    # calculated according to the Eisenberg scale was used
    hydrophobicity = {'A':  0.620,'R': -2.530,'N': -0.780,'D': -0.900,'C':  0.290,'Q': -0.850,'E': -0.740,'G':  0.480,'H': -0.400,'Y':  0.260,
                          'I':  1.380,'L':  1.060,'K': -1.500,'M':  0.640,'F':  1.190,'P':  0.120,'S': -0.180,'T': -0.050,'W':  0.810,'V':  1.080}
    sequence = row['aa_sequence']
    # check of there are any non of the letters within the aa seq
    if re.search(r'[^ARNDCQEGHYILKMFPSTWVarndcqeghyilkmfpstwv]', sequence):
        row['hydrophobicity'] = 'None'
    else:
        valid_values = []
        for resi in sequence:
            valid_values.append(hydrophobicity[resi])
            row['hydrophobicity'] = round(float(sum(valid_values)), 3)
    return row

########################################
#  FUNCTION: structure prediction
#########################################
def functionality(df):
    table = df.apply(mwt, axis=1)
    table = table.apply(structure, axis=1)
    table = table.apply(ip, axis=1)
    table = table.apply(hydrophobicity, axis=1)
    return table    
