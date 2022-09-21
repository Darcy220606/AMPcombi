import os
import pandas as pd
import argparse
from Bio import SeqIO

#test_df = pd.read_csv('../test_faa/sample_1_testsummary.csv')
test_faa = '../test_faa/sample_1.faa'
faa_seq = SeqIO.parse(open(test_faa), 'fasta')
fasta_df=pd.DataFrame(columns=['contig_id', 'aa_sequence'])
for contig in faa_seq:
    contig_id, sequence = contig.id, str(contig.seq)
    fasta_df = fasta_df.append({'contig_id':contig_id, 'aa_sequence':sequence}, ignore_index=True)
print(fasta_df.shape)
print(fasta_df.head())


def arg_contigs(data, fasta, out_name):
    #filter contigs with antibiotic resistance genes
    arg_contigs = data['contig_name'].drop_duplicates().to_list()
    # filter contig sequence information from original fasta file
    #filter fasta for contigs with antibiotic resistance genes (arg) for taxonomic analysis
    fasta_sequences = SeqIO.parse(open(fasta),'fasta')
    with open(out_name, 'w') as out_file:
        for fasta in fasta_sequences:
            #name, sequence = fasta.id, fasta.seq.tostring() #tostring() should be replaced by str(fasta.seq), but is not working on my computer
            name, sequence = fasta.id, str(fasta.seq) 
            for c in arg_contigs:
                if (name==c):
                    out_file.write('>'+ name + '\n' + sequence + '\n')