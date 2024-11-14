#!/bin/python3

from Bio import SeqIO

#Filter assembly faa for contigs of interest (summary_df) and save to out_path.faa
# return path to out_path.faa as input for taxonomic analysis 
def amp_fasta(summary_df, faa, out_path):
    #get contig ids from summary_df 
    arg_contigs = summary_df['contig_id'].drop_duplicates().to_list()
    #filter faa for AMP-contigs of interest for AMP protein analysis (MMSeq)
    faa_sequences = SeqIO.parse(open(faa),'fasta')
    with open(out_path, 'w') as out_file:
        for faa in faa_sequences:
            id, sequence = faa.id, str(faa.seq) 
            for contig in arg_contigs:
                if (id==contig):
                    out_file.write('>'+ id + '\n' + sequence + '\n')
    return out_path # path for input to MMSeq