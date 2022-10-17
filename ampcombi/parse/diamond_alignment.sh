#!/bin/bash

#########################################
#  Diamond alignment
#########################################

INPUT_FASTA=$1
OUTPUT_DIR=$2
REF_DIR=$3

read INPUT_FASTA
read OUTPUT_DIR
read REF_DIR

# Adjust path according to the input folder with the ist of fasta files
IN=$INPUT_FASTA
OUT=$OUTPUT_DIR 
REF_DB=$REF_DIR

diamond blastp \
-p 28 -d $REF_DB/amp_ref -q $IN --quiet \
--outfmt 6 qseqid sseqid pident evalue nident full_qseq full_sseq qseq sseq qcovhsp scovhsp --max-target-seqs 1 --ultra-sensitive -e10000 --id2 1 -s1 -c1 --masking 0 --gapped-filter-evalue 0 --algo 0 --min-score 0 --shape-mask 1111 \
-o $OUT/diamond_matches.txt

echo -e "contig_id\ttarget_id\tpident\tevalue\tnident\tfull_qseq\tfull_sseq\tqseq\tsseq\tqcovhsp\tscovhsp" | cat - $OUT/diamond_matches.txt > $OUT/diamond_matches.tsv

#--log