#!/bin/bash

#########################################
#  Diamond alignment
#########################################

INPUT_FASTA=$1
#OUTPUT_DIR=$2

read INPUT_FASTA
#read OUTPUT_DIR

IN=$INPUT_FASTA
#OUT=$OUTPUT_DIR 

#cd $OUT
diamond makedb --in $IN -p 28 -d amp_ref --quiet