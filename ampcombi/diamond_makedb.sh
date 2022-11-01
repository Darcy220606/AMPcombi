#!/bin/bash

#########################################
#  Diamond alignment
#########################################

INPUT_FASTA=$1
THREADS=$2
#OUTPUT_DIR=$2

read INPUT_FASTA
read THREADS
#read OUTPUT_DIR

IN=$INPUT_FASTA
P=$THREADS
#OUT=$OUTPUT_DIR 

#cd $OUT
diamond makedb --in $IN -p $P -d amp_ref --quiet