#!/bin/python3

#### Title: Extract alignment stats from HMM output txt files ####
# Description: The script parses all HMM files (ending with txt) in the current (or specified) working directory.
# Author: Adapted and modified from Jasmin Frangenberg @jasmezz original script
# Date: October 2023

### Set-up
## Import modules
import argparse, os, re, sys

## Program variables
records = {} # All HMM search records
out = "hmm_stats.csv"
i = 0 # Iterate over records

## Set program description and commandline arguments
d = "This script parses all HMM files (ending with txt). The following alignment stats will be extracted: Query name and accession number, E-value, score and bias (for both full sequence and best domain), exp, N, Sequence, Description.\n To run the script do this: python3 ./script.py -i ./input_dir -o ./output_dir/output.csv"
parser = argparse.ArgumentParser(description=d, add_help=True)
parser.add_argument("-i", "--input", dest="input", nargs='?', help="Specify input directory (relative path). (default: %(default)s", default='./', type=str)
parser.add_argument("-o", "--output", dest="output", nargs='?', help="Specify output file name. Default: " + out)
args = parser.parse_args()

## Get commandline arguments
#input
inp = args.input
#output
out = args.output
outdir = out.split("/")
if len(outdir) > 1:
	outdir = outdir[:len(outdir) - 1]
	outdir = "/".join(outdir)
	os.makedirs(outdir, exist_ok=True)

### Loop over all files
sample = os.path.basename(inp).rsplit('/', 1)[-1].rstrip(".txt")
acc = "NA" # In case no accession number is provided, write NA
with open(inp, "r") as f:

### Parse the current file and extract the following infos: Query, Accession, Scores table with E-values
		for line in f:
			
			line = line.strip() # Get rid of any whitespace characters
			if line.startswith("Query: "): # Get query name of current record
				i += 1
				que = line.split(" ")[7]
				stats = [] # All alignment stats belonging to the record go into this array
			elif line.startswith("Accession: "): # Get accession number of record
				acc = line.split(" ")[3]
			elif line.startswith("//"):
				if stats != []: # Add record only if at least one hit was found
					records[str(i)] = stats
			elif not line.startswith("#"): # Get alignment stats of record
				pat = re.search("\d*\.?\d*e?-?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*e?-?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\s+.+\s+.+", line)
				if pat:
					sta = re.split("\s+", pat.group(0)) # Get list of alignment stats by splitting the match group of the regular expression search group by consecutive whitespaces
					sta = sta[:9] + [" ".join(sta[9:])]
					sta.insert(0, acc)
					sta.insert(0, que)
					stats.append([sample] + sta)

		# Don't forget the last record in file
		if stats != []: # Add record only if at least one hit was found
					records[str(i)] = stats
					
### Write the whole bunch into a csv file
## The format and header should be:
#Query,Accession,E-value,score,bias,E-value,score,bias,exp,N,Sequence,Description,More_description
#,,|---,full sequence,---|,|---,best domain,---|,#domains,---|,,

with open(out, "w") as o:
	o.write("Sample,Query,Accession,Evalue,score,bias,E-value_domain,score,bias,exp,N,Sequence,Description,More_description\n")
	#o.write(",,,|-----,full sequence,-----|,|-----,best domain,-----|,|--- #domains,-----|,,\n")

	for key in sorted(records.keys(), key = int): # Go through all records
		for subject in records[key]: # Go through all subject sequences for this record
			output = (",".join(subject) + "\n")
			o.write(output)
