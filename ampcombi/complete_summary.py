#!/usr/bin/env python3


#### TOD: no create subcommand line :
# EG.:
# ampcombi parse_tables
# ampcombi dreplicate or combine_tables
# ampcombi visualise

#### check this  NOT TOO DIFFICULT :
If you want to create a command-line tool with two subcommands, "mmseqs cluster" and "mmseqs align," you can use the `argparse` module in Python. Here's an example of how you can implement it:

```python
import argparse

def cluster(args):
    # Implementation of mmseqs cluster command
    print("Running mmseqs cluster...")
    print("Arguments:", args)

def align(args):
    # Implementation of mmseqs align command
    print("Running mmseqs align...")
    print("Arguments:", args)

# Create the top-level parser
parser = argparse.ArgumentParser(prog='pipi')

# Create subparsers for the cluster and align commands
subparsers = parser.add_subparsers(dest='command')

# Create parser for the 'cluster' command
cluster_parser = subparsers.add_parser('cluster')
cluster_parser.add_argument('input', help='Input file')

# Create parser for the 'align' command
align_parser = subparsers.add_parser('align')
align_parser.add_argument('input', help='Input file')
align_parser.add_argument('--param', help='Additional parameter')

# Parse the command-line arguments
args = parser.parse_args()

# Call the appropriate subcommand function based on the 'command' argument
if args.command == 'cluster':
    cluster(args)
elif args.command == 'align':
    align(args)
```

In this example, the `argparse` module is used to define the command-line interface. The `cluster` and `align` functions represent the implementations of the "mmseqs cluster" and "mmseqs align" commands, respectively. You can replace the print statements in these functions with the actual logic you want to execute for each command.

To run the tool, you can execute the Python script and provide the command and corresponding arguments. For example:

```
$ python pipi.py cluster input_file.fasta
```

This will execute the `cluster` function with the provided input file.

```
$ python pipi.py align input_file.fasta --param value
```

This will execute the `align` function with the provided input file and additional parameter.

Make sure to replace the placeholders (`cluster(args)` and `align(args)`) with the actual implementation logic for each command.




def complete_summary(param_combine_results_pipeline_mode,sample_name,complete_summary_df, sample_summary_df):
    # if pipeline mode is OFF:
    if param_combine_results_pipeline_mode = False:
        #if (complete_summary):
        # concatenate the sample summary to the complete summary and overwrite it
        complete_summary_df = pd.concat([complete_summary_df, sample_summary_df], ignore_index=True)
        complete_summary_df.to_csv('AMPcombi_summary.tsv', sep='\t', index=False)
        #html_generator() #remove now that we will have the server
        #move this to ampcombi main.py:
        #print(f'\n FINISHED: Appended {sample_name} summary file to complete AMPcombi_summary.tsv (and folder AMPcombi_interactive_summary/ were) and saved to your current working directory.')
        print(f'\n FINISHED: Added {sample_name} to complete AMPcombi_summary.tsv.)')
        # CLUSTER
    # if pipeline mode is ON: e.g. in the case of FUNCSCAN
    else:
        

