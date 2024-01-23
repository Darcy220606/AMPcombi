#!/usr/bin/env python3

import os
import sys
import argparse
import warnings
from contextlib import redirect_stdout
from version import __version__
import json
import os.path
import shutil

# import functions from sub-scripts to main:
from reformat_tables import *
from amp_fasta import *
from check_input import *
from amp_database import *
from clustering_hits import *
from print_header import *
from functionality import *
from optional_inputs import *
from parse_gbks import *

# Define input arguments:
parser = argparse.ArgumentParser(prog = 'ampcombi', formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=('''\
    .............................................................................
                                    *AMP-combi*
    .............................................................................
                This tool parses the results of amp prediction tools 
    and aligns the hits against reference databases for functional classification.
            For detailed usage documentation please refer to <github_repo>
    .............................................................................'''),
                                epilog='''Thank you for running AMP-combi!''',
                                add_help=True)

parser.add_argument("--amp_results", dest="amp", nargs='?', help="Enter the path to the folder that contains the different tool's output files in sub-folders named by sample name. \n If paths are to be inferred, sub-folders in this results-directory have to be organized like '/amp_results/toolsubdir/samplesubdir/tool.sample.filetype' \n (default: %(default)s)",
                    type=str, default='./test_files/')
parser.add_argument("--sample_list", dest="samples", nargs='*', help="Enter a list of sample-names, e.g. sample_1 sample_2 sample_n. \n If not given, the sample-names will be inferred from the folder structure",
                    default=[])
parser.add_argument("--path_list", dest="files", nargs='*', action='append', help="Enter the list of paths to the files to be summarized as a list of lists, e.g. --path_list path/to/my/sample1.ampir.tsv path/to/my/sample1.amplify.tsv --path_list path/to/my/sample2.ampir.tsv path/to/my/sample2.amplify.tsv. \n If not given, the file-paths will be inferred from the folder structure",
                    default=[])
parser.add_argument("--amp_cutoff", dest="p", help="Enter the probability cutoff for AMPs for all tools except for HMMsearch \n (default: %(default)s)",
                    type=float, default=0)
parser.add_argument("--hmm_evalue", dest="hmmevalue", help="Enter the evalue cutoff for AMPs for HMMsearch)  \n (default: %(default)s)",
                    type=float, default=None)
parser.add_argument("--db_evalue", dest="dbevalue", help="Enter the evalue cutoff for AMPs for the database diamond alignment. Any evalue below this value will only remove the DRAMP classification and not the entire hit \n (default: %(default)s)",
                    type=float, default=0.05)
parser.add_argument("--aminoacid_length", dest="length", help="Enter the length of the aa sequences required. Any hits below that cutoff will be removed \n (default: %(default)s)",
                    type=int, default=100)
parser.add_argument("--window_size_stop_codon", dest="stopwindowsize", help="Enter the length of the window size required to look for stop codons downstream and upstream of the CDS hits. \n (default: %(default)s)",
                    type=int, default=60)
parser.add_argument("--window_size_transporter", dest="transporterwindowsize", help="Enter the length of the window size required to look for a 'transporter' e.g. ABC transporter downstream and upstream of the CDS hits. \n (default: %(default)s)",
                    type=int, default=11)
parser.add_argument("--remove_stop_codons", dest="removestops", help="Removes any hits/CDSs that dont have a stop codon found in the window below or upstream of the CDS assigned by '--window_size_stop_codon'. Must be turned on if hits are to be removed. \n (default: %(default)s)",
                    type=bool, default=False)
parser.add_argument("--faa", dest="faa", help="Enter the path to the folder containing the reference .faa files or to one .faa file (running only one sample). Filenames have to contain the corresponding sample-name, i.e. sample_1.faa \n (default: %(default)s)",
                    type=str, default='./test_faa/')
parser.add_argument("--gbk", dest="gbk", help="Enter the path to the folder containing the reference .gbk/.gbff files or to one .gbk/.gbff file (running only one sample). Filenames have to contain the corresponding sample-name, i.e. sample_1.gbk/ sample_1.gbff \n (default: %(default)s)",
                    type=str, default='./test_gbff/')
parser.add_argument("--ampir_file", dest="ampir", nargs='?', help="If AMPir was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'ampir.tsv'",
                    type=str, default=None)
parser.add_argument("--amplify_file", dest="amplify", nargs='?', help="If AMPlify was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'amplify.tsv'",
                    type=str, default=None)
parser.add_argument("--macrel_file", dest="macrel", nargs='?', help="If Macrel was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'macrel.tsv'",
                    type=str, default=None)
parser.add_argument("--neubi_file", dest="neubi", nargs='?', help="If Neubi was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'neubi.fasta'",
                    type=str, default=None)
parser.add_argument("--hmmsearch_file", dest="hmmsearch", nargs='?', help="If HMMer/HMMsearch was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'hmmsearch.txt'",
                    type=str, default=None)
parser.add_argument("--ensemblamppred_file", dest="amppred", nargs='?', help="If EnsemblAMPpred was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'ensembleamppred.txt'",
                    type=str, default=None)
parser.add_argument("--ampgram_file", dest="ampgram", nargs='?', help="If AMPgram was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'ampgram.txt'",
                    type=str, default=None)
parser.add_argument("--amptransformer_file", dest="amptransformer", nargs='?', help="If AMPtransformer was used, enter the ending of the input files (as they appear in the directory tree), e.g. 'amptransformer.txt'",
                    type=str, default=None)
#parser.add_argument("--tooldict", dest="tools", help="Enter a dictionary of the AMP-tools used with their output file endings (as they appear in the directory tree), \n Tool-names have to be written as in default: \n (default: %(default)s)",
#                    type=str, default='{"ampir":"ampir.tsv", "amplify":"amplify.tsv", "macrel":"macrel.tsv", "neubi":"neubi.fasta", "ampgram":"ampgram.tsv", "amptransformer":"amptransformer.txt", "hmmer_hmmsearch":"hmmsearch.txt", "ensembleamppred":"ensembleamppred.txt"}')
parser.add_argument("--amp_database", dest="ref_db", nargs='?', help="Enter the path to the folder containing the reference database files (.fa and .tsv); a fasta file and the corresponding table with functional and taxonomic classifications. \n (default: DRAMP database)",
                    type=str, default=None)
parser.add_argument("--complete_summary", dest="complete", nargs='?', help="Concatenates all sample summaries to one final summary and outputs a tsv file",
                    type=bool, default=False)
parser.add_argument("--cluster_hits", dest="cluster", nargs='?', help="Clusters the amp hits ONLY if complete summary is activated. Clusztering is done by MMseqs2 cluster",
                    type=bool, default=False)
parser.add_argument("--cluster_cov_mode", dest="mmseqscovmode", nargs='?', help="This assigns the cov. mode to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs at https://mmseqs.com/latest/userguide.pdf",
                    type=int, default=0)
parser.add_argument("--cluster_mode", dest="mmseqsclustermode", nargs='?', help="This assigns the cluster mode to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs at https://mmseqs.com/latest/userguide.pdf",
                    type=int, default=1)
parser.add_argument("--cluster_coverage", dest="mmseqscoverage", nargs='?', help="This assigns the coverage to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs at https://mmseqs.com/latest/userguide.pdf",
                    type=float, default=0.8)
parser.add_argument("--cluster_seq_id", dest="mmseqsseqid", nargs='?', help="This assigns the seqsid to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs at https://mmseqs.com/latest/userguide.pdf",
                    type=float, default=0.4)
parser.add_argument("--cluster_sensitivity", dest="mmseqssensitivity", nargs='?', help="This assigns sensitivity of alignment to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs at https://mmseqs.com/latest/userguide.pdf",
                    type=float, default=4.0)
parser.add_argument("--cluster_remove_singletons", dest="removesingletons", nargs='?', help="This removes any hits that did not form a cluster",
                    type=bool, default=True)
parser.add_argument("--cluster_retain_label", dest="retainlabels", nargs='?', help="This removes any cluster that only has a certain label in the sample name. For example if you have samples labels with 'S1_metaspades' and 'S1_megahit', you can retain clusters that have samples with suffix '_megahit' by running '--retain_clusters_label megahit'",
                    type=str, default='megahit')
parser.add_argument("--cluster_min_member", dest="minnumber", nargs='?', help="This removes any cluster that has a hit number lower than this",
                    type=int, default=3)
parser.add_argument("--sample_metadata", dest="samplemetadata", help="Path to a tsv-file containing sample metadata, e,g, 'path/to/sample_metadata.tsv'. The metadata table can have more information for sample identification that will be added to the output summary. The table needs to contain the sample names in the first column. \n (default: %(default)s)",
                    type=str, default=None)
parser.add_argument("--contig_metadata", dest="contigmetadata", help="Path to a tsv-file containing contig metadata, e,g, 'path/to/contig_metadata.tsv'. The metadata table can have more information for contig classification that will be added to the output summary. The table needs to contain the sample names in the first column and the contig_ID in the second column. This can be the output from MMseqs2, pydamage and MetaWrap. \n (default: %(default)s)",
                    type=str, default=None)
parser.add_argument("--log", dest="log_file", nargs='?', help="Silences the standard output and captures it in a log file)",
                    type=bool, default=False)
parser.add_argument("--threads", dest="cores", nargs='?', help="Changes the threads used for DIAMOND alignment (default: %(default)s)",
                    type=int, default=4)
parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

# get command line arguments
args = parser.parse_args()

# assign input arguments to variables
path = args.amp
samplelist_in = args.samples
filepaths_in = args.files
p = args.p
aa_len = args.length
dbevalue = args.dbevalue
hmmevalue = args.hmmevalue
faa_path = args.faa
gbk_dir = args.gbk
stop_codon_window = args.stopwindowsize 
filter_stop_codon = args.removestops
transporter_window = args.transporterwindowsize
#tooldict = json.loads(args.tools)
ampir_file = args.ampir
amplify_file = args.amplify
macrel_file = args.macrel
neubi_file = args.neubi
ampgram_file = args.ampgram
amptransformer_file = args.amptransformer
hmmer_file = args.hmmsearch
amppred_file = args.amppred
database = args.ref_db
complete_summary = args.complete
clustering = args.cluster
cov_mod = args.mmseqscovmode
cluster_mode = args.mmseqsclustermode
coverage = args.mmseqscoverage
seq_id = args.mmseqsseqid
sensitivity = args.mmseqssensitivity
remove_singletons = args.removesingletons
min_cluster_members = args.minnumber
retain_clusters_with = args.retainlabels
add_samplemetadata = args.samplemetadata
add_contigmetadata = args.contigmetadata
threads = args.cores

# additional variables
tooldict = dict()
#check fileending input (at least one has to be given)
if ampir_file is not None:
    tooldict['ampir'] = ampir_file
if amplify_file is not None:
    tooldict['amplify'] = amplify_file
if macrel_file is not None:
    tooldict['macrel'] = macrel_file
if neubi_file is not None:
    tooldict['neubi'] = neubi_file
if ampgram_file is not None:
    tooldict['ampgram'] = ampgram_file
if amptransformer_file is not None:
    tooldict['amptransformer'] = amptransformer_file
if hmmer_file is not None:
    tooldict['hmmer_hmmsearch'] = hmmer_file
if amppred_file is not None:
    tooldict['ensembleamppred'] = amppred_file
# extract list of tools from input dictionary. If not given, default dict contains all possible tools
tools = [key for key in tooldict]
# extract list of tool-output file-endings. If not given, default dict contains default endings.
fileending = [val for val in tooldict.values()]
#print(tooldict)
#print(tools, fileending)

# supress panda warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
# supress bipython warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

#########################################
# MAIN FUNCTION
#########################################
def main_workflow():
    # print AMPcombi header
    print_header()
    # check input sample-list and create sample-list if input empty
    samplelist = check_samplelist(samplelist_in, tools, path)
    # check input parameters
    check_input_complete(path, samplelist, filepaths_in, tools)
    # check input filepaths and create list of list of filepaths per sample if input empty
    filepaths = check_pathlist(filepaths_in, samplelist, fileending, path)
    # check amp_ref_database filepaths and create a directory if input empty
    db = check_ref_database(database)
    # initiate a final_summary dataframe to concatenate each new sample-summary
    if (complete_summary):
        complete_summary_df = pd.DataFrame([])

    # generate summary for each sample
    amp_faa_paths = []
    create_diamond_ref_db(db,threads)
    for i in range(0, len(samplelist)):
        # create log for each sample if true
        if args.log_file:
            with open(f'{samplelist[i]}_ampcombi.log', 'w') as f:
                with redirect_stdout(f):
                    print_header()
                    main_list = []
                    print('\n ########################################################## ')
                    print(f'Processing AMP-files from sample: {samplelist[i]}')
                    os.makedirs(samplelist[i], exist_ok=True)
                    # fill main_list with tool-output filepaths for sample i
                    read_path(main_list, filepaths[i], p, hmmevalue, tooldict, faa_path, samplelist[i])
                    # get the path to the samples' corresponding faa file
                    faa_name = check_faa_path(faa_path, samplelist[i])
                    # use main_list to create the summary file for sample i
                    summary_df = summary(main_list, samplelist[i], faa_name, aa_len)
                    # Generate the AMP-faa.fasta for sample i
                    out_path = samplelist[i] +'/'+samplelist[i]+'_amp.faa'
                    amp_fasta(summary_df, faa_name, out_path)
                    amp_faa_paths.append(out_path)
                    print(f'The fasta containing AMP sequences for {samplelist[i]} was saved to {samplelist[i]}/ \n')
                    amp_matches = samplelist[i] +'/'+samplelist[i]+'_diamond_matches.txt'
                    print(f'The diamond alignment for {samplelist[i]} in progress ....')
                    diamond_df = diamond_alignment(db, amp_faa_paths, amp_matches, threads, dbevalue)
                    print(f'The diamond alignment for {samplelist[i]} was saved to {samplelist[i]}/.')
                    # Merge summary_df and diamond_df
                    sample_summary_df = pd.merge(summary_df, diamond_df, on = 'contig_id', how='left')
                    # Insert column with sample name on position 0
                    sample_summary_df.insert(0, 'name', samplelist[i])
                    # Remove the temp directory
                    shutil.rmtree('./temp')
                    # Estimate the aa functions: chemical and physical
                    print(f'The estimation of functional and structural properties for {samplelist[i]} in progress ....')
                    sample_summary_df_functions = functionality(sample_summary_df)
                    sample_summary_df = sample_summary_df_functions
                    # Add contig_ids and filter by stop codon presence and extract new gbks
                    outgbk = samplelist[i] + '/contig_gbks'
                    # Create the new gbks dir
                    os.makedirs(outgbk, exist_ok=True)
                    print(f'Parsing of the corresponding genebank file for {samplelist[i]} in progress ....')
                    sample_summary_df = gbkparsing(sample_summary_df, gbk_dir, stop_codon_window, transporter_window, filter_stop_codon, outgbk)
                    # Merge sample metadata if present
                    sample_metadata_df = sample_metadata_addition(sample_summary_df, add_samplemetadata)
                    sample_summary_df = sample_metadata_df
                    # Merge contig metadata if present
                    contig_metadata_df = contig_metadata_addition(sample_summary_df, add_contigmetadata)
                    sample_summary_df = contig_metadata_df
                    # Fix the column names to match other summary files 
                    sample_summary_df.rename(columns={'name': 'sample_id', 'contig_id':'CDS_id', 'contig_name':'contig_id' }, inplace=True)
                    # Write sample summary into sample output folder
                    sample_summary_df.to_csv(samplelist[i] +'/'+samplelist[i]+'_ampcombi.tsv', sep='\t', index=False)
                    print(f'The summary file for {samplelist[i]} was saved to {samplelist[i]}/.')
                    # Write the log file in the respective sample directory
                    shutil.move(f'{samplelist[i]}_ampcombi.log', samplelist[i] + '/' + samplelist[i]+'_ampcombi.log')
        else:
            main_list = []
            print('\n ########################################################## ')
            print(f'Processing AMP-files from sample: {samplelist[i]}')
            os.makedirs(samplelist[i], exist_ok=True)
            # fill main_list with tool-output filepaths for sample i
            read_path(main_list, filepaths[i], p, hmmevalue, tooldict, faa_path, samplelist[i])
            # get the path to the samples' corresponding faa file
            faa_name = check_faa_path(faa_path, samplelist[i])
            # use main_list to create the summary file for sample i
            summary_df = summary(main_list, samplelist[i], faa_name, aa_len)
            # Generate the AMP-faa.fasta for sample i
            out_path = samplelist[i] +'/'+samplelist[i]+'_amp.faa'
            amp_fasta(summary_df, faa_name, out_path)
            amp_faa_paths.append(out_path)
            print(f'The fasta containing AMP sequences for {samplelist[i]} was saved to {samplelist[i]}/ \n')
            amp_matches = samplelist[i] +'/'+samplelist[i]+'_diamond_matches.txt'
            print(f'The diamond alignment for {samplelist[i]} in progress ....')
            diamond_df = diamond_alignment(db, amp_faa_paths, amp_matches, threads, dbevalue)
            print(f'The diamond alignment for {samplelist[i]} was saved to {samplelist[i]}/.')
            # Merge summary_df and diamond_df
            sample_summary_df = pd.merge(summary_df, diamond_df, on = 'contig_id', how='left')
            # Insert column with sample name on position 0
            sample_summary_df.insert(0, 'name', samplelist[i])
            # Remove the temp directory
            shutil.rmtree('./temp')
            # Estimate the aa functions: chemical and physical
            print(f'The estimation of functional and structural properties for {samplelist[i]} in progress ....')
            sample_summary_df_functions = functionality(sample_summary_df)
            sample_summary_df = sample_summary_df_functions
            # Add contig_ids and filter by stop codon presence and extract new gbks
            outgbk = samplelist[i] + '/contig_gbks'
            # Create the new gbks dir
            os.makedirs(outgbk, exist_ok=True)
            print(f'Parsing of the corresponding genebank file for {samplelist[i]} in progress ....')
            sample_summary_df = gbkparsing(sample_summary_df, gbk_dir, stop_codon_window, transporter_window, filter_stop_codon, outgbk)
            # Merge sample metadata if present
            metadata_df = sample_metadata_addition(sample_summary_df, add_samplemetadata)
            sample_summary_df = sample_metadata_df
            # Merge contig metadata if present
            contig_metadata_df = contig_metadata_addition(sample_summary_df, add_contigmetadata)
            sample_summary_df = contig_metadata_df
            # Fix the column names to match other summary files 
            sample_summary_df.rename(columns={'name': 'sample_id', 'contig_id':'CDS_id', 'contig_name':'contig_id' }, inplace=True)
            # Write sample summary into sample output folder
            sample_summary_df.to_csv(samplelist[i] +'/'+samplelist[i]+'_ampcombi.tsv', sep='\t', index=False)
            print(f'The summary file for {samplelist[i]} was saved to {samplelist[i]}/.tsv')
        
        if (complete_summary):
        # concatenate the sample summary to the complete summary and overwrite it
            complete_summary_df = pd.concat([complete_summary_df, sample_summary_df], ignore_index=True)
            complete_summary_df.to_csv('AMPcombi_summary.tsv', sep='\t', index=False)
            print(f'\n DONE: Appended {samplelist[i]} summary file to complete AMPcombi_summary.tsv and saved to your current working directory.')
        else: 
            print(f'\n DONE: AMPcombi created summaries for all input samples.')
    
    # Only if complete summary and clustering is activated:
    if (complete_summary and clustering):
        print(f'\n All hits in the AMPcombi_summary.tsv will now be clustered by MMSeqs2')
        merged_df = parsing_input_for_cluster(complete_summary_df)
        mmseqs_cluster(cov_mod,cluster_mode,coverage,seq_id,sensitivity,threads)
        complete_summary_clusters_df = compile_clusters(merged_df,retain_clusters_with,remove_singletons,min_cluster_members)
        complete_summary_clusters_df.to_csv(f'AMPcombi_summary_clusters.tsv', sep='\t', index=False)
        print(f'\n DONE: AMPcombi_summary_clusters.tsv was saved to your current working directory.')
        print('\n ########################################################## ')

#########################################
# LOG FUNCTION
#########################################
def main():
    if (args.log_file == True and not os.path.exists('Ampcombi.log')):
        with open(f'Ampcombi.log', 'w') as f:
            with redirect_stdout(f):
                main_workflow()
    elif(args.log_file == True and os.path.exists('Ampcombi.log')):
        with open(f'Ampcombi.log', 'a') as f:
            with redirect_stdout(f):
                main_workflow()
    else: main_workflow()
    
if __name__ == "__main__":
    main()