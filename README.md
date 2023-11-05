# AMPcombi : AntiMicrobial Peptides parsing and functional classification tool

<img src="https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/docs/amp-combi-logo.png" width="620" height="200" />

This tool parses the results of antimicrobial peptide (AMP) prediction tools into a single table and aligns the hits against a reference AMP database for functional classifications.

For parsing: AMPcombi is developed to parse the output of these **AMP prediction tools**:
 
| Tool | Version | Link |
| ------------- | ------------- | ------------- |
| Ampir  | 1.1.0  | https://github.com/Legana/ampir |
| AMPlify  | 1.0.3  | https://github.com/bcgsc/AMPlify |
| Macrel  | 1.1.0  | https://github.com/BigDataBiology/macrel |
| HMMsearch  | 3.3.2  | https://github.com/EddyRivasLab/hmmer |
| EnsembleAMPpred  | - | https://pubmed.ncbi.nlm.nih.gov/33494403/ |
| NeuBI  | -  | https://github.com/nafizh/NeuBI |
| AMPgram  | -  | https://github.com/michbur/AmpGram |
| AMPtransformer  | -  | https://github.com/Brendan-P-Moore/AMPTransformer |

For classification: AMPcombi is developed to offer functional annotation of the detected AMPs by alignment to an **AMP reference databases**, for e.g.,:

| Tool | Version | Link |
| ------------- | ------------- | ------------- |
| DRAMP  | 3.0 | https://github.com/CPU-DRAMP/DRAMP-3.0 |

<span style="color:yellow">**Please Note:** In the latest version of AMPcombi, if no database is provided by the user, it will automatically download the [DRAMP db](https://github.com/CPU-DRAMP/DRAMP-3.0)</span>

Alignment to the reference database is done using [diamond blastp v.2.0.15](https://www.nature.com/articles/s41592-021-01101-x)

======================
## Installation
======================

To install AMPcombi:

Install the dependencies of the tool:
- `python` > 3.0
- `biopython`
- `pandas`
- `diamond`
  
Installation can be done using:

 - pip installation
```
pip install AMPcombi
```
 - git repository
 ```
 git clone https://github.com/Darcy220606/AMPcombi.git
 ```
 - conda
```
conda env create -f ampcombi/environment.yml
```
or
```
 conda install -c bioconda AMPcombi
```

======================
## Usage:
======================

There are two basic commands to run AMPcombi:

1. Using `--amp_results`
```console
ampcombi \
--amp_results path/to/my/result_folder/ \
--faa path/to/sample_faa_files/ \
--<tool>_file '.tsv'
```
*<tool> can be ampir, macrel, amplify, hmmsearch, amppred, ampgram, amptransformer*

Here the head folder containing output files has to be given. AMPcombi finds and summarizes the output files from different tools, if the folder is structured  and named as: `/result_folder/toolsubdir/samplesubdir/sample.filetype`. 
 - Note that the filetype ending might vary and can be specified with `--<tool>_file`. For example if the user wants to run ampcombi on results from both ampir and amplify, `--ampir_file '.tsv' --amplify_file '.tsv'` . This creates a dictionary with keys and items. i.e. `'{"key1":"item1", "key2":"item2"}'`
- Note that `--sample_list` can also be given if only specfic samples are needed from the driectory.

The path to the folder containing the respective protein fasta files has to be provided with `--faa`. The files have to be named with `<samplename>.faa`.

Structure of the results folder:

```console
amp_results/
├── tool_1/
|   ├── sample_1/
|   |   └── sample_1.tsv
|   └── sample_2/
|   |   └── sample_2.tsv
├── tool_2/
|   ├── sample_1/
|   |   └── sample_1.txt
|   └── sample_2/
|   |   └── sample_2.txt
├── tool_3/
    ├── sample_1/
    |   └── sample_1.predict
    └── sample_2/
        └── sample_2.predict
```

1. Using `--path_list` and `--sample_list`

```console
ampcombi \
--path_list path_to_sample_1_tool_1.csv path_to_sample_1.csv \
--path_list path_to_sample_2_tool_1.csv path_to_sample_2.csv \
--sample_list sample_1 sample_2 \
--faa path/to/sample_faa_files/
```

Here the paths to the output-files to be summarized can be given by `--path_list` for each sample. Together with this option a list of sample-names has to be supplied.
Either the path to the folder containing the respective protein fasta files has to be provided with `--faa` or, in case of only one sample, the path to the corresponding `.faa` file. The files have to be named with `<samplename>.faa`.


### Input options:
| command | definition | default | example |
| ------------- | ------------- | ------------- | ------------- |
| --amp_results | path to the folder containing different tool's output files | ./test_files/ | ../amp_results/ |
| --sample_list  | list of samples' names | - | sample_1 sample_2 |
| --path_list  | list of paths to output files | - | path_to_sample_1_tool_1.csv path_to_sample_1_tool_1.csv |
| --amp_cutoff  | probability cutoff to filter AMPs by probability (not applicable for hmmsearch) | 0 | 0.5 |
| --hmm_evalue  | probability cutoff to filter AMPs by evalue (only applicable for HMMsearch) | None | 0.05 |
| --db_evalue  | probability cutoff to filter database classifications by evalue -any hit with value below this will have it's database classification removed-| None | 0.05 |
| --aminoacid_length  | probability cutoff to filter AMP hits by the length of the amino acid sequence| 100 | 60 |
| --faa  | path to the folder containing the samples`.faa` files or, in case of only one sample, the path to the corresponding `.faa` file. Filenames have to contain the corresponding sample-name, i.e. sample_1.faa | ./test_faa/ | ./faa_files/|
| --metadata  | path to the folder containing samples metadata `.tsv`. The first column must contain the sample names. | None | ./metadata.tsv/ |
| --tooldict | dictionary of AMP-tools and their respective output file endings | '{"ampir":"ampir.tsv", "amplify":"amplify.tsv", "macrel":"macrel.tsv", "hmmer_hmmsearch":"hmmsearch.txt", "ensembleamppred":"ensembleamppred.txt", "ampgram":"ampgram.tsv", "amptransformer":"amptransformer.txt"}' | - |
| --amp_database | path to the folder containing the reference database files: (1) a fasta file with <.fasta> file extension and (2) the corresponding table with functional and taxonomic classifications in <.tsv> file extension | [DRAMP 'general amps'](http://dramp.cpu-bioinfor.org/downloads/) database | ./amp_ref_database/ |
| --complete_summary | concatenates all samples' summarized tables into one and generates both 'csv' and interactive 'html' files | False | True |
| --log  | print messages into log file instead of stdout | False | True |
| --threads  | adjust the number of threads required for DIAMOND alignemnt depending on the computing resources available  | 4 | 32 |
| --version  | print the version number into stdout | - | 0.1.4 |

 - Note: The fasta file corresponding to the AMP database should not contain any characters other than ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
  - Note: The reference database table should be tab delimited.
  - Note: Further details on teh samples (i.e. metadata) can be optionally added.

### Output:
The output will be written into your working directory, containing the following files and folders:
```console
<pwd>/
├── amp_ref_database/
|   ├── amp_ref.dmnd
|   ├── general_amps_<DATE>_clean.fasta
|   └── general_amps_<DATE>.tsv
├── sample_1/
|   ├── sample_1_amp.faa
|   ├── sample_1_ampcombi.tsv
|   ├── sample_1_diamond_matches.txt
|   └── sample_1_ampcombi.log
├── sample_2/
|   ├── sample_2_amp.faa
|   ├── sample_2_ampcombi.tsv
|   ├── sample_2_diamond_matches.txt
|   └── sample_2_ampcombi.log
├── AMPcombi_summary.tsv
├── AMPcombi_summary.html
└── Ampcombi.log
```

======================
## Contribution:
======================

AMPcombi is a tool developed for parsing results from published AMP prediction tools. We therefore welcome fellow contributors who would like to add new AMP prediction tools results for parsing and alignment.

### Adding a new tool to AMPcombi
In `ampcombi/reformat_tables.py`
- add a new tool function to read the output to a pandas dataframe and return two columns named `contig_id` and `prob_<toolname>`
- add the new function to the `read_path` function


In `ampcombi/ampcombi.py`
- add a new parameter equivalent to the tool `<--tool_file>`.


======================


**Authors**: @louperelo and @darcy220606
