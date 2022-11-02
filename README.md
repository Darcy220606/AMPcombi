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

For classification: AMPcombi is developed to offer functional annotation of the detected AMPs by alignment to an **AMP reference databases**, for e.g.,:

| Tool | Version | Link |
| ------------- | ------------- | ------------- |
| DRAMP  | 3.0 | https://github.com/CPU-DRAMP/DRAMP-3.0 |

Alignment to the reference database is done using [diamond blastp v.2.0.15](https://www.nature.com/articles/s41592-021-01101-x)

======================
## Installation
======================

To install AMPcombi:

Add dependencies of the tool; `python` > 3.0, `biopython`, `pandas` and `diamond`.
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
--faa_folder path/to/sample_faa_files/
```

Here the head folder containing output files has to be given. AMPcombi finds and summarizes the output files from different tools, if the folder is structured  and named as: `/result_folder/toolsubdir/samplesubdir/sample.tool.filetype`. 
 - Note that the filetype ending might vary and can be specified with `--tooldict`, if it is different from the default. When passing a dictionary via command line, this has to be done as a string with single quotes `' '` and the dictionary keys and items with double quotes `" "`. i.e. `'{"key1":"item1", "key2":"item2"}'`
- Note that `--sample_list` can also be given if only specfic samples are needed from the driectory.

The path to the folder containing the respective protein fasta files has to be provided with `--faa_folder`. The files have to be named with `<samplename>.faa`.

Structure of the results folder:

```console
amp_results/
├── tool_1/
|   ├── sample_1/
|   |   └── sample_1.tool_1.tsv
|   └── sample_2/
|   |   └── sample_2.tool_1.tsv
├── tool_2/
|   ├── sample_1/
|   |   └── sample_1.tool_2.txt
|   └── sample_2/
|   |   └── sample_2.tool_2.txt
├── tool_3/
    ├── sample_1/
    |   └── sample_1.tool_3.predict
    └── sample_2/
        └── sample_2.tool_3.predict
```

2. Using `--path_list` and `--sample_list`

```console
ampcombi \
--path_list path_to_sample_1_tool_1.csv path_to_sample_1_tool_1.csv \
--path_list path_to_sample_2_tool_1.csv path_to_sample_2_tool_1.csv \
--sample_list sample_1 sample_2 \
--faa_folder path/to/sample_faa_files/
```

Here the paths to the output-files to be summarized can be given by `--path_list` for each sample. Together with this option a list of sample-names has to be supplied.
The path to the folder containing the respective protein fasta files has to be provided with `--faa_folder`. The files have to be named with `<samplename>.faa`.


### Input options:
| command | definition | default | example |
| ------------- | ------------- | ------------- | ------------- |
| --amp_results | path to the folder containing different tool's output files | ./test_files/ | ../amp_results/ |
| --sample_list  | list of samples' names | - | sample_1 sample_2 |
| --path_list  | list of paths to output files | - | path_to_sample_1_tool_1.csv path_to_sample_1_tool_1.csv |
| --cutoff  | probability cutoff to filter AMPs | 0 | 0.5 |
| --faa_folder  | path to the folder containing the samples` .faa files, Filenames have to contain the corresponding sample-name, i.e. sample_1.faa | ./test_faa/ | ./faa_files/|
| --tooldict | dictionary of AMP-tools and their respective output file endings | '{"ampir":"ampir.tsv", "amplify":"amplify.tsv", "macrel":"macrel.tsv", "hmmer_hmmsearch":"hmmsearch.txt", "ensembleamppred":"ensembleamppred.txt"}' | - |
| --amp_database | path to the folder containing the reference database files: (1) a fasta file with <.fasta> file extension and (2) the corresponding table with with functional and taxonomic classifications in <.tsv> file extension | [DRAMP 'general amps'](http://dramp.cpu-bioinfor.org/downloads/) database | ./amp_ref_database/ |
| --complete_summary | concatenates all samples' summarized tables into one and generates both 'csv' and interactive 'html' files | False | True |
| --log  | print messages into log file instead of stdout | False | True |
| --threads  | adjust the number of threads required for DIAMOND alignemnt depending on the computing resources available  | 4 | 32 |
| --version  | print the version number into stdout | - | 0.1.4 |

 - Note: The fasta file corresponding to the AMP database should not contain any characters other than ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
  - Note: The reference database table should be tab delimited.

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
|   ├── sample_1_ampcombi.csv
|   └── sample_1_diamond_matches.txt
├── sample_2/
|   ├── sample_2_amp.faa
|   ├── sample_2_ampcombi.csv
|   └── sample_2_diamond_matches.txt
├── AMPcombi_summary.csv
├── AMPcombi_summary.html
└── ampcombi.log
```

======================
## Contribution:
======================

AMPcombi is a tool developed for parsing results from published AMP prediction tools. We therefore welcome fellow contributors who would like to add new AMP prediction tools results for parsing and alignment.

### Adding a new tool to AMPcombi
In `ampcombi/reformat_tables.py`
- add a new tool function to read the output to a pandas dataframe and return two columns named `contig_id` and `prob_<toolname>`
- add the new function to the `read_path` function


In `ampcombi/main.py`
- add your default `tool:tool.fileending` to the default of `--tooldict`


======================


**Authors**: @louperelo and @darcy220606
