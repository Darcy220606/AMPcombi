# AMPcombi : Antimicorbial peptides parsing and functional classification tool

# ![Logo](docs/ampcombi_logo_update.png)

This tool parses the results of amp prediction tools into a single table and aligns the hits against a reference database of antimicrobial peptides for functional classifications.

For parsing: AMpcombi is developed to parse the output of these **amp prediction tools**:
 
| Tool | Version | Link |
| ------------- | ------------- | ------------- |
| Ampir  | 1.1.0  | https://github.com/Legana/ampir |
| AMPlify  | 1.0.3  | https://github.com/bcgsc/AMPlify |
| Macrel  | 1.1.0  | https://github.com/BigDataBiology/macrel |
| HMMsearch  | 3.3.2  | https://github.com/EddyRivasLab/hmmer |
| EnsembleAMPpred  | - | https://pubmed.ncbi.nlm.nih.gov/33494403/ |
| NeuBI  | -  | https://github.com/nafizh/NeuBI |

For classification: AMPcombi is developed to offer functional annotation of the detcted AMPs by alignemnt to **AMP reference databases**, for e.g.,:

| Tool | Version | Link |
| ------------- | ------------- | ------------- |
| DRAMP  | 3.0 | https://github.com/CPU-DRAMP/DRAMP-3.0 |

Alignment to the reference database is done using [diamond blastp v.2.0.15](https://www.nature.com/articles/s41592-021-01101-x)

======================
## Installation
======================

To install AMPcombi: 

Add dependencies of the tool; python > 3.0, biopython, pandas and diamond.
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
 conda install AMPcombi
``` 

 and the latest release can be installed directly from pip, conda, docker, this repository, or from the galaxy toolshed:

pip install hAMRonization
PyPI version PyPI downloads

Or

conda create --name hamronization --channel conda-forge --channel bioconda --channel defaults hamronization
version-on-conda conda-download last-update-on-conda

Or to install using docker:

docker pull finlaymaguire/hamronization:latest
Or to install the latest development version:

git clone https://github.com/pha4ge/hAMRonization
pip install hAMRonization
### For pip installation





### Create conda environment


======================
## Usage:
======================

There are two basic commands to run AMPcombi:

1. Using `--amp_results`
```console
ampcombi --amp_results path/to/my/result_folder/
```

Here the head folder containing output files has to be given. AMPcombi finds and summarizes the output files from different tools, if the folder is structured  and named as: `/result_folder/toolsubdir/samplesubdir/sample.tool.filetype`.
 - Note that the filetype ending might vary and can be specified with `--tooldict`, if it is different from the default.


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
ampcombi --path_list [[list of paths to sample_1-outputs][list of paths to sample_2-outputs]] --sample_list [sample_1, sample_2] 
```

Here the paths to the output-files to be summarized can be given as a list for each sample. Together with this option a list of sample-names has to be supplied.


### Input options:
| command | definition | default | example |
| ------------- | ------------- | ------------- | ------------- |
| --amp_results | path to the folder containing different tool's output files | ./test_files/ | ../amp_results |
| --sample_list  | list of samples' names | [] | [sample_1, sample_2] |
| --path_list  | list of paths to output files | [] | [[paths to sample_1 output], [paths to sample_2 outputs]] |
| --outdir  | name of the output directory | ./ampcombi_results/ | ./ampcombi_results/ |
| --cutoff  | probability cutoff to filter AMPs | 0 | 0.5 |
| --faa_folder  | path to the folder containing the samples` .faa files, Filenames have to contain the corresponding sample-name, i.e. sample_1.faa | ./test_faa/ | ./faa_files/|
| --tooldict | dictionary of AMP-tools and their respective output file endings | {'ampir':'ampir.tsv', 'amplify':'amplify.tsv', 'macrel':'macrel.tsv', 'hmmer_hmmsearch':'hmmsearch.txt', 'ensembleamppred':'ensembleamppred.txt'} | - |
| --amp_database | path to the folder containing the reference database files: (1) a fasta file with <.fasta> file extension and (2) the corresponding table with with functional and taxonomic classifications in <.tsv> file extension | [DRAMP 'general amps'](http://dramp.cpu-bioinfor.org/downloads/) database | ./amp_ref_database/ |

 - Note: The fasta file corresponding to the AMP database should not contain any characters other than ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W',',Y']
  - Note: The refernce database table should be tab delimited.


======================
## Contribution:
======================

AMPcombi is a tool developed for parsing results from published AMP prediction tools. We therfore welcome fellow contributers who would like to add new AMP prediction tools results for parsing and alignment.

### Adding a new tool to AMPcombi
In `ampcombi/reformat_tables.py` 
- add a new tool function to read the output to a pandas dataframe
- add the new function to the `read_path` function


In `ampcombi/main.py`
- add your default `tool:tool.fileending`to the default of `--tooldict`


======================


**Authors**: @louperelo and @darcy220606
