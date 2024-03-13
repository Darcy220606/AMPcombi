# AMPcombi : A parsing, filtering and annotation workflow for tools predicting AntiMicrobial Peptide genes

<img src="https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/docs/amp-combi-logo.png" width="300" height="250" /> [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) [![Install with Bioconda](https://anaconda.org/bioconda/ampcombi/badges/downloads.svg)](https://anaconda.org/bioconda/ampcombi) [![Install with Bioconda](https://anaconda.org/bioconda/ampcombi/badges/version.svg)](https://anaconda.org/bioconda/ampcombi)


AMPcombi and its submodules provide a command-line interface to parse the results of antimicrobial peptide (AMP) prediction tools into a single table, aligns the AMP hits against a reference AMP database for functional classifications, filters the AMP hits according to their physicochemical properties, clusters the filtered hits and predicts signaling peptides if present.

---
## Introduction

<span style="color:green"> (A) For parsing and filtering; AMPcombi is developed to parse and filter the output of only the following **AMP prediction tools**:
 
| Tool | Version |
| ------------- | ------------- |
| [Ampir](https://github.com/Legana/ampir) | 1.1.0  |
| [AMPlify](https://github.com/bcgsc/AMPlify) | 1.0.3  |
| [Macrel](https://github.com/BigDataBiology/macrel)  | 1.1.0  | 
| [HMMsearch](https://github.com/EddyRivasLab/hmmer)  | 3.3.2  |
| [EnsembleAMPpred](https://pubmed.ncbi.nlm.nih.gov/33494403/)  | - |
| [NeuBI](https://github.com/nafizh/NeuBI) | -  | 
| [AMPgram](https://github.com/michbur/AmpGram) | -  | 
| [AMPtransformer](https://github.com/Brendan-P-Moore/AMPTransformer) | -  |

<span style="color:green"> (B) For classification; AMPcombi is developed to provide the functional annotation of the detected AMPs by alignment to an **AMP reference databases** using [diamond blastp v.2.0.15](https://www.nature.com/articles/s41592-021-01101-x):

| Tool | Version |
| ------------- | ------------- |
| [DRAMP](https://github.com/CPU-DRAMP/DRAMP-3.0) | 3.0 | 
| [Diamond](https://github.com/bbuchfink/diamond) | 2.0.15 |

⚠️ If no database is provided by the user, AMPcombi will automatically download the [DRAMP db](https://github.com/CPU-DRAMP/DRAMP-3.0) and use the files for classification.

<span style="color:green"> (C) For structural and physical annotations; corresponding to the molecular weight, isoelectric point, hydrophobicity, pH and the fraction of helix turns and beta sheets were calculated using:

| Tool | Version |
| ------------- | ------------- |
| BioPython - ProteinAnalysis | 1.80 | 

 Any transporter gene if present in the vicinity of the AMP on the contig is also reported.

<span style="color:green"> (D) For filtering AMP hits: AMPcombi is developed to filter hits based on their AMP probabilities given by each AMP tool - E-values in case of HMMsearch - and the presence of stop-codons downstream and upstream of the AMP hit on the contig.

 ⚠️ If contig or sample details/information are available this can be added here, for example, metadata, pydamage, contig taxonomic classification, etc.

<span style="color:green"> (E) For clustering AMP hits: AMPcombi is developed to cluster the AMPs predicted and remove singletons and clusters with a minium number of members or clusters with only a specific string in their name. This is done using:

| Tool | Version |
| ------------- | ------------- |
| [MMseqs2](https://github.com/soedinglab/MMseqs2) | 15.6f452 | 

<span style="color:green"> (F) For signal peptide detection: Signal peptides can be predicted from clustered AMPs  given the user installs SignalP separately. For licensing issues, SignalP can only be downloaded and used by academic users; other users are requested to contact DTU Health Technology Software Package before using it. Please refer to [SignalP documentation](https://services.healthtech.dtu.dk/services/SignalP-6.0/). For obtaining SignalP please follow:

| Tool | Version |
| ------------- | ------------- |
| [SignalP](https://services.healthtech.dtu.dk/services/SignalP-6.0/) | 6.0 | 

<span style="color:green"> (F) For visualization: a shiny app is accessible through `./pyshiny`. This is a user interface that renders the AMPcombi summary table and genrates a number of figures for data analysis.

| Tool | Version |
| ------------- | ------------- |
| [Shiny](https://shiny.posit.co/py/) | 0.7.1 | 


---
## Installation

To install AMPcombi:

- Using **conda**:
```
conda create -n ampcombi python==3.11 diamond==2.0.15 mmseqs==15.6f452 ampcombi
```
or 
```
conda env create -f ./ampcombi/environment.yml
```

- Using **singularity and docker**:
```

```

- From git repository:
 ```
 git clone https://github.com/Darcy220606/AMPcombi.git
 ```
---
## Usage and Output:

For full usage documentaion of ampcombi and it's submodules please refer to the help documentation:

```
ampcombi --help
```
 
### Submodules
---

`ampcombi parse_tables`

This subcommand parses and filters the output files generated by different AMP prediction tools, aligns the amino-acid sequences to the reference database and estimates the physicochemical annotations .

To get a full list of options available and their defaults please refer to the help documentation of the submodule

```
ampcombi parse_tables --help
```

Example usage (1): 

```console
ampcombi parse_tables \
--amp_results path/to/my/result_folder/ \
--faa path/to/sample_faa_files/ \
--gbk path/to/sample_gbk_or_gbff_files/ \
--sample_list sample_1 sample_2 \
--contig_metadata path/to/contig_metadata.tsv
--<tool_1>_file '.tsv' 
--<tool_2>_file '.txt' 
--log true
--threads 10
```
In this case we use the `--amp_results` option, to supply AMP tool prediction results from **many samples** in a folder format. The folder must follow this structure:

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
    |   └── sample_1.fasta
    └── sample_2/
        └── sample_2.fasta
```
- **--<tool>_file** the <tool> should be changed to ampir, macrel, amplify, neubi, hmmsearch, ensemblamppred, ampgram, amptransformer. The argument value should be a suffix of the files generated by that tool. We have assigned defaults to each tool, however the user can change those defaults according to the user input file extensions.
- **--contig_metadata** a `*.tsv` file that must contain the sample name in the first column and the contig ID/name in the second column. The column headers are not important.
- **--faa** a folder containing annotated files of the AMP hits with a suffix `*.faa`. This can be generated by any annotation tool *i.e.*, [PROKKA](https://github.com/tseemann/prokka), [PYRODIGAL](https://github.com/althonos/pyrodigal), etc.  **NOTE** The files have to include the sample name, for example, `<samplename>.faa`.
- **--gbk** a folder containing annotated files of the AMP hits with a suffix `*.gbk` or `*.gbff`. This can be generated by any annotation tool *i.e.*, [PROKKA](https://github.com/tseemann/prokka), [PYRODIGAL](https://github.com/althonos/pyrodigal), etc. **NOTE** The files have to include the sample name, for example, `<samplename>.gbk/gbff`.
- **--log** captures all standard output and standard errors in a log file
- **--threads** the number of cores to parallelize the job

Example usage (2):

```console
ampcombi parse_tables \
--path_list path_to_sample_1_tool_1.csv path_to_sample_1_tool_2.txt \
--sample_list sample_1 \
--faa path/to/sample_faa_files/sample_1.faa \
--gbk path/to/sample_gbk_or_gbff_files/sample_1.<gbk><gbff> \
--<tool_1>_file '.tsv' 
--<tool_2>_file '.txt' 
```

In this case we use the `--path_list` option, to supply AMP tool prediction results from a **single sample** in a list format. 

Optional arguments:

| command | definition | default | example |
| ------------- | ------------- | ------------- | ------------- |
| --amp_cutoff  | Probability cutoff to filter AMPs by probability (not applicable for hmmsearch) | 0.0 | 0.5 |
| --hmm_evalue  | Probability cutoff to filter AMPs by E-value (only applicable for HMMsearch) | None | 0.05 |
| --db_evalue  | Probability cutoff to filter database classifications by E-value - any hit with an E-value below this will have it's database classification removed | None | 0.05 |
| --aminoacid_length  | Probability cutoff to filter AMP hits by the length of the amino acid sequence| 100 | 60 |
| --window_size_stop_codon  | The length of the window size required to look for stop codons downstream and upstream of the CDS hits | 60 | 40 |
| --window_size_transporter  | The length of the window size required to look for a 'transporter' e.g. ABC transporter downstream and upstream of the CDS hits | 11 | 20 |
| --remove_stop_codons  | Removes any AMP hits that don't have a stop codon found in the window downstream or upstream of the CDS assigned by '--window_size_stop_codon'. Must be turned on if hits are to be removed | False | True |
| --sample_metadata  | Path to a tsv-file containing sample metadata, e.g. 'path/to/sample_metadata.tsv'. The metadata table can have more information for sample identification that will be added to the output summary. The table needs to contain the sample names in the first column. | None | ./sample_metadata.tsv/ |
| --contig_metadata  | Path to a tsv-file containing contig metadata, e.g. 'path/to/contig_metadata.tsv'. The metadata table can have more information for contig classification that will be added to the output summary. The table needs to contain the sample names in the first column and the contig_ID in the second column. The metadata table can be the output from MMseqs2, pydamage and MetaWrap. | None | ./contig_metadata.tsv/ |
| --amp_database | Path to the folder containing the reference database files: (1) a fasta file with <.fasta> file extension and (2) the corresponding table with functional and taxonomic classifications in <.tsv> file extension | [DRAMP 'general amps'](http://dramp.cpu-bioinfor.org/downloads/) | ./custom_amp_ref_database/ |

 ⚠️ With regards to the reference database supplied to `--amp_database`: 
  - The fasta file corresponding to the AMP database should not contain any characters other than ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'].
  - The reference database table should be tab delimited.

Output:

The output will be written into your working directory, containing the following files and folders:

```console
<pwd>/
├── amp_ref_database/
|   ├── amp_ref.dmnd
|   ├── general_amps_<DATE>_clean.fasta
|   └── general_amps_<DATE>.tsv
├── sample_1/
|   ├── contig_gbks/
|   ├── sample_1_amp.faa
|   ├── sample_1_ampcombi.tsv
|   ├── sample_1_diamond_matches.txt
|   └── sample_1_ampcombi.log
├── sample_2/
|   ├── contig_gbks/
|   ├── sample_2_amp.faa
|   ├── sample_2_ampcombi.tsv
|   ├── sample_2_diamond_matches.txt
|   └── sample_2_ampcombi.log
└── Ampcombi_parse_tables.log
```
---

`ampcombi complete`

This subcommand concatenates the ampcombi summaries generated by running `ampcombi parse_tables`, in a single complete summary.

To get a full list of options available and their defaults please refer to the help documentation of the submodule

```
ampcombi complete --help
```

Example usage (1): 

```console
ampcombi complete \
--summaries_directory path/to/ampcombi_parse_tables_results_folder/ 
```
In this case we use the `--summaries_directory` option, to supply the samples' result folder from `--ampcombi parse_tables` which should contain the folder structure from `ampcombi parse_tables` in a  parent folder for example named `./ampcombi/...`.

Example usage (2): 

```console
ampcombi complete \
--summaries_files path/to/ampcombi_parse_tables/sample_1_ampcombi.tsv path/to/ampcombi_parse_tables/sample_2_ampcombi.tsv/ 
```
In this case we use the `--summaries_files` option, to supply the `ampcombi_parse_tables` AMPcombi summary files in a  list format.

Output:

The output will be written into your working directory, containing the following files and folders:

```console
<pwd>/
└── Ampcombi_summary.tsv
└── Ampcombi_complete.log
```
---
`ampcombi cluster`

This subcommand clusters the AMPcombi complete summary generated by running `ampcombi complete`. As this uses `mmseqs cluster` in the background some parameters that were deemed important for the purpose of AMPcombi were incorporated as optional arguments.

To get a full list of options available and their defaults please refer to the help documentation of the submodule

```
ampcombi cluster --help
```

Example usage: 
```console
ampcombi cluster \
--ampcombi_summary path/to/Ampcombi_summary.tsv 
```
The input option `--ampcombi_summary` requires the 'Ampcombi_summary.tsv' which is generated by the `ampcombi complete` submodule.

Optional arguments:

| command | definition | default | example |
| ------------- | ------------- | ------------- | ------------- |
| --cluster_cov_mode | This assigns the cov. mode to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs [here](https://mmseqs.com/latest/userguide.pdf). | 0 | 2 |
| --cluster_mode | This assigns the cluster mode to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs [here](https://mmseqs.com/latest/userguide.pdf). | 1 | 2 |
| --cluster_coverage | This assigns the coverage to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs[here](https://mmseqs.com/latest/userguide.pdf).  | 0.8 | 0.9 |
| --cluster_seq_id | This assigns the seqsID to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs [here](https://mmseqs.com/latest/userguide.pdf).  | 0.4 | 0.7 |
| --cluster_sensitivity | This assigns sensitivity of alignment to the mmseqs2 cluster module- More information can be obtained in mmseqs2 docs [here](https://mmseqs.com/latest/userguide.pdf.) | 4.0 | 7.0 |
| --cluster_remove_singletons | This removes any hits that did not form a cluster. | True | False |
| --cluster_retain_label | This removes any cluster that only has a certain label in the sample name. For example if you have sample labels with 'S1_metaspades' and 'S1_megahit', you can retain clusters that have samples with suffix '_megahit' by running '--retain_clusters_label megahit'. | '' | 'megahit' |
| --cluster_min_member | This removes any cluster that has a hit number lower than assigned here. | 3 | 1 |

Output:

The output will be written into your working directory, containing the following files and folders:


```console
<pwd>/
└── Ampcombi_summary_cluster.tsv
├── Ampcombi_summary_cluster_representative_seq.tsv
└── Ampcombi_cluster.log
```

- **Ampcombi_summary_cluster.tsv** includes the contents of the complete summary plus a column with cluster IDs. 
- **Ampcombi_summary_cluster_representative_seq.tsv** contains a table with all the representative hits from each cluster.

---
`ampcombi signal_peptide`

This subcommand predicts whether a signal peptide was found on the filtered and clustered AMP hits. This only works if the user installs SignalP separately. For licensing issues, SignalP can only be downloaded and used by academic users; other users are requested to contact DTU Health Technology Software Package before using it. Please refer to [SignalP documentation](https://services.healthtech.dtu.dk/services/SignalP-6.0/).

To get a full list of options available and their defaults please refer to the help documentation of the submodule

```
ampcombi signal_peptide --help
```

Example usage: 
```console
ampcombi signal_peptide \
--signalp_model path/to/signalp_model/ \
--ampcombi_cluster path/to/Ampcombi_summary_cluster.tsv \
--log true
```
The input option `--ampcombi_cluster` requires the 'Ampcombi_summary_cluster.tsv' which is generated by the `ampcombi cluster` submodule. 

Output:

The output will be written into your working directory, containing the following files and folders:

```console
<pwd>/
└── Ampcombi_summary_cluster_SP.tsv
├── Ampcombi_summary_cluster_SP_onlyclusterswithSP.tsv
├── Ampcombi_summary_cluster_SP_onlyclusterswithSP.tsv
├── signalp
|   ├── output_*.png/
|   ├── prediction_results_index.tsv
|   ├── prediction_results.tsv
|   ├── representative_seq.txt
└── Ampcombi_signalpeptide.log
```

- **Ampcombi_summary_cluster_SP.tsv** includes the contents of the cluster summary plus a column with yes/no indicating the presence of a signal peptide sequence. 
- **Ampcombi_summary_cluster_SP_onlyclusterswithSP.tsv** includes the contents of the cluster summary plus a column with yes/no indicating the presence of a signal peptide sequence. But in this case clusters are retained only if they contain a hit or more with a signaling peptide.
- **signalp** contains the results from the tool signalp in `*.png` format showing the location of the predicted signaling peptide and the `prediction_results.tsv` which contains a table with the location of the signaling peptide and the identity. The `prediction_results_index.tsv` contains a table that gives an index number to every hit found in `./AMPcombi_summary_ao_human_nonhuman_clusters_SP_onlyclusterswithSP.tsv`. This can be used to rename the [localcolabfold](https://github.com/YoshitakaMo/localcolabfold) results downstream in the workflow. 

---

## Example runs using test files:

To test the function and output for AMPcombi, we provide test files that can be downloaded as described below, distributed in the three packages named `test_faa`, `test_gbk` and `test_files`.

Step1: Download the test directories and unzip :
```
wget https://github.com/Darcy220606/AMPcombi/tree/main/test_faa.tar.gz 
wget https://github.com/Darcy220606/AMPcombi/tree/main/test_gbk.tar.gz 
wget https://github.com/Darcy220606/AMPcombi/tree/main/test_files.tar.gz

tar -xzvf test_faa.tar.gz
tar -xzvf test_gbk.tar.gz
tar -xzvf test_files.tar.gz
```

Step2: Parse tables from all AMP tools. This can be produced by running the AMP workflow from [FUNCSCAN](https://github.com/nf-core/funcscan) - a pipeline for predicting functional genes in metagenomes.
```
ampcombi parse_tables --amp_results ./test_files/ --faa ./test_faa/ --gbk ./test_gbk/ --sample_list sample_1 sample_2 --ampir_file '.tsv' --amplify_file '.tsv' --macrel_file '.tsv' --neubi_file '.fasta' --hmmsearch_file '.txt' --ampgram_file '.tsv' --amptransformer_file '.txt' --threads 28 --log true
```

Step3: Concatenate all AMPcombi summary files from all samples:
```
ampcombi complete --summaries_directory ./test_ampcombi/ --log true
```

Step4: Cluster the hits and remove singletons:
```
ampcombi cluster --ampcombi_summary Ampcombi_summary.tsv --log true
```

---

## Visualization:

To visualize the result tables from AMPcombi, a pyshiny app can be rendered by running:

```
cd ./pyshiny
pip install -r requirements.txt
shiny run --port 36317 --reload app.py
```
 ⚠️ the port can be changed accordingly

The user can upload the `Ampcombi_summary_cluster_SP.tsv` to generate tables and figures ready for publication. 3D structures in PDB format can also be uploaded to generate an overlay structure.

---

## References for tools, packages and databases used in AMPcombi:

- Steinegger M and Soeding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature Biotechnology, doi: 10.1038/nbt.3988 (2017).

- Steinegger M and Soeding J. Clustering huge protein sequence sets in linear time. Nature Communications, doi: 10.1038/s41467-018-04964-5 (2018).

- Mirdita M, Steinegger M and Soeding J. MMseqs2 desktop and local web server app for fast, interactive sequence searches. Bioinformatics, doi: 10.1093/bioinformatics/bty1057 (2019).

- Mirdita M, Steinegger M, Breitwieser F, Soding J, Levy Karin E: Fast and sensitive taxonomic assignment to metagenomic contigs. Bioinformatics, doi: 10.1093/bioinformatics/btab184 (2021). 

- Teufel, F., Almagro Armenteros, J.J., Johansen, A.R. et al. SignalP 6.0 predicts all five types of signal peptides using protein language models. Nat Biotechnol 40, 1023–1025 doi: 10.1038/s41587-021-01156-3 (2022). 

- Buchfink B, Reuter K, Drost HG, Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature Methods 18, 366–368 doi:10.1038/s41592-021-01101-x  (2021).

- Shi G., Kang X., Dong F., Liu Y., Zhu N., Hu Y., Xu H., Lao X., Zheng H., DRAMP 3.0: an enhanced comprehensive data repository of antimicrobial peptides, Nucleic Acids Research, 50,D1, doi: 10.1093/nar/gkab651 (2022).

- The Shiny development team, Shiny for Python, https://shiny.posit.co/py/, license:https://github.com/posit-dev/py-shiny/blob/main/LICENSE v.0.8.0

- Inc., P. T. Collaborative data science. Montreal, QC: Plotly Technologies Inc. Retrieved from https://plot.ly (2015)

- Upsetplot, https://github.com/jnothman/UpSetPlot license: https://github.com/jnothman/UpSetPlot/blob/master/LICENSE v.0.9.0

- py3Dmol, https://github.com/avirshup/py3dmol license: https://github.com/avirshup/py3dmol/blob/master/LICENSE.txt

---

## Contribution:

AMPcombi is a tool developed for parsing results from published AMP prediction tools. We therefore welcome fellow contributors who would like to add new AMP prediction tools results for parsing and alignment. 

**Adding a new tool to AMPcombi:**

In `ampcombi/reformat_tables.py`
- add a new tool function to read the output to a pandas dataframe and return two columns named `contig_id` and `prob_<toolname>`
- add the new function to the `read_path` function

In `ampcombi/ampcombi.py`
- add a new parameter equivalent to the tool `<--tool_file>`.

---

**Authors and credits**: 

The tool was written mainly by @louperelo and @darcy220606 with major scientific contributions from @RosaLuzia.

---

**Funding**:

*This project was funded by Werner Siemens Foundation grant 'Palaeobiotechnology'*

---
