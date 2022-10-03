# Ampcombi
![Logo](https://github.com/Darcy220606/Ampcombi/tree/amp_fasta/docs/outputname.png)

This tool parses the results of amp prediction tools into a single table and aligns the hits against a reference database of antimicrobial peptides for functional classifications.

For parsing: AMpcombi is developed to parse the output of these **amp prediction tools**:
 
| Tool | Version | Link |
| ------------- | ------------- | ------------- |
| Ampir  | ?  | https://github.com/Legana/ampir |
| AMPlify  | ?  | https://github.com/bcgsc/AMPlify |
| Macrel  | ?  | https://github.com/BigDataBiology/macrel |
| HMMsearch  | ?  | ? |
| EnsembleAMPpred  | ? | ? |
| NeuBI  | ?  | https://github.com/nafizh/NeuBI |

For classification: AMPcombi is developed to offer functional annotation of the detcted AMPs by alignemnt to **AMP reference databases**, for e.g.,:

| Tool | Version | Link |
| ------------- | ------------- | ------------- |
| DRAMP  | 3.0 | https://github.com/CPU-DRAMP/DRAMP-3.0 |

======================
## Installation
======================

======================
## Usage:
======================

Note: AMP database sequences shoud not contain any character other than ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W',',Y']
Note: The database files should have an extension .fasta and .tsv

======================
## Contribution:
======================

As tool is for parsing results for published AMP prediction tools, we welcome fellow contributers who would like to add new AMP prediction tools results.

### Adding a new tool to AMPcombi
in `reformat_tables.py` 
- add the tools function to read the output to a pandas dataframe
- add the new function to the `read_path` function
in `main.py``
- add your default `tool:tool.fileending`to the default of `--tooldict`



======================

### TODO : for PYPI format adapt to this style: 
[python setup.py sdist](https://the-hitchhikers-guide-to-packaging.readthedocs.io/en/latest/creation.html)


=================================================================================================
**Authors**: @lperelo and @darcy220606


