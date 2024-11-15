.. _change:

Changelog
=========

v0.1.0
--------------
Initial release.

v0.1.1
--------------
- Minor changes.

v0.1.2
--------------
- Minor changes.

v0.1.3, 10.10.2022
------------------
- PyPi package and conda-recipe / biocontainer release.

v0.1.4, 18.10.2022
------------------
- Included a new optional argument "--complete_summary" to concatenate the results from multiple samples in one table.
- Added a universal log file to append to an existing log file rather than creating multiple new ones every time a sample is run.
- The "--path_list" can be called multiple times to include a list of files from individual samples in multiple lists.

v0.1.5, 27.10.2022
------------------
- Adapted reading of HMMER hmmsearch output to deal with varying header lines.
- Fixed syntax in "if" statements in "check_input.py".
- Included "check_faa_path" function to find .faa files also in subdirectories.

v0.1.6, 02.11.2022
------------------
- Included the HTML output for the complete summary.
- Added option --threads for DIAMOND (make database and alignment).
- Included check if database was downloaded once to prevent repeated downloads.

v0.1.7, 03.11.2022
------------------
- Added the option to submit a single .faa file instead of the faa-folder path (useful for summarizing a single sample).

v0.1.8, 14.02.2023
------------------
- Linked to Zenodo archive.

v0.2.0, 09.02.2024
------------------
- Changed the output extension from `.csv` to `.tsv`.
- Added a new feature to estimate the isoelectric point, molecular weight, structure fraction, and hydrophobicity.
- Filtered the DRAMP db `.tsv` to retain only necessary columns.
- Added the AMPtransformer and AMPgram tools.
- Enhanced hmmsearch to parse both single and multi models.
- Fixed dependencies in the `environment.yml` file.
- Created a log file for each sample, in addition to the main log for the complete summary (useful for pipelines like nf-core/funcscan).
- Added optional `--sample_metadata` and `--contig_metadata` flags.
- Removed `--tooldict` parameter; added individual parameters for each tool.
- Added parameters `--hmm_evalue`, `--aminoacid_length`, and `--db_evalue` to filter results based on specific criteria.
- Renamed `--cutoff` to `--amp_cutoff`.
- Included new parameters for input gbk files, window size for stop codons, transporter searches, and output handling.
- Added submodules for AMP clustering (using MMSeqs2) and signaling peptide prediction (using SignalP-6.0).
- Replaced HTML output with a Shiny for Python app accessible via the command line.
- Updated to subcommands for standardized use.

v0.2.1, 14.03.2024
------------------
- Fixed package versions in setup.
- Updated readme for installation setup.
- Fixed `./temp` dir removal at the end of the process.
- Changed the default matrix in DIAMOND blastp to PAM250.

v0.2.2, 21.03.2024
------------------
- Added a check for ./temp dir before attempting removal to prevent issues with pipelines.

v2.0.0, <XX.11.2024>
--------------------
- Added support for using InterProScan output (`--use_interproscan` and `--interproscan_output`) to remove AMPs classified as ribosomal proteins.
- Updated README and added documentation on Read the Docs.
- Improved DRAMP download script to remove non-amino acid characters.
- Added APD and UniRef100 databases to enhance AMP classification.
- Replaced DIAMOND alignment with MMseqs2 search.
- Fixed "NoneType" error when no AMP hits are retrieved or alignments fail to meet thresholds.
- Created Read the Docs for AMPcombi.