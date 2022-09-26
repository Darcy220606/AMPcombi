# Ampcombi_test
For testing out the ampcombi tool for running amp tools.

## Instalation

## Usage

## Contribute
### Adding a new tool to AMPcombi
in `reformat_tables.py` 
- add the tools function to read the output to a pandas dataframe
- add the new function to the `read_path` function
in `main.py``
- add your default `tool:tool.fileending`to the default of `--tooldict`
