# Ampcombi
For testing out the ampcombi tool for running amp tools.

## Installation

## Usage

Note: AMP database sequences shoud not contain any character other than ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W',',Y']
Note: The database files should have an extension .fasta and .tsv

## Contribute
### Adding a new tool to AMPcombi
in `reformat_tables.py` 
- add the tools function to read the output to a pandas dataframe
- add the new function to the `read_path` function
in `main.py``
- add your default `tool:tool.fileending`to the default of `--tooldict`

# TODO : for PYPI format adapt to this style: 
===========
Towel Stuff
===========

Towel Stuff provides such and such and so and so. You might find
it most useful for tasks involving <x> and also <y>. Typical usage
often looks like this::

    #!/usr/bin/env python

    from towelstuff import location
    from towelstuff import utils

    if utils.has_towel():
        print "Your towel is located:", location.where_is_my_towel()

(Note the double-colon and 4-space indent formatting above.)

Paragraphs are separated by blank lines. *Italics*, **bold**,
and ``monospace`` look like this.


A Section
=========

Lists look like this:

* First

* Second. Can be multiple lines
  but must be indented properly.

A Sub-Section
-------------

Numbered lists look like you'd expect:

1. hi there

2. must be going

Urls are http://like.this and links can be
written `like this <http://www.example.com/foo/bar>`_.