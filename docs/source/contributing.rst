.. _contributing:

Contributing
============

Help always wanted
------------------

AMPcombi needs your help.
AMPcombi is a tool developed for parsing results from published AMP prediction tools. 
We therefore welcome fellow contributors who would like to add new AMP prediction tools results for parsing and alignment. 
Whether it's contributing patches, submitting bug reports, or just asking and
answering questions, you are helping make AMPcombi a better tool.

Adding a new tool
-----------------

We first ask you to first make an issue on `AMPcombi issues <https://github.com/Darcy220606/AMPcombi/issues>`_

Then in ``ampcombi/reformat_tables.py``:

        - add a new tool function to read the output to a pandas dataframe and return two columns labelled ``contig_id`` and ``prob_toolname``
        - add the new function to the ``read_path`` function

In ``ampcombi/ampcombi.py``
- add a new parameter equivalent to the new tool ``--tool_file``.
