AMPcombi documentation
======================
.. raw:: html

   <img src="https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/docs/amp-combi-logo.png" alt="ampcombi icon" width="300" height="250">

----

.. image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :alt: license icon
.. image:: https://anaconda.org/bioconda/ampcombi/badges/downloads.svg
   :alt: downloads icon
.. image:: https://anaconda.org/bioconda/ampcombi/badges/version.svg
   :alt: version icon
.. image:: https://img.shields.io/pypi/dm/ampcombi
   :alt: Downloads

----

**AMPcombi** and its submodules provide a command-line interface to parse the results of antimicrobial peptide (AMP) prediction tools into a single table, aligns the AMP hits against reference AMP databases for functional and structural classifications, filters the AMP hits according to their physiochemical properties, clusters the filtered hits and predicts signaling peptides if present.

The input data for AMPcombi can be generated by running the individual AMP prediction tools described in (:ref:`usage`). AMPcombi functionality is also integrated in the AMP workflow of `nf-core/funcscan <https://github.com/nf-core/funcscan>`_. However not all AMP prediction tools implemented in AMPcombi can be found in nf-core/funcscan.

AMPcombi source code is hosted on `Github repository <https://github.com/Darcy220606/AMPcombi>`_

If you plan to use AMPcombi for your research please don't forget to **cite** us 💙


.. toctree::
   :maxdepth: 2
   :caption: Table of Contents:
   :glob:

   about
   usage
   tests
   visualization
   references
   contributing
   changes
