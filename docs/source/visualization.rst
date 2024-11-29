.. _visualization:


Visualization
=============

Dashboard
---------

To explore the final summary tables obtained from the different submodules and generate publication ready figures,
a user interface can be accessed :

**Option1 : Singularity image**

    .. code-block:: bash
        
        wget https://github.com/Darcy220606/AMPcombi-interface/releases/download/v2.0.0/ampcombi_interface.sif

        singularity run ampcombi_interface.sif

**Option2 : Through CLI for backend custom editing**

    .. code-block:: bash

        git clone https://github.com/Darcy220606/AMPcombi.git
        cd AMPcombi

        conda create -n ampcombi_gui python=3.13 -y
        conda activate ampcombi_gui
        pip install -r ./pyshiny/requirements.txt

        python -m shiny run --port 37231 --reload --autoreload-port 36257 ./pyshiny/app.py

    💡 After rendering the app, feel free to upload the test file in the interface. The test file can be found in ``./pyshiny/tests/Ampcombi_summery_cluster.tsv``.

    .. warning::

        - The updated app files can be found [here](https://github.com/Darcy220606/AMPcombi-interface).

        - This interface was created with an assumption that all AMPcombi submodules are run, including clustering of AMPs and prediction of signaling peptide. Additionally, for the taxonomy tab, it assumes that the user provided a column ``mmseqs_lineage_contig`` which contains the lineage format obtained from running `MMseqs2 taxonomy module <https://mmseqs.com/latest/userguide.pdf>`_ ,which can also be generated when running `nf-core/funcscan <https://github.com/nf-core/funcscan>`_ with the AMP and taxonomy workflows.

Example table: 

.. image:: https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/docs/ampcombi_interface_screenshot.png
   :alt: ampcombi screenshot
   :width: 800
   :height: 370

Example upset plot:

.. image:: https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/docs/ampcombi_interface_screenshot2.png
   :alt: ampcombi screenshot 2
   :width: 800
   :height: 370
   
----

