.. _visualization:


Visualization
=============

Dashboard
---------

To explore the final summary tables obtained from the different submodules and generate publication ready figures to explore the datasets, we set up a user interface. 
This can be accessed by running the following from the CLI: 

    .. code-block:: bash

        git clone https://github.com/Darcy220606/AMPcombi.git
        cd AMPcombi

        conda create -n ampcombi_gui python=3.13 -y
        conda activate ampcombi_gui
        pip install -r ./pyshiny/requirements.txt

        python -m shiny run --port 37231 --reload --autoreload-port 36257 ./pyshiny/app.py

    💡 After rendering the app, feel free to upload the test file in the interface. The test file can be found in ``./pyshiny/tests/Ampcombi_summery_cluster.tsv``.

    .. warning::

        This interface was created with an assumption that all AMPcombi submodules are run, including clustering of AMPs and prediction of signaling peptide.
        Additionally, for the taxonomy tab, it assumes that the user provided a column ``mmseqs_lineage_contig`` which contains the lineage format obtained from running 
        `MMseqs2 taxonomy module <https://mmseqs.com/latest/userguide.pdf>`_ ,which can also be generated when running `nf-core/funcscan <https://github.com/nf-core/funcscan>`_ 
        with the AMP and taxonomy workflows.


.. image:: https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/docs/ampcombi_interface_screenshot.png
   :alt: ampcombi screenshot
   :width: 400
   :height: 400

.. image:: https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/docs/ampcombi_interface_screenshot2.png
   :alt: ampcombi screenshot 2
   :width: 400
   :height: 400
   
----

