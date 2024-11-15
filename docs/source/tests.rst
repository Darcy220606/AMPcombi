.. _test:

Test runs
=========

InterProScan
------------

Example run for interproscan to generate input for ``--interproscan_output``.  
Please refer to the `InterProScan documentation <https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html>`_ for full description of how to use it.

Download the interproscan datasets as described `in the user documents <https://interproscan-docs.readthedocs.io/en/latest/UserDocs.html#obtaining-a-copy-of-interproscan>`_:

    .. code-block:: console

        singularity pull docker://interpro/interproscan:latest

        # Beaware! If the sequences contain any non-alphabet characters, it will crash!
        faa=path/to/faa/files 
        mkdir output temp

        for FN in $faa/*.faa; do
        N=$(basename $FN .faa) 
        singularity exec \
            -B $db:/opt/interproscan/data \
            -B $faa/output:/output \
            -B $faa/temp:/temp \
            -B $FN:$FN \
            interproscan_latest.sif \
            /opt/interproscan/interproscan.sh \
            --input $FN \
            --disable-precalc \
            --output-dir /output \
            --tempdir /temp \
            --cpu 16 \
            --applications PANTHER,ProSiteProfiles,ProSitePatterns,Pfam \
            --disable-residue-annot \
            --enable-tsv-residue-annot \
            --formats tsv  ; done

Example run
-----------

To test the functionality of AMPcombi, we provide test files for the required and optional inputs. 
Those can be found in the `tests directory <https://raw.githubusercontent.com/Darcy220606/AMPcombi/main/tests/>`_.

**Step1**: 
Download the test files and untar:

    .. code-block:: console

        git clone https://github.com/Darcy220606/AMPcombi.git

        tar -xzvf ./tests/test_faa.tar.gz
        tar -xzvf ./tests/test_gbk.tar.gz  
        tar -xzvf ./tests/test_files.tar.gz
        tar -xzvf ./tests/test_optional_files.tar.gz

📍 These input files can be generated ina  streamlined approach using `nf-core/funcscan <https://github.com/nf-core/funcscan>`_ - a pipeline for predicting functional genes in metagenomes.

**Step2**: 
Parse the tables and filter the AMP hits recovered:

    .. code-block:: console

        ampcombi parse_tables \
        --amp_results ./tests/test_files/ \
        --faa ./tests/test_faa/ \
        --gbk ./tests/test_gbk/ \
        --interproscan_output ./tests/test_optional_files/interproscan_output/ \
        --sample_metadata ./tests/test_optional_files/sample_metadata.tsv \
        --contig_metadata ./tests/test_optional_files/contig_metadata.tsv \
        --sample_list sample_1 sample_2 \
        --amp_database 'DRAMP' \
        --aminoacid_length 100 --db_evalue 100 --amp_cutoff 0.7 \
        --ampir_file '.tsv' --amplify_file '.tsv' --macrel_file '.tsv' --neubi_file '.fasta' --hmmsearch_file '.txt' --ampgram_file '.tsv' --amptransformer_file '.txt' \
        --log 'true' --threads 16

**Step3**: 
Concatenate the summary files:

    .. code-block:: console

        mv sample_1 test_ampsummaries/
        mv sample_2 test_ampsummaries/        
        
        ampcombi complete --summaries_directory ./test_ampsummaries --log 'true'

**Step4**: 
Cluster the the filtered AMPs into families:

    .. code-block:: console

        ampcombi cluster --ampcombi_summary Ampcombi_summary.tsv --log 'true' --threads

**Step5**: 
Predict signal peptides:

    .. code-block:: console

        ampcombi signal_peptide \
        --ampcombi_cluster Ampcombi_summary_cluster.tsv \
        --signalp_model ./signalpv6.0h-slowsequential/models --log 'true'


