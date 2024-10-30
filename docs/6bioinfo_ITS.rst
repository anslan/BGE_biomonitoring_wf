.. |logo_BGE_alpha| image:: _static/logo_BGE_alpha.png
  :width: 400
  :alt: Alternative text
  :target: https://biodiversitygenomics.eu/

.. |eufund| image:: _static/eu_co-funded.png
  :width: 220
  :alt: Alternative text

.. |chfund| image:: _static/ch-logo-200x50.png
  :width: 210
  :alt: Alternative text

.. |ukrifund| image:: _static/ukri-logo-200x59.png
  :width: 150
  :alt: Alternative text

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


|logo_BGE_alpha|

Fungi/ITS
*********

| This is executable step-by-step pipeline for **rRNA ITS1/ITS2** amplicon data analyses from Illumina sequencing machine.
|  
| The **full bioinformatics workflow can be automatically run through** `PipeCraft2 <https://pipecraft2-manual.readthedocs.io/en/latest/>`_ (v1.1.0; releasing this soon, with a tutorial),
| which implemets also various error handling processes and sequence summary statistics (lacking here in step-by-step code). 
| 
| The bioinformatic workflow results in amplicon sequence variants (ASVs) and well as operational taxonomic units (OTUs).

+----------------------------------------------+--------------+---------+
| Process                                      | Software     | Version |
+==============================================+==============+=========+
| :ref:`Remove primers <remove_primers>`       | cutadapt     | 4.4     |
+----------------------------------------------+--------------+---------+
| :ref:`Quality filtering <quality_filtering>` | DADA2        | 2.28    |
+----------------------------------------------+--------------+---------+
| :ref:`Denoise <denoise>`                     | DADA2        | 2.28    |
+----------------------------------------------+--------------+---------+
| :ref:`Merge paired-end reads <denoise>`      | DADA2        | 2.28    |
+----------------------------------------------+--------------+---------+
| :ref:`Chimera filtering <remove_chimeras>`   | DADA2        | 2.28    |
+----------------------------------------------+--------------+---------+
| :ref:`Remove tag-jumps <tagjumps>`           | UNCROSS2     | x       |
+----------------------------------------------+--------------+---------+
| :ref:`Extract ITS1/2 region <itsx>`          | ITSx         | 1.1.3   |
+----------------------------------------------+--------------+---------+
| :ref:`Merge sequencing runs* <mergeRuns>`    |              |         |
+----------------------------------------------+--------------+---------+
| :ref:`Taxonomy assignment <taxAssign>`       | PROTAX/BLAST | x       |
+----------------------------------------------+--------------+---------+
| :ref:`Clustering ASVs to OTUs <clustering>`  | optimOTU     | x       |
+----------------------------------------------+--------------+---------+

\*only applicable when there are multiple sequencing runs per study. 

The bioinformatic workflow results in amplicon sequence variants (ASVs) and well as 
operational taxonomic units (OTUs). 

.. _remove_primers:

Remove primers
~~~~~~~~~~~~~~

.. _quality_filtering:

Quality filtering 
~~~~~~~~~~~~~~~~~


.. _denoise:

Denoise and merge paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _remove_chimeras:

Chimera filtering 
~~~~~~~~~~~~~~~~~

.. _tagjumps:

Remove tag-jumps
~~~~~~~~~~~~~~~~

.. _mergeRuns:

Merge sequencing runs
~~~~~~~~~~~~~~~~~~~~~

.. _itsx:

Extract ITS1/ITS2 region 
~~~~~~~~~~~~~~~~~~~~~~~~

Extract the ITS1/ITS2 region, i.e., clip conservative primer binding sites (18S, 5.8S, 28S) from the ASVs, for
making ASVs that differ only withing the ITS1/ITS2 part. 

.. code-block:: bash
   :caption: Extract ITS with ITSx
   :linenos:

    #output dir
    output_dir=$"ITSx_out"
    mkdir $output_dir
    mkdir -p tempdir2

    #load variables
    organisms=$"-t all"
    regions=$"--save_regions ITS2"
    partial=$"--partial 50"
    cores=$"--cpu 20"
    eval=$"-E 0.01"
    score=$"-S 0"
    domains=$"-N 2" 
    complement=$"--complement F"
    only_full=$"--only_full F"
    truncate=$"--truncate T"

    $itsx_path -i tempdir/$input.unique.$extension \
    -o tempdir/$input. \
    --preserve T \
    --graphical F \
    $organisms \
    $partial \
    $regions \
    $cores \
    $eval \
    $score \
    $domains \
    $complement \
    $only_full \
    $truncate




.. _taxAssign:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

.. _clustering:

Clustering ASVs to OTUs
~~~~~~~~~~~~~~~~~~~~~~~

A

____________________________________________________

|eufund| |chfund| |ukrifund|
