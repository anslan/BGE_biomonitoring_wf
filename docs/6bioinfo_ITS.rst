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

The below **full bioinformatics workflow can be run through XXX**, 
but see also subsections below for step-by-setp scripts.

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
| :ref:`Extract ITS2 region <itsx>`            | ITSx         | 1.1.3   |
+----------------------------------------------+--------------+---------+
| :ref:`Merge sequencing runs* <mergeRuns>`    |              |         |
+----------------------------------------------+--------------+---------+
| :ref:`Taxonomy assignment <taxAssign>`       | PROTAX/BLAST | x       |
+----------------------------------------------+--------------+---------+
| :ref:`Clustering ASVs to OTUs <clustering>`  | optimOTU     | x       |
+----------------------------------------------+--------------+---------+

\*only applicable when there are multiople sequencing runs per study. 

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

Extract ITS2 region 
~~~~~~~~~~~~~~~~~~~


.. _taxAssign:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

.. _clustering:

Clustering ASVs to OTUs
~~~~~~~~~~~~~~~~~~~~~~~

A

____________________________________________________

|eufund| |chfund| |ukrifund|
