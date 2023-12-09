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


Arthropods/COI
**************

The below **full bioinformatics workflow can be run through XXX**, 
but see also subsections below for step-by-setp scripts.

+-------------------------------------------------+--------------+---------+
| Process                                         | Software     | Version |
+=================================================+==============+=========+
| :ref:`Remove primers <remove_primersCOI>`       | cutadapt     | 4.4     |
+-------------------------------------------------+--------------+---------+
| :ref:`Quality filtering <quality_filteringCOI>` | DADA2        | 2.28    |
+-------------------------------------------------+--------------+---------+
| :ref:`Denoise <denoiseCOI>`                     | DADA2        | 2.28    |
+-------------------------------------------------+--------------+---------+
| :ref:`Merge paired-end reads <denoiseCOI>`      | DADA2        | 2.28    |
+-------------------------------------------------+--------------+---------+
| :ref:`Remove chimeras <remove_chimerasCOI>`     | DADA2        | 2.28    |
+-------------------------------------------------+--------------+---------+
| :ref:`Remove tag-jumps <tagjumpsCOI>`           | UNCROSS2     | x       |
+-------------------------------------------------+--------------+---------+
| :ref:`Merge sequencing runs* <mergeRunsCOI>`    |              |         |
+-------------------------------------------------+--------------+---------+
| :ref:`Remove NUMTs <numtsCOI>`                  | metaMATE     | x       |
+-------------------------------------------------+--------------+---------+
| :ref:`Taxonomy assignment <taxAssignCOI>`       | PROTAX/BLAST | x       |
+-------------------------------------------------+--------------+---------+
| :ref:`Remove non-Metazoa <nonMetazoaCOI>`       | bash         | x       |
+-------------------------------------------------+--------------+---------+
| :ref:`Clustering <clusteringCOI>`               | optimOTU     | x       |
+-------------------------------------------------+--------------+---------+

\*only applicable when there are multiople sequencing runs per study. 

The bioinformatic workflow results in amplicon sequence variants (ASVs) and well as 
operational taxonomic units (OTUs). 

.. _remove_primersCOI:

Remove primers
~~~~~~~~~~~~~~

.. _quality_filteringCOI:

Quality filtering 
~~~~~~~~~~~~~~~~~


.. _denoiseCOI:

Denoise and merge paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. _remove_chimerasCOI:

Remove chimeras 
~~~~~~~~~~~~~~~

.. _tagjumpsCOI:

Remove tag-jumps
~~~~~~~~~~~~~~~~

.. _mergeRunsCOI:

Merge sequencing runs
~~~~~~~~~~~~~~~~~~~~~

.. _numtsCOI:

Remove NUMTs
~~~~~~~~~~~~


.. _taxAssignCOI:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

.. _nonMetazoaCOI:

Remove non-Metazoa
~~~~~~~~~~~~~~~~~~

.. _clusteringCOI:

Clustering
~~~~~~~~~~

A

____________________________________________________

|eufund| |chfund| |ukrifund|
