.. |logo_BGE_alpha| image:: _static/logo_BGE_alpha.png
  :width: 300
  :target: https://biodiversitygenomics.eu/

.. |plutoFgo_phone| image:: _static/plutof/plutoFgo_phone.png
  :width: 200

.. |upload| image:: _static/plutof/upload.png
  :width: 400

.. |new_project| image:: _static/plutof/new_project.png
  :width: 500

.. |eufund| image:: _static/eu_co-funded.png
  :width: 200

.. |chfund| image:: _static/ch-logo-200x50.png
  :width: 210

.. |ukrifund| image:: _static/ukri-logo-200x59.png
  :width: 150

.. |logo_BGE_small| image:: _static/logo_BGE_alpha.png
  :width: 120
  :target: https://biodiversitygenomics.eu/

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red

.. |to_settings| image:: _static/plutof/to_settings.png
  :width: 200
  :align: top

.. |settings| image:: _static/plutof/settings.png
  :width: 200
  :align: top

.. |enable_material_sample| image:: _static/plutof/enable_material_sample.png
  :width: 200
  :align: top

.. |arrow| raw:: html

   <span style="font-size: 2em; color: green; vertical-align: middle;">&#8594;</span>
  
.. |arrow_down| raw:: html

   <span style="font-size: 4em; color: green;">&#8595;</span>

.. |material_sample_tab| image:: _static/plutof/material_sample_tab.png
  :width: 250


|logo_BGE_alpha|

.. _registering_samples_in_plutof:

Registering samples in PlutoF and PlutoF GO 
*******************************************

The `PlutoF Biodiversity Data Management Platform <https://plutof.ut.ee/en>`_ 
(Abarenkov et al., `2010 <https://doi.org/10.4137/EBO.S6271>`_) is an open access data management platform for biodiversity data 
including observations, specimens, material samples, sequences and related metadata. 

PlutoF supports streamlined, automated submission of curated sample and 
occurrence data to major international repositories such as `ENA (European Nucleotide Archive) <https://www.ebi.ac.uk/ena/browser/home>`_ 
and `GBIF (Global Biodiversity Information Facility) <https://www.gbif.org/>`_, 
helping to ensure that datasets are **findable, accessible, interoperable and 
reusable (FAIR)**. 
PlutoF centralises data storage, enforces consistent data structures and 
vocabularies, facilitates collaboration among project partners, 
and reduces the risk of data loss or duplication compared to maintaining 
separate local spreadsheets or databases.


`PlutoF GO <https://plutof.ut.ee/go>`_ is data collection tool 
for biodiversity data in a handy **phone application** that 
can be used to record samples during fieldwork, and upload 
them to the PlutoF platform.

|plutoFgo_phone|


PlutoF **video tutorials** available at: https://www.youtube.com/@PlutoFplatform/videos 

___________________________________________________

Creating a project in PlutoF
----------------------------

Before using the PlutoF GO app, `become a user <https://app.plutof.ut.ee/register>`_ in 
the `PlutoF website <https://plutof.ut.ee/en>`_.

**Registering samples in PlutoF GO requires a Project ID**. 
Generate one if needed in `PlutoF <https://plutof.ut.ee/en>`_ 
(Main menu -> Projects -> New).

|new_project|

Projects allow to group and organise data for 
bulk data management, set up pre-defined sampling areas, manage 
user access rights and roles within a project, and prepare data 
for publishing. Projects can be public, private, or shared with collaborators.

____________________________________________________

Add samples
-----------

| **This is an example how to add 'material sample' to PlutoF data management platform.**
|     *Material samples are for example soil, water, malaise trap samples.*


.. admonition:: Open PlutoF GO

  **1.** Open the PlutoF GO application on your phone/tablet, 
     and enable material sample gathering through 'Settings'. 

|to_settings| |arrow| |settings| |arrow| |enable_material_sample| 

``Add material sample`` is now displayed on the main screen.           

|material_sample_tab|

.. admonition:: Add material sample
  
  | **2.** Go to ``Add material sample``. 
  | GPS coordinates are captured automatically (but can be edited in the ``Location`` box)

.. admonition:: Choose Project (mandatory)
  
  | **3.** Choose Project by typing in your Project ID.
  | Project ID where the recorded samples are **allocated by default** can be added in via 'Settings'.
  | GPS coordinates are captured automatically (but can be edited in the ``Location`` box)

.. admonition:: Add Sample ID (mandatory)
  
  | **4.** Add Sample ID or scan the `QR-code of a pre-registered sample <https://www.youtube.com/watch?v=1My4Vn10YkA>`_ *(click to open the link to tutorial video "QR Codes for Project" in youtube)*.

.. figure:: _static/malaise_bottle.jpg
   :width: 250
   :align: center

   *Example of a QR code sticker on a sample, ready for scanning with the PlutoF GO app.*

.. admonition:: Fill other optional fields
  
  | **5.** Fill other optional fields and ``save`` the record. 

Note that **images, videos, audio and other sample associated files** can be also added.

.. image:: _static/plutof/multimedia.jpg
  :width: 500
  :align: center

.. admonition:: Upload 
  
  | **6.** Upload samples to PlutoF
  | After pushing ``save`` the records are only **locally saved**; and can be edited. 
  | **Press the cloud icon** on the top right corner to export records to PlutoF platform. 

| |upload|


| Once the records are in PlutoF, they can be further edited only in `PlutoF web platform <https://plutof.ut.ee/en>`_.

___________________________________________________

.. admonition:: Note

  Examples of "material sample" data entry using PlutoF 
  GO for soil, and malaise trap samples are included also in the SOPs in  
  the :ref:`BGE case studies section <casestudies>`.

____________________________________________________

Managing sample data on the PlutoF workbench
--------------------------------------------

Once the data has been uploaded from the 
PlutoF GO app (or imported via CSV files), it can be validated, 
organised, and further edited on the `PlutoF workbench <https://plutof.ut.ee/en>`_.
All records can be updated individually or in bulk (via functionality in Clipboard).

.. figure:: _static/plutof/project_view.png
  :width: 600
  :align: center

  *Project view; red circle indicate (bottom left) links to project-related search and clipboard functionalities*

___________________________________________________

Uploading sample records to ENA
-------------------------------

`ENA (European Nucleotide Archive) <https://www.ebi.ac.uk/ena/browser/home>`_ 
is an internationally recognized public repository for nucleotide 
sequence data and associated sample metadata, 
ensuring your data is **findable, accessible, and reusable**.

While PlutoF facilitates biodiversity data management, 
the raw metabarcoding sequence data cannot be deposited to PlutoF. 
Therefore, PlutoF enables automated submission of sample records to ENA,
where samples **can be linked to raw sequence data**. 


Using the **Publishing module** in PlutoF, 
users can submit sample records 
to the ENA database. 
The PlutoF platform acts as a broker for ENA, 
utilising its programmatic Webin submission 
service for sample data submission. 
The resulting ENA and BioSamples identifiers 
are stored in PlutoF alongside the original sample records.

Publishing is project-based: all samples within a 
selected study are submitted together, and the dataset 
can be updated later by re-publishing.

.. figure:: _static/plutof/to_ENA.png
  :width: 650
  :align: center

*To publish your dataset in ENA, go to Main menu -> Laboratories -> Publishing Lab -> ENA Datasets -> New*


**Steps to publish:**

| **1.** Select the project name from the autocomplete list (project **moderator rights** are required).
| **2.** Fill in any missing mandatory field values as required by ENA.
| **3.** Save the dataset.
| **4.** After administrator approval, publish the dataset to ENA.

.. note::
  
  Samples uploaded to ENA are treated as independent samples - 
  that is, they are not linked to a BioProject. 
  **Samples will be linked to a BioProject when the raw sequence data is associated 
  with the samples.**
  See 'How to Submit Raw Reads' `tutorial here <https://ena-docs.readthedocs.io/en/latest/submit/reads.html>`_.



___________________________________________________

|logo_BGE_small| |eufund| |chfund| |ukrifund|
