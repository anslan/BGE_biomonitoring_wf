.. |eufund| image:: _static/eu_co-funded.png
  :width: 200
  :alt: Alternative text

.. |chfund| image:: _static/ch-logo-200x50.png
  :width: 210
  :alt: Alternative text

.. |ukrifund| image:: _static/ukri-logo-200x59.png
  :width: 150
  :alt: Alternative text

.. |logo_BGE_small| image:: _static/logo_BGE_alpha.png
  :width: 120
  :alt: Alternative text
  :target: https://biodiversitygenomics.eu/

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


VAL_OTU_ID: eDNA taxonomy validation
************************************

VAL_OTU_ID is a post-assignment validation pipeline for eDNA metabarcoding taxonomy.
It does not assign taxonomy itself, but checks, validates, and corrects taxonomic 
assignments based on geographic context and reference database completeness.

VAL_OTU_ID assesses the geographic plausibility of species-level assignments using 
occurrence records from global biodiversity databases (GBIF and OBIS), 
re-evaluates sequence matches against congeners in the reference database 
to account for incomplete or biased reference coverage, and considers whether 
alternative, locally occurring taxa represent more plausible identifications. 
The pipeline flags assignments that may warrant closer scrutiny 
by cross-referencing validated and original identifications against invasive 
alien species databases and the IUCN Red List. This framework reduces 
implausible taxonomic assignments and improves the ecological 
interpretability of metabarcoding results.

The full tutorial  
is hosted externally. You can either open it in a new tab or view it embedded below. 
(the embedded view loads the live page).

**Open in new tab:** `VAL_OTU_ID tutorial <https://gdunshea.github.io/VAL_OTU_ID/>`_

**Embedded tutorial** (content is loaded from the URL above; scroll inside the frame):

.. raw:: html

   <iframe src="https://gdunshea.github.io/VAL_OTU_ID/#" 
   title="VAL_OTU_ID: eDNA Taxonomy Validation Pipeline" 
   width="100%" height="950" style="border: 1px solid #333; 
   border-radius: 8px;"></iframe>


.. hide: 

  The code is hosted on GitHub `here <https://github.com/gdunshea/VAL_OTU_ID>`_.

____________________________________________________

|logo_BGE_small| |eufund| |chfund| |ukrifund|
