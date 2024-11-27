.. |eufund| image:: _static/eu_co-funded.png
  :width: 220
  :alt: Alternative text

.. |chfund| image:: _static/ch-logo-200x50.png
  :width: 210
  :alt: Alternative text

.. |ukrifund| image:: _static/ukri-logo-200x59.png
  :width: 150
  :alt: Alternative text

.. |hmsc1| image:: _static/hmsc1.png
  :width: 650
  :alt: Alternative text

.. raw:: html

    <style> .red {color:#ff0000; font-weight:bold; font-size:16px} </style>

.. role:: red


Hierarchical Modelling of Species Communities
*********************************************

Hierarchical Modelling of Species Communities (HMSC) is a framework for Joint Species Distribution Modelling; 
a model-based approach for analyzing community 
ecological data (`Ovaskainen et al.2017 <https://doi.org/10.1111/ele.12757>`_).


The **obligatory data** for HMSC-analyses includes a matrix of **species occurrences or abundances** and a 
matrix of **environmental covariates** (sampling units as **rows**). Optionally additional input data include species 
**traits** and **phylogeny**, and information about the spatiotemporal context of the 
sampling design. HMSC yields inference both at species and community levels. 

|hmsc1|


Starting point 
~~~~~~~~~~~~~~

The input for HMSC analysis consists of species/OTU community matrix (Y matrix) accompanied with 
the matrix of environmental covariates (X matrix), and optionally species 
traits (T matrix) and phylogeny (C matrix).

.. admonition:: note that samples are rows in the input Y and X matrices0

  *example of input Y (community) matrix*:

  +-------------+--------+--------+--------+--------+-----+
  |             | OTU_01 | OTU_02 | OTU_03 | OTU_04 | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample1** | 579    | 405    | 0      | 0      | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample2** | 0      | 345    | 0      | 62     | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample3** | 0      | 449    | 231    | 345    | ... |
  +-------------+--------+--------+--------+--------+-----+
  | **sample4** | 0      | 430    | 69     | 0      | ... |
  +-------------+--------+--------+--------+--------+-----+

  *example of input X (environmental metadata) matrix*:

  +-------------+---------+-----------------+----------+-----------+-----+
  |             | site    | collection_date | latitude | longitude | etc |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample1** | site_01 | 2022_07_15      | 0.1      | 0.2       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample2** | site_01 | 2022_07_21      | 0.1      | 0.2       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample3** | site_02 | 2022_07_28      | 1.0      | 1.0       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+
  | **sample4** | site_02 | 2022_08_05      | 1.0      | 1.0       | ... |
  +-------------+---------+-----------------+----------+-----------+-----+


.. admonition:: Traits

  Including the traits matrix helps to understand how species-specific traits influence the community assebly processes. 
  Traits matrix may include data for example about body size, shape, feeding type, etc. 

  +-------------+-----------+----------+---------------+
  |             | body_size | shape    | trophic_guild |
  +-------------+-----------+----------+---------------+
  | **OTU_01**  | 2         | square   | 1             |
  +-------------+-----------+----------+---------------+
  | **OTU_02**  | 2         | round    | 1             |
  +-------------+-----------+----------+---------------+
  | **OTU_03**  | 0.1       | round    | 2             |
  +-------------+-----------+----------+---------------+
  | **OTU_04**  | 0.2       | variable | 3             |
  +-------------+-----------+----------+---------------+


.. admonition:: Phylogeny

  Including the phylogeny matrix helps to understand if the species responses to the environmental covariates are phylogenetically structured, i.e., do similar species  respond similarly.

  The phylogeny matrix may be presented as Newick tree. 
  Alternatively, data on taxonomic idenity may used as a proxy of phylogenetic relatedness. 
  In the latter case, the UNCLASSIFIED taxonomic ranks (that are common in the metabarcoding data) should be informative 
  in a sense that **not all e.g. family level unclassified OTUs are closely related**. 
  For example, the distance between unclassified OTUs could be calculated and then OTUs falling within the user defined distance threshold could be classified as various levels of *pseudotaxa*. 

  *Example of the taxonomy table:*

  +------------+-----+------------+------------------+-----------------+---------------+
  |            | ... | Class      | Order            | Family          | Genus         |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_01** | ... | Collembola | Entomobryomorpha | Entomobryidae   | Entomobrya    |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_02** | ... | Collembola | Entomobryomorpha | **Family_0032** | **Genus_001** |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_03** | ... | Collembola | Entomobryomorpha | Isotomidae      | Parisotoma    |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **OTU_04** | ... | Collembola | Entomobryomorpha | **Family_0032** | **Genus_022** |
  +------------+-----+------------+------------------+-----------------+---------------+

  Here, for example the sequences of **OTU_02 and OTU_04 differ 9%**; so we consider those are originating from different genera but likely from the same family. 

___________________________________________________

Install HMSC
~~~~~~~~~~~~

.. code-block:: R
   :caption: install hmsc R package  
   :linenos:

   #!/usr/bin/Rscript

   # install 'devtools'; if not yet installed
   install.packages("devtools") 
   library(devtools)
   
   # install hmsc package
   install_github("hmsc-r/HMSC")

   # load hmsc
   library(Hmsc)

   # check the version
   packageVersion("Hmsc")

___________________________________________________

Define model
~~~~~~~~~~~~

Here, we are defining the model for HMSC. 
For example, our input Y and X matrixes are tab delimited text files where 
samples are in rows; and taxonomy table format follows the above example. 

.. code-block:: R
   :caption: load data and define model
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)

    # load community matrix (ASV/OTU/species table); tab-delimited txt file
    Y = read.table("OTU_table.txt", sep = "\t", 
                      check.names = F, header = T, row.names = 1)

    # load environmental covariates; tab-delimited txt file
    X = read.table("env_meta.txt", sep = "\t", 
                        check.names = F, header = T, row.names = 1)

    # load taxonomy; tab-delimited txt file; 
      # so we can use it as a proxy for phylogeny 
    taxonomy = read.table("taxonomy.txt", sep = "\t", 
                          check.names = F, header = T, row.names = 1)

    
    ### reducing a dataset by sub-selection OTUs that are present at least 99 samples
    prev = colSums(Y)
    sel.sp = prev>=100
    Y = Y[, sel.sp]
    taxonomy = taxonomy[sel.sp, ]

    # creating a phylogenetic tree from a set of nested taxonomic variables in the taxonomy table
    taxonomy$Phylum = as.factor(taxonomy$Phylum)
    taxonomy$Class = as.factor(taxonomy$Class)
    taxonomy$Order = as.factor(taxonomy$Order)
    taxonomy$Family = as.factor(taxonomy$Family)
    taxonomy$Genus = as.factor(taxonomy$Genus)
    tax.tree = as.phylo(~Phylum/Class/Order/Family/Genus, 
                      data = taxonomy, collapse = FALSE)



    # defining our study design; structure of the data
    studyDesign = data.frame(
                    sample = as.factor(rownames(X)), 
                    site = as.factor(X$site)
                    )

    # convert sample collection dates into Julian days relative to a specific start date 
    da = as.Date(X$collection_date)
    jday =  1 + julian(da) - julian(as.Date("2022-07-15"))

    XData = data.frame(seqdepth = log(read_counts$raw_nread), jday)


    #we use 3.141593 instead of pi because otherwise R will think that pi is a variable in the model and script S7_make_predictions will yield an error message
    XFormula = ~cos(2*3.141593*jday/365) +  sin(2*3.141593*jday/365) + cos(4*3.141593*jday/365) + sin(4*3.141593*jday/365) + seqdepth

    for(i in 4:12) taxonomy[,i] = as.factor(taxonomy[,i])
    phy.tree = as.phylo.formula(~kingdom/phylum/class/order/family/subfamily/tribe/genus/species,data=taxonomy)
    phy.tree$tip.label = colnames(Y)
    plot(phy.tree,cex=0.2)
    rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
    rL.site = HmscRandomLevel(units = levels(studyDesign$site))

    m = Hmsc(Y = Y, distr = "probit",
              XData = XData, XFormula = XFormula,
              phyloTree = phy.tree,
              studyDesign = studyDesign,
              ranLevels = list(sample = rL.sample,site = rL.site))

    models = list(m)
    names(models) = c("Swedish_seasonal_model")
    save(models, file = paste0("models/unfitted_models.RData"))
    models

 
Fit model 
~~~~~~~~~

The following HMSC pipeline is a modified version of the pipeline presented at the `ISEC 2024 Hmsc workshop <https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc>`_.




____________________________________________________

|eufund| |chfund| |ukrifund|
