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

.. |output_icon| image:: _static/output_icon.png
  :width: 50
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
traits (T) matrix and phylogeny (C) matrix.

.. admonition:: Community and Environment (metadata)

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
  |             | site    | collection_date | latitude | longitude | ... |
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

  The phylogeny matrix may be presented as **Newick tree** file. 
  Alternatively, data on **taxonomic ranks** may used as a proxy of phylogenetic relatedness. 
  In the latter case, the UNCLASSIFIED taxonomic ranks *(that are common in the metabarcoding data)* 
  should be informative   in a sense that 
  **not all e.g. family level unclassified OTUs would be considered as closely related** because 
  they have the 'same label'.   For example, the distance between unclassified OTUs could 
  be calculated and then OTUs falling within the user defined distance threshold could be 
  classified as various levels of *pseudotaxa*. 

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

  Here, for example the sequences of **OTU_02 and OTU_04 differ 9.8%**; so we consider those are originating from different genera but likely from the same family. 

___________________________________________________

Install HMSC
~~~~~~~~~~~~

Install hmsc R package (if already not installed).

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

Select data
~~~~~~~~~~~

In this 'select data' section, we are assuming that 
our input Y and X matrixes are tab delimited text files where 
samples are in rows; and taxonomy table format follows the above example.

Here, we are also deciding if we want to proceed with the **presence-absence or abundance** (read count)
community matrix.

.. note:: 

  Before using to the full dataset, try fitting the model with a **small subset**
  for faster model testing and validation.


.. code-block:: R
   :caption: select data and prep. data
   :linenos:

   #!/usr/bin/Rscript

   # input matrices file names; according to the above HMSC figure.
      # here, expecting all files to be tab-delimited.
   Y_file = "OTU_table.txt" # Community; samples are rows
   X_file = "env_meta.txt"  # Environment; samples are rows
   C_file = "taxonomy.txt"  # Phylogeny (optional) 
                            # (herein a 'pseudo-phylogeny' based on the 
                            # assigned taxonomic ranks; species are rows)
   T_file = ""              # Traits (optional)
   sp_prevalence = 100      # set a species prevalence threshold; 
                            # meaning that perform HMSC for species 
                            # in a Y-matrix, that occur >= 100 samples in this case. 
    #-------------------------------------------------------------------------------------------#

    # load community matrix Y
    Y = read.table(Y_file, sep = "\t", 
                      check.names = F, 
                      header = T, 
                      row.names = 1)
    # load metadata matrix X
    X = read.table(X_file, sep = "\t", 
                        check.names = F, 
                        header = T, 
                        row.names = 1)
    # load taxonomy matrix C; 
      # so we can use it as a proxy for phylogeny 
    taxonomy = read.table(C_file, sep = "\t", 
                          check.names = F, 
                          header = T, 
                          row.names = 1)

    # herein, converting Y matrix to presence-absence
    Y = 1*(Y>0)

    # select OTUs/species that are present at least in $sp_prevalence (specified above) samples
    prevalence = colSums(Y != 0)
    select.sp = prevalence >= sp_prevalence
    Y = Y[, select.sp]
    taxonomy = taxonomy[select.sp, ]

    ### creating a phylogenetic tree from a set of nested taxonomic ranks in the taxonomy table
     # ranks as.factors; assuming that ranks start from the 2nd column in the taxonomy table 
        # and we have 7 ranks
    for(i in 2:8) taxonomy[,i] = as.factor(taxonomy[,i])
     # convert tax to phylo tree
    phy.tree = as.phylo(~Phylum/Class/Order/Family/Genus/Species, 
                    data = taxonomy, collapse = FALSE)
     
    # this "pseudo-tree" does not have any branch lengths, which are needed for the model;
    # assign arbitrary branch lengths
      if (is.null(phy.tree$edge.length)) {
        phy.tree$edge.length = rep(1, nrow(phy.tree$edge))
      }
     # rename tree tip lables according to the labels in the Y matrix.
     phy.tree$tip.label = colnames(Y)

     # check the tree 
     plot(phy.tree, cex=0.6)

___________________________________________________

Define model
~~~~~~~~~~~~

According to our dataset, we are defining the model for HMSC.
That is, we specify the structure, including the response variable (community data), 
covariates (environmental predictors), random effects, and phylogenetic relationships. 

.. code-block:: R
   :caption: define model
   :linenos:

   #!/usr/bin/Rscript

   library(Hmsc)

    # defining our study design; structure of the data
    studyDesign = data.frame(
                    sample = as.factor(rownames(X)),    # rownames(X) = sample names
                    site = as.factor(X$Site)
                    )

    ### incorporating random effects into the  HMSC model. 
     # (to capture the influence of unmeasured factors that vary across different 
     # levels of the data; e.g., among sites, samples). 
    # sampling units
    rL.sample = HmscRandomLevel(units = levels(studyDesign$sample))
    # sampling sites
    rL.site = HmscRandomLevel(units = levels(studyDesign$site)) 

    # convert sample collection dates into Julian days relative to a specific start date 
    da = as.Date(meta$CollectionDate)
    jday = 1 + julian(da) - julian(as.Date("2024-01-01"))

    # not covered here, but DOWNLOAD RELEVANT COVARIATES (E.G., CLIMATE, WEATHER, LANDCOVER) 
    # FROM DATABASES BASED ON COORDINATES AND SAMPLING TIMES

    # create a data frame of covariates (predictor variables) that will be included in the model
    XData = data.frame(seqdepth = log(meta$seq_depth),  # number of sequences per sample
                        jday)                           # Julian days

    ### specify the formula for the fixed effects
    # using 3.141593 instead of pi to prevent issues when 'pi' is considered as a variable; 
    XFormula = ~cos(2*3.141593*jday/365) +  # model seasonal effects; annual cycles
                sin(2*3.141593*jday/365) +  # model seasonal effects; annual cycles
                cos(4*3.141593*jday/365) +  # model seasonal effects; semiannual cycles
                sin(4*3.141593*jday/365) +  # model seasonal effects; semiannual cycles
                seqdepth                    # number of sequences per sample

    ### define a model 
    m = Hmsc(Y = Y,             # response matrix
          distr = "probit",      # distribution model for the response variable ('probit' for PA)
          XData = XData,          # predictor variables 
          XFormula = XFormula,     # fixed effects in the model
          phyloTree = phy.tree,     # phylogenetic tree object
          studyDesign = studyDesign, # study design object
          ranLevels = list(sample = rL.sample, site = rL.site)) # random level objects
    
    # organize, name, and save your HMSC models (to easily manage multiple models if needed)
    models = list(m)
    names(models) = c("model_1")
    save(models, file = paste0("models/unfitted_models.RData"))

    # check models
    models

___________________________________________________

Fit model
~~~~~~~~~

The following HMSC pipeline is a modified version of the pipeline presented at the `ISEC 2024 Hmsc workshop <https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/hmsc>`_.

| ``input:`` unfitted_models file (**unfitted_models.RData**)
| ``output:`` fitted models, with fitting done for multiple RUNs:

...


____________________________________________________

|eufund| |chfund| |ukrifund|
