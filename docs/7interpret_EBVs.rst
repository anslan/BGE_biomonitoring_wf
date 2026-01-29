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


Essential Biodiversity Variables of genetic composition
*******************************************************

Workflow to calculate `Essential Biodiversity Variables (EBV) <https://geobon.org/ebvs/what-are-ebvs/>`_ of genetic composition.

___________________________________________________

Starting point 
==============

Required input data: 

1. ASVs table

2. ASV taxonomy table
   
3. ASVs fasta file
   
4. OTU taxonomy table
   
5. Sample metadata


.. admonition:: example of an ASV table

  +-------------+---------+---------+---------+---------+-----+
  |             | sample1 | sample2 | sample3 | sample4 | ... |
  +-------------+---------+---------+---------+---------+-----+
  | **ASV_01**  | 579     | 0       | 0       | 0       | ... |
  +-------------+---------+---------+---------+---------+-----+
  | **ASV_02**  | 405     | 345     | 449     | 430     | ... |
  +-------------+---------+---------+---------+---------+-----+
  | **ASV_03**  | 0       | 0       | 231     | 69      | ... |
  +-------------+---------+---------+---------+---------+-----+
  | **ASV_04**  | 0       | 62      | 345     | 0       | ... |
  +-------------+---------+---------+---------+---------+-----+

.. admonition:: example of an ASV taxonomy table

  +------------+-----+------------+------------------+-----------------+---------------+
  |            | ... | Class      | Order            | Family          | Genus         |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **ASV_01** | ... | Collembola | Entomobryomorpha | Entomobryidae   | Entomobrya    |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **ASV_02** | ... | Collembola | Entomobryomorpha | **Family_0032** | **Genus_001** |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **ASV_03** | ... | Collembola | Entomobryomorpha | Isotomidae      | Parisotoma    |
  +------------+-----+------------+------------------+-----------------+---------------+
  | **ASV_04** | ... | Collembola | Entomobryomorpha | **Family_0032** | **Genus_022** |
  +------------+-----+------------+------------------+-----------------+---------------+


.. admonition:: example of a fasta file 

    .. code-block:: text
       :class: small-font

       >ASV_01
       ACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACTTCTCTTAGTTTATTAATTCG...
       >ASV_02
       ATAGTAGGAACATCTCTTAGTTTATTAATTCGAACTGAACTAGGAAATCCAGGTTCACTTATTGG...
       >ASV_03
       TCTTTACCTTTTATTCGGTGCCTGAGCTGGCATGGTGGGGACTGCTCTTAGTCTTCTAATCCGGG...
       >ASV_04
       TACTTTGTATTTTGTTTTTGGGGTGTGATCTGGTATGTTGGGGACTAGGTTCAGAAGACTAATTC...
       ...

.. admonition:: example of an OTU taxonomy table

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


.. admonition:: sample metadata

  +-------------+---------+-------------------+
  |             | site    | population        |
  +-------------+---------+-------------------+
  | **sample1** | site_01 | population_s01_01 |
  +-------------+---------+-------------------+
  | **sample2** | site_01 | population_s01_02 |
  +-------------+---------+-------------------+
  | **sample3** | site_02 | population_s02_05 |
  +-------------+---------+-------------------+
  | **sample4** | site_02 | population_s02_01 |
  +-------------+---------+-------------------+

___________________________________________________

Calculate EBVs of genetic composition: genetic richness, nucleotide diversity and genetic differentiation
=========================================================================================================


.. code-block:: R
    :caption: Preparing the input
    :linenos:

    #!/usr/bin/Rscript

    ## Script to calculate Essential Biodiversity Variables (EBV) of genetic composition
    ##  - haplotype richness
    ##  - nucleotide diversity
    ##  - genetic differentiation between populations

    ##-----------------------Settings-----------------------##
    ## Working directory and results directory paths
    wd = "E:/BGE/EBV_genetic composition"
    results = file.path(wd, "results")

    ## Input ASV table
    ASV.table = file.path(wd, "data/ASV_table_Arthropoda.csv")

    ## ASV taxonomy table. Information of ASV clustering into OTUs (by unique taxonomy, output of OptimOTU)
    ASV.tax = file.path(wd, "data/asv_taxonomy.tsv")

    ## OTU taxonomy table (all OTUs, the filtered and discarded sets, output of OptimOTU)
    all.OTU = file.path(wd, "data/otu_taxonomy.tsv")

    ## Sample metadata
    sampleMetadata = file.path(wd, "data/Sample_metadata.csv") #latest shared sampling metadata with "site_name", where we can see "A" or "1" for the main trap and "B" or "2" for the back up trap.

    ## ASV fasta file (ASV sequences)
    fasta_asv = file.path(wd, "data/ASVs_Arthropoda.fasta")

    ## Taxonomic filtering: specify the taxonomic level and target group
    ## Taxonomic level options: "phylum", "class", "order", "family", "genus", "species"
    target_taxonomic_level = "phylum"
    ## Target value for the specified taxonomic level
    target_taxonomic_group = "Arthropoda"
    ##------------------------------------------------------##



    library(data.table)
    library(stringr)
    library(iNEXT)
    library(ggplot2)
    library(pegas)
    library(ape)
    library(adegenet)
    library(mmod) #Jost's D & Gst



    ##------------------------##
    ## 1. Preparing the input ##
    ##------------------------##

    ## Load ASV and OTU data
    ##----------------------

    ## Load input ASV table:
    ## This table was obtained after a strict filtering with metaMATE 
    #                	    (to obtain reliable haplotype information)
    ASV_table = fread(file = ASV.table, header = TRUE, sep = ";")

    ## Load ASV taxonomy table (ASVs with taxonomic identification)
    ASVtax = fread(file = ASV.tax, header = TRUE)

    ## Load the OTU taxonomy table
    ## This list includes all OTUs detected, 
    #            including those discarded by the strict metaMATE filtering
    OTU_all = fread(file = all.OTU, header = TRUE, sep = "\t")

    ## ASVs presence checks between ASVtax and ASV_table
    ##--------------------------------------------------
    ## Check if all ASVs from ASV_table are present in ASVtax
    missing_in_ASVtax <- setdiff(ASV_table$ASV, ASVtax$seq_id)
    if(length(missing_in_ASVtax) > 0){
      warning(paste("The following", length(missing_in_ASVtax), 
                    "ASVs from ASV_table are not present in ASVtax:",
                    paste(head(missing_in_ASVtax, 10), collapse = ", "),
                    ifelse(length(missing_in_ASVtax) > 10, "...", "")))
    } else {
      message("All ASVs from ASV_table are present in ASVtax")
    }

    ## Check if all ASVs from ASVtax are present in ASV_table
    missing_in_ASV_table <- setdiff(ASVtax$seq_id, ASV_table$ASV)
    if(length(missing_in_ASV_table) > 0){
      warning(paste("The following", length(missing_in_ASV_table), 
                    "ASVs from ASVtax are not present in ASV_table:",
                    paste(head(missing_in_ASV_table, 10), collapse = ", "),
                    ifelse(length(missing_in_ASV_table) > 10, "...", "")))
    } else {
      message("All ASVs from ASVtax are present in ASV_table")
    }

    ## Find the intersection of ASVs (ASVs present in both matrices)
    common_ASVs <- intersect(ASV_table$ASV, ASVtax$seq_id)
    message(paste("Number of common ASVs in both matrices:", 
                  length(common_ASVs)))

    ## Filter both matrices to contain only matching ASVs
    ASVtax <- ASVtax[seq_id %in% common_ASVs]
    ASV_table <- ASV_table[ASV %in% common_ASVs]

    ## Add OTU id to ASVtax
    ##---------------------
    ASV.OTUtax <- merge(x=ASVtax, y=OTU_all[,.(OTU, species)], all.x = TRUE)

    ## Add OTU id to the ASV_table
    names(ASV.OTUtax)[which(names(ASV.OTUtax)=="seq_id")] <- "ASV"
    ASV_table = merge(x = ASV_table, y = ASV.OTUtax[, .(ASV, OTU)], 
                      by = "ASV", all.x = TRUE)

    ## Keeping only the specified taxonomic group (OTU and ASV)
    ##---------------------------------------------------------
    # Check if the specified taxonomic level exists in the data
    if(!target_taxonomic_level %in% names(ASV.OTUtax)){
      stop(paste("Error: Taxonomic level '", target_taxonomic_level, 
                "' not found in ASV.OTUtax. Available levels:",
                paste(names(ASV.OTUtax), collapse = ", ")))
    }

    ### Filter by the specified taxonomic level and group
    targetTax <- ASV.OTUtax[get(target_taxonomic_level) == 
                              target_taxonomic_group]
    if(nrow(targetTax) == 0){
      warning(paste("No ASVs found with", target_taxonomic_level, "=", 
                    target_taxonomic_group))
    } else {
      message(paste("Filtered to", nrow(targetTax), "ASVs with", 
                    target_taxonomic_level, "=", target_taxonomic_group))
    }

    ASV_table_targetTax <- ASV_table[ASV %in% targetTax$ASV]

    ## Convert ASV table to long format
    ASV = melt(data = ASV_table_targetTax,
              id.vars = colnames(ASV_table_targetTax)
              [c(1, ncol(ASV_table_targetTax))], 
              variable.name = "sample", value.name = "Abundance")

    # Remove ASV-sample combinations with zero abundance
    ASV = ASV[Abundance > 0]
    gc()


    ## Sample metadata - loading file
    ##-------------------------------
    sample_metadata = fread(file = sampleMetadata, header = TRUE, sep = ",")

    ## Adding "site & country" to ASV dataset
    ## Here, countries are used as 'populations' to compute Jost's D and Gst
    #          (for pairwise genetic distance comparisons between countries)
    # (Keeping the site column is optional)
    ASV_site <- merge(x=ASV, y=sample_metadata[, .(site, country, sample)], 
                      by="sample", all.x = TRUE)



OTU selection
-------------

We use `iNEXT R package <https://cran.r-project.org/web/packages/iNEXT/index.html>`__ to assess sample completeness (coverage) and expected diversity per OTU.
All OTUs with a coverage ≥ 0.95 and an expected/observed diversity ratio ≥ 0.8 were kept for calculating EBV of genetic composition.  


.. code-block:: R
    :caption: OTU selection
    :linenos:

    #!/usr/bin/Rscript

    ## Selecting OTUs (species):
    ##--------------------------
    ## Here, we work only with OTU that were well sampled. 
    # To identify them, we use iNEXT

    ## 1) Creating the input for iNEXT, it can be "incidence_raw",
    ## which is a presence-absence matrix, 
    ## or it can be "incidence_freq" which is a vector of frequencies,
    ## where the first position should be the number of sites 
    ## where the OTU was found.

    ## "incidence_raw"
    ## Here we also calculate the nucleotide diversity, using pegas::nuc.div()
    ## it calculates the average proportion of nucleotide differences 
    ## across all pairwise sequence comparisons.
    ## (sum(proportion of sites that differ between two sequences = 
    ##    no.differences/total length))/no.comparisons(=no.seqs*(no.seqs-1)/2)
    asv.seqs.all <- ape::read.dna(fasta_asv, format = "fasta")

    OTU_hap <- list()
    nucleotide.div <- list()

    ## OTUs in our sample:
    OTU2keep <- ASV_table_targetTax[, unique(OTU)]

    # Build a presence/absence matrix for each OTU and compute nucleotide diversity
    #  (this may take a minute or two)
    for(i in OTU2keep){
      asv_otu <- ASV_site[OTU==i, .(sample, ASV)]
      asv_otu[, presence := 1] # new variable fill of 1 that indicates presence
      ## To wide format: rows = ASV, column = samples, 0=absence, 1=presence
      asv_table <- dcast(asv_otu, ASV~sample, value.var = "presence")
      asv_table[is.na(asv_table)] <- 0
      namesASV <- asv_table[, ASV]
      asv_table[, ASV := NULL]
      asv_matrix <- as.matrix(asv_table)
      rownames(asv_matrix) <- namesASV
      OTU_hap[[i]] <- asv_matrix
      
      ## Nucleotide richness
      asv.seqs.subset <- asv.seqs.all[rownames(asv.seqs.all) %in% namesASV,]
      nucleotide.div[[i]] <- pegas::nuc.div(asv.seqs.subset)
    }

    ## "incidence_freq"
    hap_freq <- list()

    # Build a incidence frequency format for iNEXT rarefaction analysis
    #  (this may take a minute or two)
    for(i in OTU2keep){
      asv_otu <- ASV_site[OTU==i, .(sample, ASV)]
      asv_otu[, presence := 1] # new variable fill of 1 that indicates presence
      ## To wide format: rows = ASV, column = samples, 0=absence, 1=presence
      asv_table <- dcast(asv_otu, ASV~sample, value.var = "presence")
      asv_table[is.na(asv_table)] <- 0
      asv_table[, ASV := NULL]
      hap_freq[[i]] <- c(ncol(asv_table), apply(asv_table, 1, sum)) 
    }

    ## 2) Rarefaction using iNEXT
    ## --------------------------
    #out.raw <- iNEXT(hap_freq, q=0, datatype = "incidence_freq") #q=0 -> haplotype richness
    ## We get warnings and errors when trying to run all species at the same time.
    ## Instead, we do it one by one in a loop:

    results_list2 <- list()
    ## Set a progress bar to check the progress of the most demanding part of the loop
    pb <- txtProgressBar(min = 1, max = length(OTU2keep), style = 3)
    count=0

    for(sp in OTU2keep){
      count=count+1
      tryCatch({ #we added this because for some OTU is not possible to run iNEXT() and this way the loop does not stop executing
        results_list2[[sp]] <- iNEXT(hap_freq[sp], q=0, datatype = "incidence_freq")
      }, error = function(e){
        warning(paste("Problem with", sp))
      })
      
      
      setTxtProgressBar(pb, count) #progress bar
    }
    close(pb)

    results_list <- results_list2

    totry <- setdiff(OTU2keep, names(results_list)) #5770 OTUs did not run with iNEXT()

    ## However, some OTU that did not run with iNEXT() do run in a second trial (I don't know why)
    ##https://github.com/JohnsonHsieh/iNEXT/issues/81 (not solved)
    ##https://github.com/JohnsonHsieh/iNEXT/issues/75
    for(sp in totry){
      tryCatch({ #I added this because for some OTU is not possible to run iNEXT() and this way the loop does not stop executing
        results_list2[[sp]] <- iNEXT(hap_freq[sp], q=0, datatype = "incidence_freq")
      }, error = function(e){
        warning(paste("Problem with", sp))
      })
    }
    while(length(results_list2) > length(results_list)){
      print(paste("totry length", length(totry)))
      results_list <- results_list2
      for(sp in totry){
        tryCatch({ #I added this because for some OTU is not possible to run iNEXT() and this way the loop does not stop executing
          results_list2[[sp]] <- iNEXT(hap_freq[sp], q=0, datatype = "incidence_freq")
        }, error = function(e){
          warning(paste("Problem with", sp))
        })
      }
      totry <- setdiff(OTU2keep, names(results_list2))
    }

    ## Checking coverage (=completedness) and comparing observed and estimated number of 
    ##haplotypes per OTU
    results.data <- as.data.frame(matrix(NA, nrow=length(results_list), ncol = 2))
    names(results.data) <- c("coverage", "obs_est") #storing SC (sample coverage) and observed/estimated data
    rownames(results.data) <- names(results_list)

    for(i in 1:nrow(results.data)){
      otu = rownames(results.data)[i]
      results.data$coverage[i] <- results_list[[otu]]$DataInfo$SC
      results.data$obs_est[i] <- results_list[[otu]]$AsyEst$Observed[1]/results_list[[otu]]$AsyEst$Estimator[1] #Observed and estimated value of Species Richness
    }
    length(which(results.data$coverage == 1))/nrow(results.data)
    length(which(results.data$coverage > 0.97))/nrow(results.data)
    length(which(results.data$coverage > 0.95))/nrow(results.data)
    length(which(results.data$coverage > 0.90))/nrow(results.data)
    coverage_plot <- ggplot(results.data, aes(x=coverage)) +
      geom_histogram(binwidth = 0.01) +
      theme_minimal() +
      labs(title = "OTU frequency based on haplotype sampling completedness (coverage)", x = "coverage", y = "no. OTUs")
    coverage_plot

    length(which(results.data$obs_est == 1))/nrow(results.data)
    length(which(results.data$obs_est > 0.95))/nrow(results.data)
    length(which(results.data$obs_est > 0.90))/nrow(results.data)
    length(which(results.data$obs_est > 0.85))/nrow(results.data)
    length(which(results.data$obs_est > 0.80))/nrow(results.data)
    obsest_plot <- ggplot(results.data, aes(x=obs_est)) +
      geom_histogram(binwidth = 0.02) +
      theme_minimal() +
      labs(title = "OTU frequency based on no.observed haplotypes / estimated no.haplotypes", x = "obs/est", y = "no. OTUs")
    obsest_plot

    ## Selecting a threshold
    ## Here we select all OTUs with obs/est>=0.8 & coverage>=0.95
    OTUselected1 <- rownames(results.data)[intersect(which(results.data$coverage>=0.95), which(results.data$obs_est>=0.8))]



Haplotype richness and nucleotide diversity
-------------------------------------------

.. code-block:: R
    :caption: Haplotype richness and nucleotide diversity
    :linenos:

    #!/usr/bin/Rscript

    ##-----------------------------------------------------------------##
    ## 2. Calculating EBV: haplotype richness and nucleotide diversity ##
    ##-----------------------------------------------------------------##

    ## Adding to the OTU taxonomy table haplotype richness, nucleotide diversity per OTU,
    ##number of samples where the OTU was found and number of countries where each OTU was found:
    ebv_taxonomy <- OTU_all[OTU %in% OTUselected1][, `:=`(n.samples.final = 0, n.countries.final = 0, nucleotide.diversity = 0, haplotype.richness = 0)]

    for(otu in OTUselected1){
      samples.found.otu <- ASV[OTU==otu, sample]
      no.countries.found.otu <- ASV_site[sample %in% samples.found.otu, uniqueN(country)]
      ebv_taxonomy[OTU==otu, `:=`(haplotype.richness = length(hap_freq[[otu]])-1, 
                                  nucleotide.diversity = nucleotide.div[[otu]],
                                  n.samples.final = length(unique(samples.found.otu)),
                                  n.countries.final = no.countries.found.otu)]
    }



Genetic differentiation
-----------------------

Using the `BGE case study of high mountain systems <https://bioscanflow.readthedocs.io/en/latest/1_1HMS.html>`__ as an example, we calculate the genetic differentiation between countries.

In each country, we sampled five elevation points over a 20-week period, resulting in approximately 100 samples per country. 
For calculating Jost's D and Gst indices, each detection of a haplotype in a sample is treated as a single individual. For instance, if haplotype A is detected in samples 1 and 5 from Germany and in sample 29 from Spain, this corresponds to two individuals of haplotype A in Germany and one in Spain.


.. code-block:: R
    :caption: Genetic differentiation
    :linenos:

    #!/usr/bin/Rscript

    ##---------------------------------------------##
    ## 3. Calculating EBV: genetic differentiation ##
    ##---------------------------------------------##

    ## Here we calculate the genetic differentiation between countries
    ## 7 countries, one mountain gradient sampled per country
    ## The main factor grouping samples is "country".
    ## We calculate the overall genetic differentiation between countries per OTU
    ## and an average pairwise dissimilarity matrix across all OTUs using Jost's D.

    ## We work with all well sampled OTUs, removing those not shared between countries
    OTU2work <- ebv_taxonomy[n.countries.final>1, OTU]

    ## Countries are considered as populations
    pops <- ASV_site[OTU %in% OTU2work, unique(country)]

    ## Loop per OTU to calculate overall and pairwise genetic distance:
    D.matricesALL_list <- list()
    Gst.matricesALL_list <- list()
    ## New column in ebv_taxonomy to include the overall genetic differentiation
    ebv_taxonomy[, `:=`(JostD = 0, Gst =0)]
    ## Set a progress bar to check the progress of the most demanding part of the loop
    pb <- txtProgressBar(min = 1, max = length(OTU2work), style = 3)
    count=0
    for(otu in OTU2work){
      count=count+1
      ## (1) getting the ASVs of the OTU & adding individual ID
      asv_otu <- ASV_site[OTU==otu]
      asv_otu <- asv_otu[, indivID := paste0("indiv_", seq(1:nrow(asv_otu)))]
      ## (2) filtering the fasta file
      asv.seqs.subset <- asv.seqs.all[rownames(asv.seqs.all) %in% asv_otu$ASV,]
      otu.seqs.list <- as.list(asv.seqs.subset)
      ## (3) creating a sequence file following asv_otu (each indiv - one sequence --> duplicated asvs)
      indiv.seqs.list <- lapply(asv_otu$ASV, function(n) otu.seqs.list[[n]])
      names(indiv.seqs.list) <- paste0(asv_otu$indivID, "_", asv_otu$ASV)
      indiv.seqs <- do.call(rbind, indiv.seqs.list)
      class(indiv.seqs) <- "DNAbin"
      ## (3) alignment
      indiv.alignment <- muscle(indiv.seqs) #if not found: exec = "/usr/bin" (although this should work as it is in $PATH)
      ## (3) convert the alignment (DNAbin) in genind object (adegenet package)
      indiv.alignment.gd <- DNAbin2genind(indiv.alignment)
      ## (4) adding population & country information to the genind object 
      indiv.alignment.gd@pop <- as.factor(asv_otu$country)
      ##checking that individual order and population/country order is the same as in asv_otu
      for(check in 1:nrow(asv_otu)){
        stopifnot(str_detect(rownames(indiv.alignment.gd@tab)[check], asv_otu$indivID[check]))
        stopifnot(indiv.alignment.gd@pop[check] == asv_otu$country[check])
      }
      ## (5) Genetic differentiation
      ## (5.1) Overall Jost's D & Gst
      overall.D <- D_Jost(indiv.alignment.gd)
      overall.Gst <- Gst_Hedrick(indiv.alignment.gd)
      ebv_taxonomy[OTU==otu, `:=`(JostD = overall.D$global.het, Gst = overall.Gst$global)]
      ## (5.2) Pairwise Jost's D & Gst
      if(length(unique(asv_otu$country))==length(pops)){ #total number of populations/sites/elevation sampling points
        D.dist <- pairwise_D(indiv.alignment.gd)
        D.matrix <- as.matrix(D.dist)
        D.matricesALL_list[[otu]] <- D.matrix[pops, pops] #all matrices must have the same column and row order
        Gst.dist <- pairwise_Gst_Hedrick(indiv.alignment.gd)
        Gst.matrix <- as.matrix(Gst.dist)
        Gst.matricesALL_list[[otu]] <- Gst.matrix[pops, pops]
      }else{ #we need to build the whole popxpop matrix to later average all of them
        dist.obj.D <- pairwise_D(indiv.alignment.gd)
        dist.m.D <- as.matrix(dist.obj.D)
        dist.obj.Gst <- pairwise_Gst_Hedrick(indiv.alignment.gd)
        dist.m.Gst <- as.matrix(dist.obj.Gst)
        pop.matrix.D <- matrix(NA, nrow=length(pops), ncol=length(pops))
        colnames(pop.matrix.D) = rownames(pop.matrix.D) = pops
        pop.matrix.Gst <- matrix(NA, nrow=length(pops), ncol=length(pops))
        colnames(pop.matrix.Gst) = rownames(pop.matrix.Gst) = pops
        for(r in rownames(dist.m.D)){
          for(c in colnames(dist.m.D)){
            pop.matrix.D[which(rownames(pop.matrix.D)==r), which(colnames(pop.matrix.D)==c)] = dist.m.D[r,c]
            pop.matrix.Gst[which(rownames(pop.matrix.Gst)==r), which(colnames(pop.matrix.Gst)==c)] = dist.m.Gst[r,c]
          }
        }
        D.matricesALL_list[[otu]] <- pop.matrix.D
        Gst.matricesALL_list[[otu]] <- pop.matrix.Gst
      }
      setTxtProgressBar(pb, count) #progress bar
    }
    close(pb)

    ## Just in case, apply the same order in all matrices on the list (this would not be necessary as it is done inside the loop)
    D.matricesALL_list <- lapply(D.matricesALL_list, function(x) x[pops, pops])
    Gst.matricesALL_list <- lapply(Gst.matricesALL_list, function(x) x[pops, pops])

    ## Saving the two lists of pairwise distance matrices
    saveRDS(D.matricesALL_list, file.path(results, "OTUgeneticDisMatrices_JostD_ALLcountries.rds"))
    saveRDS(Gst.matricesALL_list, file.path(results, "OTUgeneticDisMatrices_Gst_ALLcountries.rds"))
    ## Calculate the mean across all matrices ignoring NAs
    D.mean_distanceALL <- Reduce("+", lapply(D.matricesALL_list, function(x){x[is.na(x)]<-0; x})) / Reduce("+", lapply(D.matricesALL_list, function(x) !is.na(x)))
    Gst.mean_distanceALL <- Reduce("+", lapply(Gst.matricesALL_list, function(x){x[is.na(x)]<-0; x})) / Reduce("+", lapply(Gst.matricesALL_list, function(x) !is.na(x)))
    ## Saving the two mean pairwise distance matrices between countries
    write.table(D.mean_distanceALL, file.path(results,"MeanPairwiseD_allCountries.txt"))
    write.table(Gst.mean_distanceALL, file.path(results,"MeanPairwiseGst_allCountries.txt"))


    ## Saving the final OTU table with the EBV of genetic composition:
    fwrite(ebv_taxonomy, file.path(results, "FinalOTUtable_geneticEBV.csv"))

   


____________________________________________________

|logo_BGE_small| |eufund| |chfund| |ukrifund|
