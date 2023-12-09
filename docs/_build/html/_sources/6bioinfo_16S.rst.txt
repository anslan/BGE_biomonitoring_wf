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


Bacteria/16S
************

The below **full bioinformatics workflow can be run through XXX** (coming soon), 
but see also subsections below for step-by-setp scripts.

+-------------------------------------------------+---------------------------+-------------------+
| Process                                         | Software                  | Version           |
+=================================================+===========================+===================+
| :ref:`Remove primers <remove_primers16S>`       | cutadapt                  | 4.4               |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Quality filtering <quality_filtering16S>` | DADA2                     | 2.28              |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Denoise <denoise16S>`                     | DADA2                     | 2.28              |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Merge paired-end reads <denoise16S>`      | DADA2                     | 2.28              |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Remove chimeras <remove_chimeras16S>`     | DADA2                     | 2.28              |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Remove tag-jumps <tagjumps16S>`           | UNCROSS2                  |                   |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Merge sequencing runs* <mergeRuns16S>`    |                           |                   |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Taxonomy assignment <taxAssign16S>`       | naive Bayesian classifier | as in DADA2 v2.28 |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Clustering <clustering16S>`               | vsearch                   | 2.26.0            |
+-------------------------------------------------+---------------------------+-------------------+

\*only applicable when there are multiople sequencing runs per study. 

The bioinformatic workflow results in amplicon sequence variants (ASVs) and well as 
operational taxonomic units (OTUs). 


Starting point
~~~~~~~~~~~~~~

.. important:: 

  When aiming to combine samples from multiple sequencing runs, then follow the below directory structure 

**Directory structure**

| **/my_working_folder**
| ├── **/sequencing_set01**
| │   ├── *sample1.R1.fastq*
| │   ├── *sample1.R2.fastq*
| │   ├── *sample2.R1.fastq*
| │   ├── *sample2.R2.fastq*
| │   ├── ...
| ├── **/sequencing_set02**
| │   ├── *sampleA.R1.fastq*
| │   ├── *sampleA.R2.fastq*
| │   ├── *sampleB.R1.fastq*
| │   ├── *sampleB.R2.fastq*
| │   ├── ...
| └── **/sequencing_set03**
|     ├── *sample11.R1.fastq*
|     ├── *sample11.R2.fastq*
|     ├── *sample12.R1.fastq*
|     ├── *sample12.R2.fastq*
|     ├── ...

.. note:: 
  
  Files with the **same name** will be considered as the same sample and will be merged in the "Merge sequencing runs" step.

.. _remove_primers16S:

Remove primers
~~~~~~~~~~~~~~

.. code-block:: bash
   :caption: using cutadapt to remove primers

   #Specify primers 
   fwd_primer=$"GTGYCAGCMGCCGCGGTAA"
   rev_primer=$"GGCCGYCAATTYMTTTRAGTTT"

   #remove primers across all files (samples)
   cutadapt \
        $mismatches \
        $min_length \
        $overlap \
        $indels \
        $cores \
        $untrimmed_output \
        $untrimmed_paired_output \
        --pair-filter=$pair_filter \
        -G $fwd_primer \
        -g $rev_primer \
        -o $output_dir/$inputR1.$extension \
        -p $output_dir/$inputR2.$extension \
        $inputR1.$extension $inputR2.$extension



.. _quality_filtering16S:

Quality filtering 
~~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: quality filtering in DADA2

   #in R, load dada2 library 
   library('dada2')
   
   #perfrom quality filtering
   qfilt = filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                    maxN = maxN, 
                    maxEE = c(maxEE, maxEE), 
                    truncQ = truncQ,  
                    truncLen = c(truncLen_R1, truncLen_R2),
                    maxLen = maxLen, 
                    minLen = minLen, 
                    minQ = minQ, 
                    rm.phix = TRUE, 
                    matchIDs = FALSE,
                    compress = FALSE, 
                    multithread = TRUE)


.. _denoise16S:

Denoise and merge paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: denoising in DADA2

    #Learn the error rates
    errF = learnErrors(filtFs, multithread = TRUE)
    errR = learnErrors(filtRs, multithread = TRUE)

    #dereplicate
    derepFs = derepFastq(filtFs, qualityType = qualityType)
    derepRs = derepFastq(filtRs, qualityType = qualityType)
    saveRDS(derepFs, (file.path(path_results, "derepFs.rds")))
    saveRDS(derepRs, (file.path(path_results, "derepRs.rds")))

    #denoise
    dadaFs = dada(derepFs, err = errF, pool = pool, selfConsist = selfConsist, multithread = TRUE)
    dadaRs = dada(derepRs, err = errR, pool = pool, selfConsist = selfConsist, multithread = TRUE)
    saveRDS(dadaFs, (file.path(path_results, "dadaFs.rds")))
    saveRDS(dadaRs, (file.path(path_results, "dadaRs.rds")))

    #merge paired-end reads
    merge = mergePairs(dadaFs, derepFs, dadaRs, derepRs, 
                            maxMismatch = maxMismatch,
                            minOverlap = minOverlap,
                            justConcatenate = justConcatenate,
                            trimOverhang = trimOverhang)

    #write RDS object
    saveRDS(ASV_tab, (file.path(path_results, "ASVs_table.denoised-merged.rds")))

    #seq count summary
    getN <- function(x) sum(getUniques(x))
        #remove 0 seqs samples from qfilt statistics
        row_sub = apply(qfilt, 1, function(row) all(row !=0 ))
        qfilt = qfilt[row_sub, ]

    seq_count <- cbind(qfilt, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merge, getN))
    colnames(seq_count) <- c("input", "qualFiltered", "denoised_R1", "denoised_R2", "merged")
    rownames(seq_count) <- sample_names
    write.table(seq_count, file.path(path_results, "seq_count_summary.txt"), 
                            sep = "\t", col.names = NA, 
                            row.names = TRUE, quote = FALSE)

.. _remove_chimeras16S:

Remove chimeras 
~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: remove chimeras in DADA2

    # Remove chimeras
    ASV_tab.nochim = removeBimeraDenovo(ASV_table, method="consensus", 
                                                  multithread=TRUE)
    saveRDS(ASV_tab.nochim, "ASV_table/ASV_tab.nochim.rds")

.. _tagjumps16S:

Remove tag-jumps
~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: removing putative tag-jumps with UNCROSS2 method
      
    set_f = 0.05
    set_p = 1

    library(data.table)

    ## Load OTU table
    cat("..Loading ASV table\n")

    if (file.exists("ASV_table/ASV_tab_len_filt.rds") == TRUE) {
        tab = readRDS("ASV_table/ASV_tab_len_filt.rds")
        cat("input table = ASV_table/ASV_tab_len_filt.rds\n")
    } else {
      tab = readRDS("ASV_table/raw_ASV_tab.rds")
      cat("input table = ASV_table/raw_ASV_tab.rds\n")
    }

    OTUTABW = as.data.table(t(tab), keep.rownames = TRUE)
      colnames(OTUTABW)[1] <- "ASV"
      cat("...Number of ASVs: ", nrow(OTUTABW), "\n")
      cat("...Number of samples: ", ncol(OTUTABW) - 1, "\n")

    ## Convert to long format
    cat("..Converting OTU table to long format\n")
    OTUTAB <- melt(data = OTUTABW, id.vars = "ASV",
      variable.name = "SampleID", value.name = "Abundance")

    ## Remove zero-OTUs
    OTUTAB <- OTUTAB[ Abundance > 0 ]
    cat("...Number of non-zero records: ", nrow(OTUTAB), "\n")

    ## Estimate total abundance of sequence per plate
    cat("..Estimating total OTU abundance\n")
    OTUTAB[ , Total := sum(Abundance, na.rm = TRUE), by = "ASV" ]

    ## UNCROSS score 
    ## (with original parameter - take a root from the exp in denominator, to make curves more steep)
    uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1){
      # x = OTU abundance in a sample
      # N = total OTU abundance
      # n = number of samples
      # f = expected cross-talk rate, e.g. 0.01
      # tmin = min score to be considered as cross-talk
      # p = power to rise the exponent (default, 1; use 1/2 or 1/3 to make cureves more stepp)

      z <- f * N / n               # Expected treshold
      sc <- 2 / (1 + exp(x/z)^p)   # t-score
      res <- data.table(Score = sc, TagJump = sc >= tmin)
      return(res)
    }

    ## Esimate UNCROSS score
    cat("..Estimating UNCROSS score\n")
    OTUTAB <- cbind(
      OTUTAB,
      uncross_score(
        x = OTUTAB$Abundance,
        N = OTUTAB$Total,
        n = length(unique(OTUTAB$SampleID)),
        f = as.numeric(set_f),
        p = as.numeric(set_p)
        )
      )

    ## Truncate singletons with total OTU abundance > 99 reads
    # OTUTAB[ Abundance == 1 & Total > 99  , TagJump := TRUE ]
    # OTUTAB[ Abundance == 2 & Total > 999 , TagJump := TRUE ]

    cat("...Number of tag-jumps: ", sum(OTUTAB$TagJump, na.rm = TRUE), "\n")

    ## Plot
    cat("..Making a plot\n")
    PP <- ggplot(data = OTUTAB, aes(x = Total, y = Abundance, color = TagJump)) +
        geom_point() + scale_x_log10() + scale_y_log10() +
        scale_color_manual(values = c("#0C7C59", "#D64933")) +
        labs(x = "Total abundance of ASV, reads", y = "Abundance of ASV in a sample, reads")

    cat("..Exporting a plot\n")
    pdf(file = "ASV_table/TagJump_plot.pdf", width = 12, height = 9.5, useDingbats = FALSE)
      PP
    dev.off()

    ## TJ stats
    cat("..Calculating tag-jump summary\n")
    TJ <- data.table(
        Total_reads = sum(OTUTAB$Abundance),
        Number_of_TagJump_Events = sum(OTUTAB$TagJump),
        TagJump_reads = sum(OTUTAB[ TagJump == TRUE ]$Abundance, na.rm = T)
        )

    TJ$ReadPercent_removed <- with(TJ, (TagJump_reads / Total_reads * 100))

    fwrite(x = TJ, file = "ASV_table/TagJump_stats.txt", sep = "\t")

    ## Exporting tag-jump data
    cat("..Exporting tag-jump data\n")
    JMPS <- OTUTAB[ TagJump == TRUE, .(ASV, SampleID) ]

    saveRDS(object = JMPS,
      file = "ASV_table/TagJump_OTUs.RData",
      compress = "xz")

    ## Prepare OTU tables, remove tag-jumps
    cat("..Removing tag-jumps\n")

    OTUTAB <- OTUTAB[ TagJump == FALSE ]

    ## Convert to wide format
    RES <- dcast(data = OTUTAB,
      formula = ASV ~ SampleID,
      value.var = "Abundance", fill = 0)

    ## Sort rows (by total abundance)
    clz <- colnames(RES)[-1]
    otu_sums <- rowSums(RES[, ..clz], na.rm = TRUE)
    RES <- RES[ order(otu_sums, decreasing = TRUE) ]

    cat("..Exporting tag-jump filtered table\n")

    ## Output table that is compadible with dada2
    output = as.matrix(RES, sep = "\t", header = TRUE, rownames = 1, 
                            check.names = FALSE, quote = FALSE)
    output = t(output)
    saveRDS(output, ("ASV_table/ASV_tab_TagJumpFiltered.rds"))

.. _mergeRuns16S:

Merge sequencing runs
~~~~~~~~~~~~~~~~~~~~~

.. _taxAssign16S:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

.. _clustering16S:

Clustering
~~~~~~~~~~

A

____________________________________________________

|eufund| |chfund| |ukrifund|
