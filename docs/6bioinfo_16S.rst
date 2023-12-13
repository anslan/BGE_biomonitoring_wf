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

The below **full bioinformatics workflow be automatically run through PipeCraft2** (coming soon), 
but see also subsections below for step-by-setp scripts.
[workflow follows the implementations as in PipeCraft2 (https://github.com/pipecraft2/pipecraft), cite if using ...]

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
| :ref:`Chimera filtering <remove_chimeras16S>`   | DADA2                     | 2.28              |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Remove tag-jumps <tagjumps16S>`           | UNCROSS2                  |                   |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Merge sequencing runs* <mergeRuns16S>`    |                           |                   |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Taxonomy assignment <taxAssign16S>`       | naive Bayesian classifier | as in DADA2 v2.28 |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Clustering ASVs to OTUs <clustering16S>`  | vsearch                   | 2.26.0            |
+-------------------------------------------------+---------------------------+-------------------+
| :ref:`Post-clustering <postclustering16S>`      | LULU, BLAST               | 0.1.0; 2.15.0     |
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


| #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
| When working with a **single directory** that hosts your fastq 
| files, then use the code inside these lines [see below]
| #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#

.. _remove_primers16S:

Remove primers
~~~~~~~~~~~~~~

Remove primer strings from paired-end data. 

.. note:: 
  
  Here, assuming that all sequences are in 5'-3' orientation! 
  *(3'-5' orient sequences will be discarded with this workflow)*

.. important:: 

  | - Paired-end files must contain "R1" and "R2" strings (not just _1 and _2)!
  | - Sample names must not contain "R1" and "R2" strings (i.e. not FR123_001_R1.fastq/FR123_001_R2.fastq)

.. code-block:: bash
   :caption: remove primers with cutadapt

    #!/bin/bash
    ## workflow to remove primers via cutadapt [Sten Anslan]

    # Working folder = /my_working_folder
    # specify the identifier string for the R1 files
    read_R1=".R1"

    # specify primers 
    fwd_primer=$"GTGYCAGCMGCCGCGGTAA"    #this is primer 515F
    rev_primer=$"GGCCGYCAATTYMTTTRAGTTT" #this is primer 926R

    # get directory names if working with multiple sequencing runs
    DIRS=$(ls -d *) # -> sequencing_set01 sequencing_set02 sequencing_set03

    for sequencing_run in $DIRS; do 
        printf "\nWorking with $sequencing_run \n"
        cd $sequencing_run
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
        mkdir primersCut_out
        mkdir -p primersCut_out/untrimmed

        ### Clip primers with cutadapt
        for inputR1 in *$read_R1*; do
            inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
            cutadapt --quiet \
            -e 1 \
            --minimum-length 32 \
            --overlap 19 \
            --no-indels \
            --cores 4 \
            --untrimmed-output primersCut_out/untrimmed/$inputR1 \
            --untrimmed-paired-output primersCut_out/untrimmed/$inputR2 \
            --pair-filter=both \
            -g $fwd_primer \
            -G $rev_primer \
            -o primersCut_out/$inputR1 \
            -p primersCut_out/$inputR2 \
            $inputR1 $inputR2
        done
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
        cd ..
    done


.. _quality_filtering16S:

Quality filtering 
~~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: quality filtering in DADA2 (in R)

    #!/usr/bin/Rscript
    ## workflow to perform quality filtering within DADA2 [Sten Anslan]

    #load dada2 library 
    library('dada2')

    # the identifier string for the R1/R2 files
    read_R1 = ".R1"
    read_R2 = ".R2"
    # the delimiter for sample name (e.g. in that case the 
        #string before the first . is the sample name in a file sample1.R1.fastq)    samp_ID = "\\."

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/my_working_folder"
    dirs = list.dirs(recursive = FALSE)
    for (i in 1:length(dirs)) {
        if(length(dirs) > 1) {
            setwd(dirs[i])
            print(paste0("Working with ", dirs[i]))
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            # output path
            path_results = "qualFiltered_out"
            # input and output file paths
            R1s = sort(list.files("primersCut_out", pattern = read_R1, full.names = TRUE))
            R2s = sort(list.files("primersCut_out", pattern = read_R2, full.names = TRUE))
            #sample names
            sample_names = sapply(strsplit(basename(R1s), samp_ID), `[`, 1)
            #filtered files path
            filtR1 = file.path(path_results, paste0(sample_names, ".R1.", "fastq.gz"))
            filtR2 = file.path(path_results, paste0(sample_names, ".R2.", "fastq.gz"))
            names(filtR1) = sample_names
            names(filtR2) = sample_names
            
            #quality filtering
            qfilt = filterAndTrim(R1s, filtR1, R2s, filtR2, 
                                maxN = 0, 
                                maxEE = c(2, 2), 
                                truncQ = 2,  
                                truncLen = c(0, 0),
                                maxLen = 600, 
                                minLen = 100, 
                                minQ = 2, 
                                rm.phix = TRUE, 
                                matchIDs = TRUE,
                                compress = TRUE, 
                                multithread = TRUE)
            saveRDS(qfilt, file.path(path_results, "qfilt_reads.rds"))

            # make sequence count report
            seq_count = cbind(qfilt)
            colnames(seq_count) = c("input", "qualFiltered")
            rownames(seq_count) = sample_names
            write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), 
                                row.names = TRUE, quote = FALSE)

            #save filtered R objects for denoising and merging (below)
            filtR1 = sort(list.files(path_results, pattern = ".R1.fastq.gz", full.names = TRUE))
            filtR2 = sort(list.files(path_results, pattern = ".R2.fastq.gz", full.names = TRUE))
            sample_names = sapply(strsplit(basename(filtR1), ".R1.fastq.gz"), `[`, 1)
            saveRDS(filtR1, file.path(path_results, "filtR1.rds"))
            saveRDS(filtR2, file.path(path_results, "filtR2.rds"))
            saveRDS(sample_names, file.path(path_results, "sample_names.rds"))
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            #set working directory back to "/my_working_folder"
            setwd(wd)
        i = i + 1
        }
    }

.. _denoise16S:

Denoise and merge paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: denoise and merge paired-end reads in DADA2

    #!/usr/bin/Rscript
    ## workflow to perform DADA2 denoising and merging [Sten Anslan]

    #load dada2 library 
    library('dada2')

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/my_working_folder"
    dirs = list.dirs(recursive = FALSE)
    for (i in 1:length(dirs)) {
        if(length(dirs) > 1) {
            setwd(dirs[i])
            print(paste0("Working with ", dirs[i]))
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            #load quality filtered files
            filtR1 = readRDS("qualFiltered_out/filtR1.rds")
            filtR2 = readRDS("qualFiltered_out/filtR2.rds")
            qfilt = readRDS("qualFiltered_out/qfilt_reads.rds")
            sample_names = readRDS("qualFiltered_out/sample_names.rds")

            # create output dir
            path_results = "denoised_merged"
            dir.create(path_results, showWarnings = FALSE)

            print("# Denoising ...")
            # learn the error rates
            errF = learnErrors(filtR1, multithread = TRUE)
            errR = learnErrors(filtR2, multithread = TRUE)

            # make error rate figures
            pdf(file.path(path_results, "Error_rates_R1.pdf"))
              print( plotErrors(errF) )
            dev.off()
            pdf(file.path(path_results, "Error_rates_R2.pdf"))
              print( plotErrors(errR) )
            dev.off()

            # dereplicate
            derepR1 = derepFastq(filtR1, qualityType = "Auto")
            derepR2 = derepFastq(filtR2, qualityType = "Auto")

            # denoise
            dadaR1 = dada(derepR1, err = errF, 
                            pool = FALSE, selfConsist = FALSE, 
                            multithread = TRUE)
            dadaR2 = dada(derepR2, err = errR, 
                            pool = FALSE, selfConsist = FALSE, 
                            multithread = TRUE)

            # merge paired-end reads
            print("# Merging ...")
            merge = mergePairs(dadaR1, derepR1, dadaR2, derepR2, 
                                maxMismatch = 2,
                                minOverlap = 15,
                                justConcatenate = FALSE,
                                trimOverhang = FALSE)
            #make sequence table
            ASV_tab = makeSequenceTable(merge)
            rownames(ASV_tab) = gsub("R1.fastq.gz", "", rownames(ASV_tab))
            #write RDS object
            saveRDS(ASV_tab, (file.path(path_results, "rawASV_table.rds")))

            # make sequence count report
            getN <- function(x) sum(getUniques(x))
            #remove 0 seqs samples from qfilt statistics
            row_sub = apply(qfilt, 1, function(row) all(row !=0 ))
            qfilt = qfilt[row_sub, ]
            seq_count <- cbind(qfilt, sapply(dadaR1, getN), 
                                sapply(dadaR2, getN), sapply(merge, getN))
            colnames(seq_count) <- c("input", "qualFiltered", "denoised_R1", "denoised_R2", "merged")
            rownames(seq_count) <- sample_names
            write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), 
                                    row.names = TRUE, quote = FALSE)
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            print("--------")
            setwd(wd)
        i = i + 1
        }
    }



.. _remove_chimeras16S:

Chimera filtering 
~~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: remove chimeras in DADA2

    #!/usr/bin/Rscript
    ## workflow to perform chimera filtering within DADA2 [Sten Anslan]

    # load libraries
    library('dada2')
    library('openssl')

    # chimera filtering method
    method = "consensus"  #"consensus" vs. "pooled"
    # collapse ASVs that have no mismatshes or internal indels (identical up to shifts and/or length)
    collapseNoMismatch = "true"  #true/false 

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/my_working_folder"
    dirs = list.dirs(recursive = FALSE)
    for (i in 1:length(dirs)) {
        if(length(dirs) > 1) {
            setwd(dirs[i])
            print(paste0("Working with ", dirs[i]))
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            # load denoised and merged ASVs
            rawASV_table = readRDS("denoised_merged/rawASV_table.rds")
            # create output dir
            path_results="ASV_table"
            dir.create(path_results, showWarnings = FALSE)
            # Remove chimeras
            print("Removing chimeric ASVs ...")
            chim_filt = removeBimeraDenovo(
                                rawASV_table, method = method, 
                                multithread = TRUE,
                                verbose = TRUE)
            saveRDS(chim_filt, "ASV_table/chim_filt.rds")

            ### format and save ASV table and ASVs.fasta files
            # sequence headers
            asv_seqs = colnames(chim_filt)
            asv_headers = openssl::sha1(asv_seqs)
            # transpose sequence table
            tchim_filt = t(chim_filt)
            # add sequences to 1st column
            tchim_filt = cbind(row.names(tchim_filt), tchim_filt)
            colnames(tchim_filt)[1] = "Sequence"
            # row names as sequence headers
            row.names(tchim_filt) = asv_headers
            # write ASVs.fasta to path_results
            asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
            write(asv_fasta, file.path(path_results, "ASVs.fasta"))
            # write ASVs table to path_results
            write.table(tchim_filt, file.path(path_results, "ASV_table.txt"), 
                                    sep = "\t", col.names = NA, 
                                    row.names = TRUE, quote = FALSE)

            ### collapse ASVs that have no mismatshes or internal indels 
                                # (identical up to shifts and/or length)
            if (collapseNoMismatch == "true") {
                print("Collapsing identical ASVs ...")
                ASV_tab_collapsed = collapseNoMismatch(chim_filt, 
                                    minOverlap = 20, orderBy = "abundance", 
                                    identicalOnly = FALSE, vec = TRUE, 
                                    band = -1, verbose = TRUE)
                saveRDS(ASV_tab_collapsed, file.path(path_results, "ASV_tab_collapsed.rds"))

                ### format and save ASV table and ASVs.fasta files
                # sequence headers
                asv_seqs = colnames(ASV_tab_collapsed)
                asv_headers = openssl::sha1(asv_seqs)
                # transpose sequence table
                tASV_tab_collapsed = t(ASV_tab_collapsed)
                # add sequences to 1st column
                tASV_tab_collapsed = cbind(row.names(tASV_tab_collapsed), tASV_tab_collapsed)
                colnames(tASV_tab_collapsed)[1] = "Sequence"
                #row names as sequence headers
                row.names(tASV_tab_collapsed) = asv_headers
                # write ASVs.fasta to path_results
                asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
                write(asv_fasta, file.path(path_results, "ASVs_collapsed.fasta"))
                # write ASVs table to path_results
                write.table(tASV_tab_collapsed, file.path(path_results, "ASVs_table_collapsed.txt"), 
                                        sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

                # print summary
                print(paste0("Output = ", length(colnames(ASV_tab_collapsed)), 
                                " chimera filtered (+collapsed) ASVs, with a total of ", 
                                sum(rowSums(ASV_tab_collapsed)), 
                                " sequences."))
                print("--------")
            } else {
                # print summary
                print(paste0("Output = ", length(colnames(chim_filt)), 
                                " chimera filtered ASVs, with a total of ", 
                                sum(rowSums(chim_filt)), 
                                " sequences."))
                print("--------")
            }
                    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            setwd(wd)
        i = i + 1
        }
    }


.. _tagjumps16S:

Remove tag-jumps
~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: removing putative tag-jumps with UNCROSS2 method

   #!/usr/bin/Rscript
   ## Script to perform tag-jump removal; (C) Vladimir Mikryukov,
                                        # edit, Sten Anslan

    # load libraries
    library(data.table)

    # set parameters
    set_f = 0.03 # f-parameter of UNCROSS (e.g., 0.03)
    set_p = 1    # p-parameter (e.g., 1.0)

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/my_working_folder"
    dirs = list.dirs(recursive = FALSE)
    for (i in 1:length(dirs)) {
        if(length(dirs) > 1) {
            setwd(dirs[i])
            print(paste0("Working with ", dirs[i]))
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            # load ASV table
             # loading ASV_table_collapsed if collapseNoMismatch was "true" (above)
            if (file.exists("ASV_table/ASV_table_collapsed.rds") == TRUE) {
                tab = readRDS("ASV_table/ASV_table_collapsed.rds")
                cat("input table = ASV_table/ASV_table_collapsed.rds\n")
            } else {
              tab = readRDS("ASV_table/chim_filt.rds")
              cat("input table = ASV_table/chim_filt.rds\n")
            }

            ASVTABW = as.data.table(t(tab), keep.rownames = TRUE)
            colnames(ASVTABW)[1] <- "ASV"
            # convert to long format
            ASVTAB <- melt(data = ASVTABW, id.vars = "ASV",
            variable.name = "SampleID", value.name = "Abundance")
            ## Remove zero-OTUs
            ASVTAB <- ASVTAB[ Abundance > 0 ]
            # estimate total abundance of sequence per plate
            ASVTAB[ , Total := sum(Abundance, na.rm = TRUE), by = "ASV" ]

            ## UNCROSS score
            uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1){
              z <- f * N / n               # Expected treshold
              sc <- 2 / (1 + exp(x/z)^p)   # t-score
              res <- data.table(Score = sc, TagJump = sc >= tmin)
              return(res)
            }

            # esimate UNCROSS score
            cat(" estimating UNCROSS score\n")
            ASVTAB <- cbind(
              ASVTAB,
              uncross_score(
                x = ASVTAB$Abundance,
                N = ASVTAB$Total,
                n = length(unique(ASVTAB$SampleID)),
                f = as.numeric(set_f),
                p = as.numeric(set_p)
                )
              )

            cat(" number of tag-jumps: ", sum(ASVTAB$TagJump, na.rm = TRUE), "\n")
          
            # TJ stats
            TJ <- data.table(
                Total_reads = sum(ASVTAB$Abundance),
                Number_of_TagJump_Events = sum(ASVTAB$TagJump),
                TagJump_reads = sum(ASVTAB[ TagJump == TRUE ]$Abundance, na.rm = T)
                )

            TJ$ReadPercent_removed <- with(TJ, (TagJump_reads / Total_reads * 100))
            fwrite(x = TJ, file = "ASV_table/TagJump_stats.txt", sep = "\t")

            # prepare ASV tables, remove tag-jumps
            ASVTAB <- ASVTAB[ TagJump == FALSE ]
            # convert to wide format
            RES <- dcast(data = ASVTAB,
              formula = ASV ~ SampleID,
              value.var = "Abundance", fill = 0)
            # sort rows (by total abundance)
            clz <- colnames(RES)[-1]
            otu_sums <- rowSums(RES[, ..clz], na.rm = TRUE)
            RES <- RES[ order(otu_sums, decreasing = TRUE) ]

            # output table that is compadible with dada2
            output = as.matrix(RES, sep = "\t", header = TRUE, rownames = 1, 
                                    check.names = FALSE, quote = FALSE)
            output = t(output)
            saveRDS(output, ("ASV_table/ASV_table_TagJumpFiltered.rds"))

            ### format and save ASV table and ASVs.fasta files
            # sequence headers
            asv_seqs = colnames(output)
            asv_headers = openssl::sha1(asv_seqs)
            # transpose sequence table
            toutput = t(output)
            # add sequences to 1st column
            toutput = cbind(row.names(toutput), toutput)
            colnames(toutput)[1] = "Sequence"
            #row names as sequence headers
            row.names(toutput) = asv_headers
            # write ASVs.fasta to path_results
            asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
            write(asv_fasta, file.path(path_results, "ASV_table_TagJumpFiltered.fasta"))
            # write ASVs table to path_results
            write.table(toutput, file.path(path_results, "ASV_table_TagJumpFiltered.txt"), 
                                    sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

            # print summary
            print(paste0("Output = ", length(colnames(output)), " ASVs, with a total of ", 
                                        sum(rowSums(output)), " sequences."))

            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            print("--------")
            setwd(wd)
        i = i + 1
        }
    }


.. _mergeRuns16S:

Merge sequencing runs
~~~~~~~~~~~~~~~~~~~~~

.. code-block:: R
   :caption: merge ASV tables from multiple sequencing runs

    #!/usr/bin/Rscript
    ## Merge sequencing runs, if working with multiple ones; [Sten Anslan]

    # load libraries
    library('dada2')

    # after merging multiple ASV tables ... 
        # collapse ASVs that have no mismatshes or internal indels
    collapseNoMismatch = "true"  #true/false 

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/my_working_folder"
    dirs = list.dirs(recursive = FALSE)
    tables = c()
    # load tables from multiple sequencing runs (dirs)
    for (i in 1:length(dirs)) {
        if(length(dirs) > 1) {
            setwd(dirs[i])
            tables = append(tables, print(file.path(paste0(dirs[i], "/ASV_table"), 
                                                "ASV_table_TagJumpFiltered.rds")))
            setwd(wd)
        i = i + 1
        }
    }

    # Merge multiple ASV tables
    print("# Merging multiple ASV tables ...")
    ASV_tables <- lapply(tables, readRDS)
    merged_table = mergeSequenceTables(tables = ASV_tables, repeats = "sum", tryRC = FALSE)

    ### collapse ASVs that have no mismatshes or internal indels 
    if (collapseNoMismatch == "true") {
        print("# Collapsing identical ASVs ...")
        merged_table_collapsed = collapseNoMismatch(merged_table, 
                                minOverlap = 20, orderBy = "abundance", 
                                identicalOnly = FALSE, vec = TRUE, 
                                band = -1, verbose = TRUE)
        saveRDS(merged_table_collapsed, "merged_table_collapsed.rds")

        ### format and save ASV table and ASVs.fasta files
        # sequence headers
        asv_seqs = colnames(merged_table_collapsed)
        asv_headers = openssl::sha1(asv_seqs)
        # transpose sequence table
        tmerged_table_collapsed = t(merged_table_collapsed)
        # add sequences to 1st column
        tmerged_table_collapsed = cbind(row.names(tmerged_table_collapsed), tmerged_table_collapsed)
        colnames(tmerged_table_collapsed)[1] = "Sequence"
        #row names as sequence headers
        row.names(tmerged_table_collapsed) = asv_headers
        # write ASVs.fasta
        asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
        write(asv_fasta, "ASVs.merged_collapsed.fasta")
        # write ASVs table
        write.table(tmerged_table_collapsed, "ASV_table.merged_collapsed.txt", 
                                sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

        # print summary
        print(paste0("Output = ", length(colnames(merged_table_collapsed)), 
                        " ASVs, with a total of ", 
                        sum(rowSums(merged_table_collapsed)), 
                        " sequences."))
    } else {
        saveRDS(merged_table, "merged_table.rds")
        ### format and save ASV table and ASVs.fasta files
        # sequence headers
        asv_seqs = colnames(merged_table)
        asv_headers = openssl::sha1(asv_seqs)
        # transpose sequence table
        tmerged_table = t(merged_table)
        # add sequences to 1st column
        tmerged_table = cbind(row.names(tmerged_table), tmerged_table)
        colnames(tmerged_table)[1] = "Sequence"
        #row names as sequence headers
        row.names(tmerged_table) = asv_headers
        # write ASVs.fasta to path_results
        asv_fasta <- c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
        write(asv_fasta, "ASVs.merged.fasta")
        # write ASVs table to path_results
        write.table(tmerged_table, "ASV_table.merged.txt", 
                                sep = "\t", col.names = NA, row.names = TRUE, quote = FALSE)

        # print summary
        print(paste0("Output = ", length(colnames(merged_table)), 
                        " ASVs, with a total of ", 
                        sum(rowSums(merged_table)), 
                        " sequences."))
    }


.. _taxAssign16S:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

coming soon ...

.. _clustering16S:

Clustering ASVs to OTUs
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash
   :caption: clustering ASVs to OTUs

   #!/bin/bash
    ## Cluster ASVs to OTUs with vsearch; [Sten Anslan]

    # specify working files (ASVs.fasta and ASV_table)
    fasta=$"ASVs.merged_collapsed.fasta"      # no size annotation (i.e. ;size=)
    ASV_tab=$"ASV_table.merged_collapsed.txt" # 2nd column must be "Sequence"
    export ASV_tab                            # export for R (below)

    #specify clustering threshold
    clustering_thresh="0.97"

    # get ASV size annotation (global sum seqs) from an ASV table
        # & cluster size annotated ASVs to OTUs using vsearch
    awk 'NR>1{for(i=3;i<=NF;i++) t+=$i; print ">"$1";size="t"\n"$2; t=0}' $ASV_tab | \
    vsearch --cluster_fast - \
        --id $clustering_thresh \
        --iddef 2 \
        --sizein \
        --xsize \
        --fasta_width 0 \
        --centroids OTUs.fasta \
        --uc OTUs.uc

.. code-block:: R
   :caption: generating OTU table based on the clustered ASVs

    #!/usr/bin/Rscript
    ## generate OTU table based on the clustered ASVs; (C) Vladimir Mikryukov,
                                                        # edit, Sten Anslan.
    # load library
    library('data.table')

    ## Required inputs
    inp_ASVTAB = Sys.getenv('ASV_tab') # get above specified ASV_table file
                                      # First column witout header
                                      # 2nd col is SEQUENCE COLUMN! 
                                      # NO SIZE ANNOTATION of ASV headers
    ## Load input data - ASV table
    cat("\n ASV_table to OTU_table: input = ", inp_ASVTAB, "\n")
    ASVTAB <- fread(file = inp_ASVTAB, header = TRUE, sep = "\t")
    #drop 2nd col which has seqs
    ASVTAB[[2]] <- NULL

    ## Load input data - UC mapping file
    UC = fread(file = "OTUs.uc", header = FALSE, sep = "\t")
    UC = UC[ V1 != "S" ]
    UC[, ASV := tstrsplit(V9, ";", keep = 1) ]
    UC[, OTU := tstrsplit(V10, ";", keep = 1) ]
    UC[V1 == "C", OTU := ASV ]
    UC = UC[, .(ASV, OTU)]

    #Rename V1 col to "ASV"
    colnames(ASVTAB)[1] = "ASV" 
    ## Convert ASV table to long format
    ASV = melt(data = ASVTAB,
            id.vars = "ASV",
            variable.name = "SampleID", 
            value.name = "Abundance")
    ASV = ASV[ Abundance > 0 ]

    ## Add OTU IDs
    ASV = merge(x = ASV, y = UC, by = "ASV", all.x = TRUE)
    ## Summarize
    OTU = ASV[ , .(Abundance = sum(Abundance, na.rm = TRUE)), by = c("SampleID", "OTU")]

    ## Reshape to wide format
    OTU_tab = dcast(data = ASV,
              formula = OTU ~ SampleID,
              value.var = "Abundance",
              fun.aggregate = sum, fill = 0)
    ## Export OTU table
      # OTU names correspond to most abundant ASV in an OTU cluster
    fwrite(x = OTU_tab, file = "OTU_table.txt", sep = "\t")


.. _postclustering16S:

Post-clustering
~~~~~~~~~~~~~~~

Post-clustering for merging consistently co-occuring 'daughter-OTUs' with `LULU <https://github.com/tobiasgf/lulu>`_. 

.. code-block:: bash
   :caption: make match list database for post-clustering

    #!/bin/bash
    ## make match list database for post-clustering using BLAST; [Sten Anslan] 

    # make blast database
    makeblastdb -in OTUs.fasta -parse_seqids -dbtype nucl
    # generate match list
    blastn -db OTUs.fasta \
            -outfmt '6 qseqid sseqid pident' \
            -out match_list.txt \
            -qcov_hsp_perc 75 \
            -perc_identity 89 \
            -query OTUs.fasta \
            -num_threads 8


.. code-block:: R
   :caption: post-clustering with LULU

    #!/usr/bin/Rscript
    ## run post-clustering with LULU; [Sten Anslan] 

    # load library 
    library('devtools')

    # load OTU table and match list
    otutable    = read.table("OTU_table.txt", header = T, row.names = 1)
    matchlist   = read.table("match_list.txt")
    lulu_result = lulu::lulu(otutable, matchlist, minimum_match = 90) 

    # export post-clustered OTU table
    write.table(lulu_result$curated_table, 
                file ="OTU_table_LULU.txt", 
                sep = "\t", col.names = NA,
                row.names = TRUE, quote = FALSE)
    write.table(lulu_result$discarded_otus, 
                file ="discarded_lulu_OTUs.txt", 
                sep = "\t", col.names = NA, 
                row.names = TRUE, quote = FALSE)


.. code-block:: bash
   :caption: get post-clustered fasta file

    #!/bin/bash
    ## make OTUs_LULU.fasta file corresponding to the OTUs in the LULU filtered table; [Sten Anslan] 

    # drop discarded OTUs
    awk 'NR>1{print $1}' OTU_table_LULU.txt > OTUs.list
    cat OTUs.fasta | seqkit grep -w 0 -f OTUs.list > OTUs_LULU.fasta

    # remove blast database files 
    rm OTUs.fasta.n* 
    rm OTUs.list 

____________________________________________________

|eufund| |chfund| |ukrifund|
