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

.. raw:: html

    <style>
        .yellow-background {
            background-color: #feffd0;
            color: #000000;  /* text color */
            font-size: 20px; /* text size */
            padding: 5px;   /* add  padding */
        }
    </style>

.. role:: yellow-background


|logo_BGE_alpha|


Arthropods/COI
**************

| This is executable step-by-step pipeline for **COI** amplicon data analyses.
|  
| The **full bioinformatics workflow can be automatically run through** `PipeCraft2 <https://pipecraft2-manual.readthedocs.io/en/latest/>`_ (v1.1.0; releasing this soon, with a tutorial),
| which implemets also various error handling processes and sequence summary statistics (lacking here in step-by-step code). 
| 
| The bioinformatic workflow results in amplicon sequence variants (ASVs) and well as operational taxonomic units (OTUs).

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
| :ref:`Chimera filtering <remove_chimerasCOI>`   | DADA2        | 2.28    |
+-------------------------------------------------+--------------+---------+
| :ref:`Remove tag-jumps <tagjumpsCOI>`           | UNCROSS2     | x       |
+-------------------------------------------------+--------------+---------+
| :ref:`Merge sequencing runs* <mergeRunsCOI>`    | DADA2        |   2.28  |
+-------------------------------------------------+--------------+---------+
| :ref:`Remove NUMTs <numtsCOI>`                  | metaMATE     |  0.4.3  |
+-------------------------------------------------+--------------+---------+
| :ref:`Taxonomy assignment <taxAssignCOI>`       | RDP/BLAST    | x       | 
+-------------------------------------------------+--------------+---------+
| :ref:`Clustering ASVs to OTUs <clusteringCOI>`  | vsearch, LULU| x       |
+-------------------------------------------------+--------------+---------+

\*only applicable when there are multiple sequencing runs per study. 


Data structure
~~~~~~~~~~~~~~

.. _multiRunDirCOI:

Multiple sequencing runs
------------------------

.. important:: 

  When aiming to combine samples from multiple sequencing runs, then follow the below directory structure 

**Directory structure:**

| **/multiRunDir** *(directory names can be changed)*
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
  
  Fastq files with the **same name** will be considered as the same sample and will be merged in the "Merge sequencing runs" step.

Single sequencing run
---------------------

| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`
| 

____________________________________________________

.. _remove_primersCOI:

Remove primers
~~~~~~~~~~~~~~

| Remove primer strings from paired-end data.
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`

.. note:: 
  
  Here, assuming that all sequences are in 5'-3' orientation! 
  *(3'-5' orient sequences will be discarded with this workflow)*

.. important:: 

  | - Paired-end files must contain "R1" and "R2" strings (not just _1 and _2)!
  | - Sample names must not contain "R1" and "R2" strings (i.e. not FR123_001_R1.fastq/FR123_001_R2.fastq)

.. code-block:: bash
   :caption: remove primers with cutadapt
   :emphasize-lines: 21-26, 51-52
   :linenos:

    #!/bin/bash
    ## workflow to remove primers via cutadapt

    # My working folder = /multiRunDir (see dir structure above)

    # specify the identifier string for the R1 files
    read_R1="_R1"

    # specify primers 
    fwd_primer=$"GGWACWGGWTGAACWGTWTAYCCYCC"    #this is primer mlCOIintF
    rev_primer=$"TANACYTCNGGRTGNCCRAARAAYCA"    #this is primer jgHCO2198

    # edit primer trimming settings
    maximum_error_rate="1" # Maximum error rate in primer string search;
                           # if set as 1, then allow 1 mismatch;
                           # if set as 0.1, then allow mismatch in 10% of the bases,
                           # i.e. if a primer is 20 bp then allowing 2 mismatches.
    overlap="22"           # The minimum overlap length. Keep it nearly as high
                           # as the primer length to avoid short random matches.

    # get directory names if working with multiple sequencing runs
    DIRS=$(ls -d *) # -> sequencing_set01 sequencing_set02 sequencing_set03

    for sequencing_run in $DIRS; do 
        printf "\nWorking with $sequencing_run \n"
        cd $sequencing_run
        #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
        # make output dirs
        mkdir -p primersCut_out
        mkdir -p primersCut_out/untrimmed

        ### Clip primers with cutadapt
        for inputR1 in *$read_R1*; do
            inputR2=$(echo $inputR1 | sed -e 's/R1/R2/')
            cutadapt --quiet \
            -e $maximum_error_rate \
            --minimum-length 32 \
            --overlap $overlap \
            --no-indels \
            --cores=0 \
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


.. _quality_filteringCOI:

Quality filtering 
~~~~~~~~~~~~~~~~~

| Quality filtering of the fastq files based on the allowed maximum error rate per sequence (as in DADA2).
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`

.. code-block:: R
   :caption: quality filtering in DADA2 (in R)
   :emphasize-lines: 13-19, 64-68
   :linenos:

    #!/usr/bin/Rscript
    ## workflow to perform quality filtering within DADA2

    #load dada2 library 
    library('dada2')

    # specify the identifier string for the R1 files
    read_R1 = ".R1"
    
    # get the identifier string for the R2 files
    read_R2 = gsub("R1", "R2", read_R1)

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/multiRunDir"
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
            sample_names = sapply(strsplit(basename(R1s), read_R1), `[`, 1)

            # filtered files path
            filtR1 = file.path(path_results, paste0(sample_names, ".R1.", "fastq.gz"))
            filtR2 = file.path(path_results, paste0(sample_names, ".R2.", "fastq.gz"))
            names(filtR1) = sample_names
            names(filtR2) = sample_names
            
            #quality filtering
            qfilt = filterAndTrim(R1s, filtR1, R2s, filtR2, 
                                maxN = 0,            # max number of allowed N bases.
                                maxEE = c(2, 2),     # max error rate per R1 and R2 read, respectively.
                                truncQ = 2,          # truncate reads at the first instance of a quality score less than or equal to specified value. 
                                truncLen = c(0, 0),  # truncate reads after specified length for R1 and R2 reads, respectively.
                                maxLen = 600,        # discard reads longer than specified.
                                minLen = 100,        # discard reads shorter than specified.
                                minQ = 2,            # discard reads (after truncation) that contain a quality score below specified value.
                                matchIDs = TRUE,     # output paired-end reads with matching IDs (for merging).
                                compress = TRUE,     # gzip the output
                                multithread = TRUE)  # use multiple threads
            saveRDS(qfilt, file.path(path_results, "qfilt_reads.rds"))

            # make sequence count report
            seq_count = cbind(qfilt)
            colnames(seq_count) = c("input", "qualFiltered")
            seq_count = as.data.frame(seq_count)
            seq_count$sample = sample_names
            # reorder columns
            seq_count = seq_count[, c("sample", "input", "qualFiltered")]
            write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), 
                                row.names = FALSE, quote = FALSE)

            # save filtered R objects for denoising and merging (below)
            filtR1 = sort(list.files(path_results, pattern = ".R1.fastq.gz", full.names = TRUE))
            filtR2 = sort(list.files(path_results, pattern = ".R2.fastq.gz", full.names = TRUE))
            sample_names = sapply(strsplit(basename(filtR1), ".R1.fastq.gz"), `[`, 1)
            saveRDS(filtR1, file.path(path_results, "filtR1.rds"))
            saveRDS(filtR2, file.path(path_results, "filtR2.rds"))
            saveRDS(sample_names, file.path(path_results, "sample_names.rds"))
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            #set working directory back to "/multiRunDir"
            setwd(wd)
        i = i + 1
        }
    }


.. _denoiseCOI:

Denoise and merge paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


| Denoise and merge paired-end Illumina reads as in DADA2.
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`


.. code-block:: R
   :caption: denoise and merge paired-end reads in DADA2
   :emphasize-lines: 7-13, 75-79
   :linenos:

    #!/usr/bin/Rscript
    ## workflow to perform DADA2 denoising and merging

    # load dada2 library 
    library('dada2')

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/multiRunDir"
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
            getN = function(x) sum(getUniques(x))
            #remove 0 seqs samples from qfilt statistics
            row_sub = apply(qfilt, 1, function(row) all(row !=0 ))
            qfilt = qfilt[row_sub, ]
            seq_count = cbind(qfilt, sapply(dadaR1, getN), 
                                sapply(dadaR2, getN), sapply(merge, getN))
            colnames(seq_count) = c("input", "qualFiltered", "denoised_R1", "denoised_R2", "merged")
            rownames(seq_count) = sample_names
            write.csv(seq_count, file.path(path_results, "seq_count_summary.csv"), 
                                    row.names = TRUE, quote = FALSE)
            #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
            print("--------")
            setwd(wd)
        i = i + 1
        }
    }



.. _remove_chimerasCOI:

Chimera filtering 
~~~~~~~~~~~~~~~~~

| Remove putative chimeras with DADA2 'consensus' mode.
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`

.. code-block:: R
   :caption: remove chimeras in DADA2
   :emphasize-lines: 14-20, 97-100
   :linenos:

    #!/usr/bin/Rscript
    ## workflow to perform chimera filtering within DADA2

    # load libraries
    library('dada2')
    library('openssl')

    # chimera filtering method
    method = "consensus" 

    # collapse ASVs that have no mismatshes or internal indels (identical up to shifts and/or length)
    collapseNoMismatch = "true"  #true/false 

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/multiRunDir"
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
            asv_fasta = c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
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
                saveRDS(ASV_tab_collapsed, file.path(path_results, "ASV_table_collapsed.rds"))

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
                asv_fasta = c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
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



.. _tagjumpsCOI:

Remove tag-jumps
~~~~~~~~~~~~~~~~

| Remove putative tag-jumps with UNCROSS2.
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`

.. code-block:: R
   :caption: removing putative tag-jumps with UNCROSS2 method
   :emphasize-lines: 12-18, 112-116
   :linenos:

   #!/usr/bin/Rscript
   ## Script to perform tag-jump removal; (C) Vladimir Mikryukov,
                                             # edit, Sten Anslan

    # load libraries
    library(data.table)

    # set parameters
    set_f = 0.03 # f-parameter of UNCROSS (e.g., 0.03)
    set_p = 1    # p-parameter (e.g., 1.0)

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/multiRunDir"
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
            } else { # loading chimera filtered ASV table
              tab = readRDS("ASV_table/chim_filt.rds")
              cat("input table = ASV_table/chim_filt.rds\n")
            }

            # format ASV table
            ASVTABW = as.data.table(t(tab), keep.rownames = TRUE)
            colnames(ASVTABW)[1] = "ASV"
            # convert to long format
            ASVTAB = melt(data = ASVTABW, id.vars = "ASV",
            variable.name = "SampleID", value.name = "Abundance")
            # remove zero-OTUs
            ASVTAB = ASVTAB[ Abundance > 0 ]
            # estimate total abundance of sequence per plate
            ASVTAB[ , Total := sum(Abundance, na.rm = TRUE), by = "ASV" ]

            ## UNCROSS score
            uncross_score = function(x, N, n, f = 0.01, tmin = 0.1, p = 1){
              z = f * N / n               # Expected treshold
              sc = 2 / (1 + exp(x/z)^p)   # t-score
              res = data.table(Score = sc, TagJump = sc >= tmin)
              return(res)
            }

            # esimate UNCROSS score
            cat(" estimating UNCROSS score\n")
            ASVTAB = cbind(
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
          
            # tag-jump stats
            TJ = data.table(
                Total_reads = sum(ASVTAB$Abundance),
                Number_of_TagJump_Events = sum(ASVTAB$TagJump),
                TagJump_reads = sum(ASVTAB[ TagJump == TRUE ]$Abundance, na.rm = T))

            TJ$ReadPercent_removed = with(TJ, (TagJump_reads / Total_reads * 100))
            fwrite(x = TJ, file = "ASV_table/TagJump_stats.txt", sep = "\t")

            # prepare ASV tables, remove tag-jumps
            ASVTAB = ASVTAB[ TagJump == FALSE ]
            # convert to wide format
            RES = dcast(data = ASVTAB,
              formula = ASV ~ SampleID,
              value.var = "Abundance", fill = 0)
            # sort rows (by total abundance)
            clz = colnames(RES)[-1]
            otu_sums = rowSums(RES[, ..clz], na.rm = TRUE)
            RES = RES[ order(otu_sums, decreasing = TRUE) ]

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
            asv_fasta = c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
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



.. _mergeRunsCOI:

Merge sequencing runs
~~~~~~~~~~~~~~~~~~~~~

| If previous processing was applied on :ref:`multiple sequencing runs <multiRunDirCOI>` , then here, 
| merge those sequenceing runs to form a single, unified ASV table. 
| Assuming that tag-jump filtering was performed (inputs = ASV_table_TagJumpFiltered.rds)

.. code-block:: R
   :caption: merge ASV tables from multiple sequencing runs
   :emphasize-lines: 1-88
   :linenos:

    #!/usr/bin/Rscript
    ## Merge sequencing runs, if working with multiple ones

    # load libraries
    library('dada2')

    # after merging multiple ASV tables ... 
        # collapse ASVs that have no mismatshes or internal indels
    collapseNoMismatch = "true"  #true/false 

    # capturing the directory structure when working with multiple runs
    wd = getwd() # -> wd is "~/multiRunDir"
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
    ASV_tables = lapply(tables, readRDS)
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
        asv_fasta = c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
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
        asv_fasta = c(rbind(paste(">", asv_headers, sep=""), asv_seqs))
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



.. _numtsCOI:

Remove NUMTs
~~~~~~~~~~~~

| Remove putative NUMTs with metaMATE. 
| This follows the workflow to automatically filter the ASVs by retaining maximum of 5% of estimated non-authentic-ASVs (nonauthentic_retained_estimate_p < 0.05).


.. important::

  1. metaMATE expects specifications file that states the filtering strategies. See `more info here. <https://github.com/tjcreedy/metamate?tab=readme-ov-file#specifications>`_ 
  Here, we will be using the metaMATE's `default specifications.txt file. <https://github.com/tjcreedy/metamate/blob/main/specifications.txt>`_ 

  1. metaMATE requires a reference COI database to determine verified-authentic ASVs. Herein using `CO1Classifier v5.1.0 database. <https://github.com/terrimporter/CO1Classifier>`_ 
  
  --- `Download the CO1Classifier v5.1.0 database here (click) <https://github.com/terrimporter/CO1Classifier/releases/download/SINTAX-COI-v5.1.0-ref/SINTAX_COIv5.1.0_ref.zip>`_ ---


Check `standard genetic codes here <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>`_ for ``genetic_code`` setting below.

.. code-block:: bash
   :caption: get required specifications file and ref database
   :linenos:

   #!/bin/bash
    
    # download the default specifications file, 
      # using this in metaMATE-find
    wget "https://raw.githubusercontent.com/tjcreedy/metamate/main/specifications.txt"
    # specify specifications file for metaMATE (with full directory path)
    specifications=$(realpath specifications.txt)


    # download the CO1Classifier reference databse
    wget "https://github.com/terrimporter/CO1Classifier/releases/download/SINTAX-COI-v5.1.0-ref/SINTAX_COIv5.1.0_ref.zip"
    # unzip the database and edit name
    unzip SINTAX_COIv5.1.0_ref.zip && mv training CO1Classifier_v5.1.0 
    mv CO1Classifier_v5.1.0/sintax.fasta CO1Classifier_v5.1.0/CO1Classifier_v5.1.0.fasta
    
    # specify reference database for metaMATE (with full directory path)
    reference_database=$(realpath CO1Classifier_v5.1.0/CO1Classifier_v5.1.0.fasta)


.. code-block:: bash
   :caption: run metaMATE-find
   :linenos:

    #!/bin/bash
    ## remove NUMTs with metaMATE
  
    ## go to the directory that hosts your ASVs.fasta and ASV table files.
    cd ASV_table # for example
    
    # specify input ASVs table and fasta
    ASV_table="ASVs_table_collapsed.txt"   # specify ASV table file 
    ASV_fasta="ASVs_collapsed.fasta"       # specify ASVs fasta file 

    # specify variables
    genetic_code="5"        # the standard genetic code. 5 is invertebrate mitochondrial code
    length="313"            # the expected length of an amplicon
    NA_abund_thresh="0.05"  # nonauthentic_retained_estimate_p threshold 
    basesvariation="3"      # allowed length variation (bp) from the expected length of an amplicon
    taxgroups="undefined"   # (optional); if sequence binning is to be performed on 
                               # a per-taxon basis (as in specifications file) 
                                  # then specify the taxon grouping file
    ## 

    # check if taxgroups is specified, if not then this var is empty.
    if [[ $taxgroups != "undefined" ]]; then
        taxgroups=$"--taxgroups $taxgroups"
    else 
        taxgroups=$""
    fi

    #output dir
    output_dir=$"metamate_out"
    echo "output_dir = $output_dir"

    # if perfoming clade binning, then WARNING when processing more than 65,536 ASVs
    ASVcount=$(grep -c "^>" $ASV_fasta)
    if (( $ASVcount > 65536 )); then
        printf '%s\n' "WARNING]: clade binning NOT performed, 
         because the input ASVs limit is 65,536 for that.
         Current input has $ASVcount ASVs."
    fi

    # quick check of the specifications file, has to contain "library" | "total" | "clade" | "taxon"
    if ! grep -q -e "library" -e "total" -e "clade" -e "taxon" $specifications; then
        printf '%s\n' "ERROR]: specifications file seems to be wrong. 
         Does not contain any of the terms (library, total, clade, taxon)."
    fi

    # remove old $output_dir if exists
    if [[ -d $output_dir ]]; then
        rm -rf $output_dir
    fi

    ### metaMATE-find
    printf "# Running metaMATE-find\n"
    metamate find \
        --asvs $ASV_fasta \
        --readmap $ASV_table \
        --specification $specifications \
        --references $reference_database \
        --expectedlength $length \
        --basesvariation $basesvariation \
        --onlyvarybycodon \
        --table $genetic_code \
        --threads 8 \
        --output $output_dir \
        --overwrite $taxgroups

    # check for the presence of "metamate_out" dir and "resultcache" file (did metaMATE-find finish)
    if [[ -d $output_dir ]] && [[ -e $output_dir/resultcache ]] && [[ -e $output_dir/results.csv ]]; then
        printf '%s\n' "metaMATE-find finished"
        # export variables for below script (Rscript)
        export NA_abund_thresh
        export output_dir
    else 
        printf '%s\n' "ERROR]: cannot find the $output_dir (metaMATE-find output) 
         to start metaMATE-dump"
    fi


.. code-block:: bash
   :caption: get the results_index from the metamate_out/results.csv file that corresponds to the specified 'NA_abund_thresh'
   :linenos:

    #!/usr/bin/env Rscript

    # NA_abund_thresh is the allowed abundance threshold of 
       # non-validated (putative artefactual) OTUs/ASVs in the filtered dataset.
    ## read results.csv
    NA_abund_thresh = as.numeric(Sys.getenv('NA_abund_thresh'))
    output_dir = Sys.getenv('output_dir')
    find_results = read.csv(file.path(output_dir, "results.csv"))

    ## filter results based on NA_abund_thresh 
    filtered_data = find_results[
                        find_results$nonauthentic_retained_estimate_p <= NA_abund_thresh, ] 

    # if no results correspond with the NA_abund_thresh, then get the next best
        # else, just select the result_index that corresponds to 
            # NA_abund_thresh with highest accuracy_score
    if (nrow(filtered_data) == 0) {
        cat(
          "\n no results correspond with the NA_abund_thresh of", NA_abund_thresh, "; 
          getting the next best setting\n"
          )
        next_best = min(find_results$nonauthentic_retained_estimate_p)
        filtered_data = find_results[
                          find_results$nonauthentic_retained_estimate_p <= next_best, ] 
        # sort based on accuracy_score
        sorted_filtered = filtered_data[order(-filtered_data$accuracy_score), ]
        # get the result with the highest accuracy_score
        metamate_selected_threshold = sorted_filtered[1,]
        write.csv(metamate_selected_threshold, file.path(output_dir, "next_best_set.csv"), 
                                              quote = F)
        # the result_index of the NA_abund_thresh with the highest accuracy_score
        result_index = metamate_selected_threshold[,1]
        write(result_index, file.path(output_dir, "selected_result_index.txt"))
    } else {
        # sort based on accuracy_score
        sorted_filtered = filtered_data[order(-filtered_data$accuracy_score), ]
        # get the result with the highest accuracy_score
        metamate_selected_threshold = sorted_filtered[1,]
        # the result_index of the NA_abund_thresh with the highest accuracy_score
        result_index = metamate_selected_threshold[,1]
        write(result_index, file.path(output_dir, "selected_result_index.txt"))
    }



.. code-block:: bash
   :caption: run metaMATE-dump to discard putative artefact ASVs
   :linenos:

    ## metaMATE-dump 
    dump_seqs=$(basename $ASV_fasta) # output file name
    # check for the presence of "metamate_out" dir and "resultcache" file (did metaMATE-find finish)
    printf "# Running metaMATE-dump\n"

    # read result_index
    read -r result_index < $output_dir/selected_result_index.txt
    printf " - selcted result_index = $result_index\n"

    # run metaMATE-dump
    metamate dump \
    --asvs $ASV_fasta \
    --resultcache $output_dir/resultcache \
    --output $output_dir/${dump_seqs%.*}_metaMATE.filt \
    --overwrite \
    --resultindex $result_index

    # generate a list of ASV IDs 
    seqkit seq -n $output_dir/${dump_seqs%.*}_metaMATE.filt.fasta > \
                        $output_dir/${dump_seqs%.*}_metaMATE.filt.list

    # filter the ASV table; include only the ASVs that are in ${dump_seqs%.*}_metaMATE.filt.list
    out_table=$(basename $ASV_table)
    awk -v var="$output_dir/${dump_seqs%.*}" 'NR==1; NR>1 {print $0 | "grep -Fwf "var"_metaMATE.filt.list"}' $ASV_table > \
                                                                              $output_dir/${out_table%.*}_metaMATE.filt.txt


.. _taxAssignCOI:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

| Assign taxonomy with **RDP-classifier and BLAST**. 
| The script below will first assign taxonomy to the metaMATE output fasta file with RDP-classifer;
| then the RDP-classifier results will be sorted to discard ASVs with a bootstrap (~assignment confidence) value < 0.8 to kingdom Metazoa. 
| The ASVs with at least 0.8 bootstrap value against Metazoa are further subjected to BLASTn search.
| **Non-Metazoan ASVs are removed** at that stage.

.. code-block:: bash
   :caption: at first, assign taxonomy with RDP-classifier
   :linenos:

    #!/bin/bash

    # USAGE: script.sh input.fasta input_OTU_table.txt BLAST
            # if $3 = BLAST, then perform BLAST
    # input.fasta = OTUs
    # input_OTU_table.txt = OTU table; may contain seqs as a 2nd col (used in metazoa sort part)

    # specify reference database. Using the same as we used for metaMATE
          # but in a format suitable for RDP-classifer 
    DB="/home/sten/Desktop/DATABASES/COI_classifier_5.1.0/rRNAClassifier.properties"
    printf "\n # Running RDP \n"
    time rdp_classifier \                  
    -Xmx12g \                              # Max MEM usage (Xmx12g = 12G)
    classify \                             # RDP classify
    -t $DB \                               # database file
    -o RDP.taxonomy.txt \                  # outout file name
    -q ${dump_seqs%.*}_metaMATE.filt.fasta # input file name (the output of metaMATE)

.. code-block:: bash
   :caption: Get only metazoa annotations (bootstrap >0.8)
   :linenos:

    # output = tax_meatazoa.csv and table_metazoa.csv
    eval "$(conda shell.bash hook)"
    conda activate dada2
    printf "\n # Sorting Metazoa \n"
    Rscript /home/sten/Desktop/HTS_data/SilvaNova/pipe_SNC/COI_postprocessing/RDP_Metazoa_sort2.R RDP.taxonomy.txt $2

    # get only metazoa fasta based on tax_metazoa.csv #
    seqkit grep -f <(awk -F ',' '{print $1}' tax_metazoa.csv) -w 0 OTUs_LULU.ORFs.fasta > OTUs_metazoa.fasta

.. code-block:: bash
   :caption: BLAST OTUs_metazoa.fasta
   :linenos:

    eval "$(conda shell.bash hook)"
    conda activate LULU #blast work tested in this env, blast = BLAST 2.11.0+

    DB="/home/sten/Desktop/DATABASES/COI_classifier_5.1.0/for_SINTAX/COIv5.1.DB"
    #query
    fasta=$"OTUs_metazoa.fasta"
    seqcount=$(grep -c "^>" $fasta)

    #BLAST
    printf "\n # Running BLAST for $seqcount seqs \n\n"
    time blastn -strand plus -num_threads 8 \
    -query $fasta \
    -db $DB \
    -out 10BestHits.txt -task blastn \
    -max_target_seqs 10 -evalue=0.001 -word_size=7 -reward=1 -penalty=-1 -gapopen=1 -gapextend=2 \
    -outfmt "6 delim=+ qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident"

    printf "#qseqid = Query Seq-id
    #qlen = Query sequence length
    #sacc = Subject accession
    #slen = Subject sequence length
    #qstart = Start of alignment in query
    #qend = End of alignment in query
    #sstart = Start of alignment in subject
    #send = End of alignment in subject
    #evalue = Expect value
    #length = Alignment length
    #pident = Percentage of identical matches
    #nident = Number of identical matches
    #mismatch = Number of mismatches
    #gapopen = Number of gap openings
    #gaps = Total number of gaps
    #1st_hit = BLAST 1st hit
    #sstrand = Subject Strand
    #qcovs = Query Coverage Per Subject" > README.txt

   # + to \t
    sed -i 's/+/\t/g' 10BestHits.txt
    #get only first occurrence of a duplicate row (1st hit)
    #awk 'BEGIN{FS="+"}''!seen[$1]++' 10BestHits.txt > 1.temphit #sep = +
    awk 'BEGIN{FS="\t"}''!seen[$1]++' 10BestHits.txt > 1.temphit #sep = \t

    #check which seqs got a hit
    gawk 'BEGIN{FS="\t"}{print $1}' < 1.temphit | uniq > gothits.names
    #add no_hits flag
    seqkit seq -n $fasta > $fasta.names
    grep -v -w -F -f gothits.names $fasta.names | sed -e 's/$/\tNo_significant_similarity_found/' >> 1.temphit
    #add header
    sed -e '1 i\qseqid+1st_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident' 
    1.temphit > BLAST_1st_hit.txt
    #+ to tab
    sed -i 's/+/\t/g' BLAST_1st_hit.txt


    ###10hits
    #do the same for next 2-10 hits
    time for i in {2..10}; do
        awk -v i="$i" 'BEGIN{FS="\t"}''++seen[$1]==i' 10BestHits.txt > $i.temphit
        gawk 'BEGIN{FS="\t"}{print $1}' < $i.temphit | uniq > gothits.names
        grep -v -w -F -f gothits.names $fasta.names | sed -e 's/$/\tNo_BLAST_hit/' >> $i.temphit && rm gothits.names
    done

    #sort
    time for file in *.temphit; do
        sort -k 1 --field-separator=\t $file > $file.temp && rm $file
    done

    #merge
    paste 1.temphit.temp 2.temphit.temp 3.temphit.temp 4.temphit.temp 5.temphit.temp 6.temphit.temp 7.temphit.temp 8.temphit.temp 
    9.temphit.temp 10.temphit.temp > BLAST_10_hits.txt
    rm *.temp

    #format 10 hits
    sed -i 's/No_significant_similarity_found.*/No_significant_similarity_found/' BLAST_10_hits.txt
    sed -i 's/No_BLAST_hit.*/No_BLAST_hit/' BLAST_10_hits.txt
    sed -i '1i\qseqid+1st_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+2nd_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+3rd_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+4th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+5th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+6th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+7th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+8th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+9th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident+qseqid+10th_hit+qlen+slen+qstart+qend+sstart+send+evalue+length+nident+mismatch+gapopen+gaps+sstrand+qcovs+pident' BLAST_10_hits.txt
    sed -i 's/+/\t/g' BLAST_10_hits.txt


    ##### BLAST 1st hit with query SEQ ######
    #fasta to oneline
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' < $fasta | sed -e 's/\r//' > $fasta.oneline

    #sort hits
    sort -k 1 --field-separator=\t BLAST_1st_hit.txt > BLAST_1st_hit_with_qseq.temp
    sed -i 's/qseqid.*//' BLAST_1st_hit_with_qseq.temp
    sed -i '/^$/d' BLAST_1st_hit_with_qseq.temp
    sort -k 1 --field-separator=\t $fasta.oneline | sed -e 's/^>//' | sed -e 's/\r//' > seqs.txt
    #merge seqs and 1st hit
    paste seqs.txt BLAST_1st_hit_with_qseq.temp > BLAST_1st_hit_with_qseq.txt && rm BLAST_1st_hit_with_qseq.temp
    sed -i '1 i\qseqid\tquery_seq\tqseqid\t1st_hit\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tlength\tnident\tmismatch\tgapopen\tgaps\tsstrand\tqcovs\tpident' BLAST_1st_hit_with_qseq.txt

    rm seqs.txt
    rm $fasta.oneline
    rm *.names

    wc -l BLAST_1st_hit_with_qseq_sep=+.txt
    grep -c "^>" $fasta


.. _nonMetazoaCOI:


.. _clusteringCOI:

Clustering ASVs to OTUs
~~~~~~~~~~~~~~~~~~~~~~~

| Clustering ASVs to OTUs with vsearch. 
| Applying also post-clustering with LULU to merge potential "daughter-OTUs".

.. code-block:: bash
   :caption: clustering
   :linenos:

    ### Get ASV size annotation (global sum seqs) from an ASV table.
    out=$(basename $fasta | awk 'BEGIN{FS=OFS="."}NF{NF -=1}1')
    awk 'NR>1{for(i=3;i<=NF;i++) t+=$i; print ">"$1";size="t"\n"$2; t=0}' $ASV_tab > \
                                                                        $out.size.fasta

    ### Cluster ASVs using vsearch.
    vsearch --cluster_fast $out.size.fasta \
        --id $clustering_thresh \
        --iddef 2 \
        --sizein \
        --xsize \
        --fasta_width 0 \
        --centroids OTUs.fasta \
        --uc OTUs.uc

.. code-block:: R
   :caption: making OTU table (based on OTUs.uc and $ASV_tab)
   :linenos:

    Rscript $pipe_path/ASVs2OTUs.R

.. code-block:: bash
   :caption: LULU post-clustering
   :linenos:

    eval "$(conda shell.bash hook)"
    conda activate LULU

    #make blast db
    makeblastdb -in OTUs.fasta -parse_seqids -dbtype nucl

    #generate match list
    blastn -db OTUs.fasta \
    -outfmt '6 qseqid sseqid pident' \
    -out match_list.txt \
    -qcov_hsp_perc 75 \
    -perc_identity 90 \
    -query OTUs.fasta \
    -num_threads 8

    #Run LULU in R
    Rscript $pipe_path/lulu.R
    wait

    # Drop discarded OTUs
    awk 'NR>1{print $1}' OTU_table_LULU.txt > OTUs.list
    cat OTUs.fasta | seqkit grep -w 0 -f OTUs.list > OTUs_LULU.fasta

    rm OTUs.fasta.n* #remove blast database


____________________________________________________

|eufund| |chfund| |ukrifund|
