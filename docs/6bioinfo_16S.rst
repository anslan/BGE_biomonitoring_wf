.. |logo_BGE_alpha| image:: _static/logo_BGE_alpha.png
  :width: 300
  :alt: Alternative text
  :target: https://biodiversitygenomics.eu/

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

.. raw:: html

    <style>
        .yellow-background {
            background-color: #feffd0;
            color: #000000;  /* text color */
            font-size: 18px; /* text size */
            padding: 5px;   /* add  padding */
        }
    </style>

.. role:: yellow-background


|logo_BGE_alpha|


Bacteria/16S
************

| This is executable step-by-step pipeline for **rRNA 16S** amplicon data from Illumina sequencing machine.
|  
| The **full bioinformatics workflow can be automatically run through** `PipeCraft2 <https://pipecraft2-manual.readthedocs.io/en/latest/>`_ (v1.1.0; releasing this soon, with a tutorial), which implemets also various error handling processes and sequence summary statistics (lacking here in step-by-step code). 
| 
| The bioinformatic workflow results in amplicon sequence variants (ASVs) and well as operational taxonomic units (OTUs).


Dependencies
~~~~~~~~~~~~

+-----------------------------------------------------+----------+----------------+
| Process                                             | Software | Version        |
+=====================================================+==========+================+
| :ref:`Remove primers <remove_primers16S>`           | cutadapt | 4.9            |
+-----------------------------------------------------+----------+----------------+
| :ref:`Quality filtering <quality_filtering16S>`     | DADA2    | 1.30           |
+-----------------------------------------------------+----------+----------------+
| :ref:`Denoise <denoise16S>`                         | DADA2    | 1.30           |
+-----------------------------------------------------+----------+----------------+
| :ref:`Merge paired-end reads <denoise16S>`          | DADA2    | 1.30           |
+-----------------------------------------------------+----------+----------------+
| :ref:`Chimera filtering <remove_chimeras16S>`       | DADA2    | 1.30           |
+-----------------------------------------------------+----------+----------------+
| :ref:`Remove tag-jumps <tagjumps16S>`               | UNCROSS2 |                |
+-----------------------------------------------------+----------+----------------+
| :ref:`Merge sequencing runs* <mergeRuns16S>`        | DADA2    | 1.30           |
+-----------------------------------------------------+----------+----------------+
| :ref:`Taxonomy assignment <taxAssign16S>`           | RDP      | in DADA2 v1.30 |
+-----------------------------------------------------+----------+----------------+
| :ref:`Clustering ASVs to OTUs <clustering16S>`      | vsearch  | 2.28.1         |
+-----------------------------------------------------+----------+----------------+
| :ref:`Post-clusteringlustering <postclustering16S>` | LULU     | 0.1.0          |
+-----------------------------------------------------+----------+----------------+

\*only applicable when there are multiple sequencing runs per study. 


.. note::

    All the dependencies/software of the pipeline are available on a `Docker image <https://hub.docker.com/r/pipecraft/bioscanflow>`_.

| Download `Docker for windows <https://www.docker.com/get-started>`_ 
| Download `Docker for Mac <https://www.docker.com/get-started>`_ 
| Install Docker for Linux - `follow the guidelines under appropriate Linux distribution <https://docs.docker.com/engine/install/ubuntu/>`_

.. code-block:: bash
   :caption: get the Docker image
   
   docker pull pipecraft/bioscanflow:1

.. code-block:: bash
   :caption: example of running the pipeline via Docker image
   
   # run docker 
    # specify the files location with -v flag  ($PWD = the current working directory)
   docker run -i --tty -v $PWD/:/Files pipecraft/bioscanflow:1 

   # inside the container, the files are accessible in the /Files dir
   cd Files

   # checking if cutadapt is available
   cutadapt -h 

   # ready to run the pipe as below ...
    ## make sure that via the shared folder (-v) path you have access also to the reference databases.



Data structure
~~~~~~~~~~~~~~

.. _multiRunDir16S:

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

.. _remove_primers16S:

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
    fwd_primer=$"GTGYCAGCMGCCGCGGTAA"    #this is primer 515F
    rev_primer=$"GGCCGYCAATTYMTTTRAGTTT" #this is primer 926R

    # edit primer trimming settings
    maximum_error_rate="2" # Maximum error rate in primer string search;
                           # if set as 1, then allow 1 mismatch;
                           # if set as 0.1, then allow mismatch in 10% of the bases,
                           # i.e. if a primer is 20 bp then allowing 2 mismatches.
    overlap="19"           # The minimum overlap length. Keep it nearly as high
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

____________________________________________________

.. _quality_filtering16S:

Quality filtering 
~~~~~~~~~~~~~~~~~

| Quality filtering of the fastq files based on the allowed maximum error rate per sequence (as in DADA2).
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`

.. code-block:: R
   :caption: quality filtering in DADA2 (in R)
   :emphasize-lines: 13-19, 67-71
   :linenos:

    #!/usr/bin/Rscript
    ## workflow to perform quality filtering within DADA2

    #load dada2 library 
    library('dada2')

    # specify the identifier string for the R1 files
    read_R1 = "_R1"
    
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

____________________________________________________

.. _denoise16S:

Denoise and merge paired-end reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| Denoise and merge paired-end Illumina reads as in DADA2.
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`


.. code-block:: R
   :caption: denoise and merge paired-end reads in DADA2
   :emphasize-lines: 7-13, 74-78
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
            cat("\n sample names = ", sample_names, "\n ")
            names(filtR1) = sample_names
            names(filtR2) = sample_names

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

            # Sample inference and merger of paired-end reads
            mergers = vector("list", length(sample_names))
            names(mergers) = sample_names
            for(sample in sample_names) {
              cat("\n -- Processing:", sample, "\n")
              derepF = derepFastq(filtR1[[sample]])
              ddF = dada(derepF, err = errF, multithread = TRUE)
              derepR = derepFastq(filtR2[[sample]])
              ddR = dada(derepR, err = errR, multithread = TRUE)
              merger = mergePairs(ddF, derepF, ddR, derepR)
              mergers[[sample]] = merger
            }
            rm(derepF); rm(derepR)
            gc()
            saveRDS(mergers, (file.path(path_results, "mergers.rds")))

            # make sequence table
            ASV_tab = makeSequenceTable(mergers)
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

____________________________________________________

.. _remove_chimeras16S:

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

____________________________________________________

.. _tagjumps16S:

Remove tag-jumps
~~~~~~~~~~~~~~~~

| Remove putative tag-jumps with UNCROSS2.
|
| When working with a **single directory** that hosts your fastq files, then
| :yellow-background:`ignore (do not execute) the script lines in yellow.`

.. code-block:: R
   :caption: removing putative tag-jumps with UNCROSS2 method
   :emphasize-lines: 15-21, 115-119
   :linenos:

   #!/usr/bin/Rscript
   ## Script to perform tag-jump removal; (C) Vladimir Mikryukov,
                                             # edit, Sten Anslan

    # load libraries
    library(data.table)

    # set parameters
    set_f = 0.03 # f-parameter of UNCROSS (e.g., 0.03)
    set_p = 1    # p-parameter (e.g., 1.0)

    # output dir
    path_results="ASV_table"

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
            write(asv_fasta, file.path(path_results, "ASVs_TagJumpFiltered.fasta"))
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

____________________________________________________

.. _mergeRuns16S:

Merge sequencing runs
~~~~~~~~~~~~~~~~~~~~~

| If previous processing was applied on :ref:`multiple sequencing runs <multiRunDir16S>` , then here, 
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

____________________________________________________

.. _taxAssign16S:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

| Assign taxonomy with **RDP classifier** (assignTaxonomy function in DADA2). 
| Herein using the `SILVA138.1 reference database <https://zenodo.org/records/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz?download=1>`_ 

.. code-block:: bash
   :caption: assignTaxonomy
   :linenos:

    #!/usr/bin/env Rscript

    # specify the query fasta file
    fasta=$"ASVs_TagJumpFiltered.fasta"
    # specify reference database 
    reference_database="silva_nr99_v138.1_wSpecies_train_set.fa.gz"

    # load fasta file
    library("seqinr")
    fasta = read.fasta(fasta, seqtype = "DNA", 
                        as.string = TRUE, 
                        forceDNAtolower = FALSE, 
                        seqonly = FALSE)
    seq_names = getName(fasta)
    seqs = unlist(getSequence(fasta, as.string = TRUE))

    # print number of ASVs
    paste("Number of input sequences = ", length(seq_names))

    # assign taxonomy with dada2 'assignTaxonomy'
    library("dada2")
    set.seed(1)
    tax = assignTaxonomy(seqs, 
            reference_database, 
            multithread = TRUE, 
            minBoot = 80, 
            tryRC = TRUE, 
            outputBootstraps = TRUE)

    # format and export taxonomy table
    tax_out = cbind(seq_names, tax$tax, tax$boot)
    colnames(tax_out)[1] = "Seq_name"
    write.table(tax_out, "taxonomy.txt", 
                sep = "\t", 
                col.names = NA, 
                row.names = TRUE, 
                quote = FALSE)

    ### output = taxonomy.txt

____________________________________________________

.. note:: 

    The final ASVs data is ``ASV_table_TagJumpFiltered.txt`` and ``ASVs_TagJumpFiltered.fasta`` in the ``ASV_table`` directory. 
    In case of Merging multiple sequencing runs, the ASVs data is ``ASV_table.merged.txt`` / ``ASV_table.merged_collapsed.txt`` and ``ASVs.merged.fasta`` / ``ASVs.merged_collapsed.fasta`` in the ``multiRunDir``.

____________________________________________________

.. _clustering16S:

Clustering ASVs to OTUs
~~~~~~~~~~~~~~~~~~~~~~~

| Clustering ASVs to OTUs with vsearch. 


.. code-block:: R
   :caption: get the size of ASVs
   :linenos:

    #!/usr/bin/env Rscript

    # specify input ASVs table and fasta
    ASV_table="ASV_table_TagJumpFiltered.txt" # specify ASV table file  
    ASV_fasta="ASVs_TagJumpFiltered.fasta"    # specify ASVs fasta file  

    ################################
    # remove "Sequence" column from the ASV table (if present)
    ASV_table = read.table(ASV_table, sep = "\t", check.names = F, 
                                header = T, row.names = 1)

    # check ASV table; if 1st col is sequence, then remove it
    if (colnames(ASV_table)[1] == "Sequence") {
        cat("## removing 'Sequence' column ... \n")
        ASV_table = ASV_table[, -1]
    }
    ################################
    # add size annototation to ASV seqs
    library(Biostrings)

    # add 'sum' column
    ASV_table$sum = rowSums(ASV_table)
    # make ASV_sums object
    ASV_sums = setNames(ASV_table$sum, rownames(ASV_table))
    # Read the FASTA file
    ASV_fasta = readDNAStringSet(ASV_fasta)
    # add ";size=*" to ASV_fasta
    names(ASV_fasta) = sapply(names(ASV_fasta), function(header) {
        paste0(header, ";size=", ASV_sums[header])
    })
    # write fasta file
    writeXStringSet(ASV_fasta, "ASVs.size.fasta",
                            width = max(width(ASV_fasta)))


.. code-block:: bash
   :caption: clustering with vsearch
   :linenos:

    #!/bin/bash 

    # specify the clustering threshold
    clustering_thresh="0.97"

    # make output dir.
    output_dir="OTU_table"
    mkdir -p $output_dir
    export output_dir

    ### cluster ASVs using vsearch.
    vsearch --cluster_fast ASVs.size.fasta \
        --id $clustering_thresh \
        --iddef 2 \
        --sizein \
        --xsize \
        --fasta_width 0 \
        --centroids $output_dir/OTUs.fasta \
        --uc $output_dir/OTUs.uc


.. code-block:: R
   :caption: generate an OTU table based on the clustered ASVs (.uc file).
   :linenos:

    #!/usr/bin/Rscript

    # specify input ASV table (the same one as for 'get the size of ASVs')
    ASV_table="ASV_table_TagJumpFiltered.txt"
    
    # read output dir
    output_dir = Sys.getenv('output_dir')

    # read output from vsearch clustering (-uc OTU.uc)
    inp_UC = file.path(output_dir, "OTUs.uc") 
    ################################
    library(data.table)
    # load input data - ASV table
    ASV_table = fread(file = ASV_table, header = TRUE, sep = "\t")

    # check ASV table; if 2nd col is 'Sequence', then remove it
    if (colnames(ASV_table)[2] == "Sequence") {
        cat("## removing 'Sequence' column ... \n")
        ASV_table = ASV_table[, -2]
    }

    ## Load input data - UC mapping file
    UC = fread(file = inp_UC, header = FALSE, sep = "\t")
    UC = UC[ V1 != "S" ]
    UC[, ASV := tstrsplit(V9, ";", keep = 1) ]
    UC[, OTU := tstrsplit(V10, ";", keep = 1) ]
    UC[V1 == "C", OTU := ASV ]
    UC = UC[, .(ASV, OTU)]

    # convert ASV table to long format
    ASV = melt(data = ASV_table,
        id.vars = colnames(ASV_table)[1],
        variable.name = "SampleID", value.name = "Abundance")
    ASV = ASV[ Abundance > 0 ]
    # add colnames, to make sure 1st is 'ASV'
    colnames(ASV) = c("ASV", "SampleID", "Abundance")

    # add OTU IDs
    ASV = merge(x = ASV, y = UC, by = "ASV", all.x = TRUE)
    # summarize
    OTU = ASV[ , .(Abundance = sum(Abundance, na.rm = TRUE)), 
                                by = c("SampleID", "OTU")]

    # reshape OTU table to wide format
    OTU_table = dcast(data = ASV,
        formula = OTU ~ SampleID,
        value.var = "Abundance",
        fun.aggregate = sum, fill = 0)

    # write OTU table
     # OTU names correspond to most abundant ASV in an OTU
    fwrite(x = OTU_table, file = file.path(output_dir, 
                                    "OTU_table.txt"), sep = "\t")

____________________________________________________

.. _postclustering16S:

Post-clustering
~~~~~~~~~~~~~~~

Post-cluster OTUs with LULU to merge consistently co-occurring 'daughter-OTUs'.

.. code-block:: bash
   :caption: generate match list for post-clustering
   :linenos:

    #!/bin/bash

    # go to directrory that contains OTUs
    cd $output_dir # 'OTU_table' in this case

    # make blast database for post-clustering
    makeblastdb -in OTUs.fasta -parse_seqids -dbtype nucl

    # generate match list for post-clustering
    blastn -db OTUs.fasta \
        -outfmt '6 qseqid sseqid pident' \
        -out match_list.txt \
        -qcov_hsp_perc 75 \
        -perc_identity 90 \
        -query OTUs.fasta \
        -num_threads 20


.. code-block:: R
   :caption: run LULU post-clustering
   :linenos:

    #!/usr/bin/Rscript

    # specify minimum threshold of sequence similarity considering any OTU as an error of another
    min_match = "90"

    # specify OTU table 
    OTU_table="OTU_table.txt"

    ################################
    library(devtools)
    # load OTU table and match list
    otutable = read.table(OTU_table, header = T, row.names = 1, sep = "\t")
    matchlist = read.table("match_list.txt")

    curated_result = lulu::lulu(otutable, matchlist, 
        minimum_match = min_match)

    # write post-clustered OTU table to file
    curated_table = curated_result$curated_table
    curated_table = cbind(OTU = rownames(curated_table), curated_table)
    write.table(curated_table, file ="OTU_table_LULU.txt", 
                sep = "\t", row.names = F, quote = FALSE)
    write.table(curated_result$discarded_otus, 
                file ="merged_units.lulu", col.names = FALSE, quote = FALSE)

.. note:: 

  Note that if the sample names start with a number, then the output OTU table may contain "X" prefix in the sample names. 


.. code-block:: bash
   :caption: match OTUs.fasta with post-clustered table (OTU_table_LULU)
   :linenos:

    #!/bin/bash

    # specify post-clustered table
    OTU_table="OTU_table_LULU.txt"
    # specify pre post-clustered OTUs fasta file
    OTUs_fasta="OTUs.fasta"

    # get matching OTUs
    awk 'NR>1{print $1}' $OTU_table > OTUs_LULU.list
    cat $OTUs_fasta | \
      seqkit grep -w 0 -f OTUs_LULU.list > OTUs_LULU.fasta

    # get matching taxonomy results for OTUs
    head -n 1 ../taxonomy.txt > taxonomy.txt
    cat ../taxonomy.txt | \
      grep -wf OTUs_LULU.list >> taxonomy.txt

    # remove unnecessary files
    rm OTUs.fasta.n*

    # move OTU_table two directories down
    cd ..
    mv $output_dir ../..

    
.. note:: 

    The final OTUs data is ``OTU_table_LULU.txt`` and ``OTUs_LULU.fasta`` in the ``OTU_table`` directory.
    The matching taxonomy file is ``taxonomy.txt`` in the ``OTU_table`` directory.

____________________________________________________

|logo_BGE_small| |eufund| |chfund| |ukrifund|
