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

| This is executable step-by-step pipeline for **COI** amplicon data analyses from Illumina sequencing machine.
|  
| The **full bioinformatics workflow can be automatically run through** `PipeCraft2 <https://pipecraft2-manual.readthedocs.io/en/latest/>`_ (v1.1.0; releasing this soon, with a tutorial),
| which implemets also various error handling processes and sequence summary statistics (lacking here in step-by-step code). 
| 
| The bioinformatic workflow results in amplicon sequence variants (ASVs) and well as operational taxonomic units (OTUs).


Dependencies
~~~~~~~~~~~~

+-----------------------------------------------------+----------+---------+
| Process                                             | Software | Version |
+=====================================================+==========+=========+
| :ref:`Remove primers <remove_primersCOI>`           | cutadapt | 4.9     |
+-----------------------------------------------------+----------+---------+
| :ref:`Quality filtering <quality_filteringCOI>`     | DADA2    | 1.30    |
+-----------------------------------------------------+----------+---------+
| :ref:`Denoise <denoiseCOI>`                         | DADA2    | 1.30    |
+-----------------------------------------------------+----------+---------+
| :ref:`Merge paired-end reads <denoiseCOI>`          | DADA2    | 1.30    |
+-----------------------------------------------------+----------+---------+
| :ref:`Chimera filtering <remove_chimerasCOI>`       | DADA2    | 1.30    |
+-----------------------------------------------------+----------+---------+
| :ref:`Remove tag-jumps <tagjumpsCOI>`               | UNCROSS2 |         |
+-----------------------------------------------------+----------+---------+
| :ref:`Merge sequencing runs* <mergeRunsCOI>`        | DADA2    | 1.30    |
+-----------------------------------------------------+----------+---------+
| :ref:`Taxonomy assignment <taxAssignCOI>`           | RDP      | 2.13    |
+-----------------------------------------------------+----------+---------+
| :ref:`Get target taxa <sorttaxaCOI>`                | R        |         |
+-----------------------------------------------------+----------+---------+
| :ref:`Remove NUMTs <numtsCOI>`                      | metaMATE | 0.4.3   |
+-----------------------------------------------------+----------+---------+
| :ref:`Clustering ASVs to OTUs <clusteringCOI>`      | vsearch  | 2.28.1  |
+-----------------------------------------------------+----------+---------+
| :ref:`Post-clusteringlustering <postclusteringCOI>` | LULU     | 0.1.0   |
+-----------------------------------------------------+----------+---------+

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
    maximum_error_rate="2" # Maximum error rate in primer string search;
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

____________________________________________________

.. _quality_filteringCOI:

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

____________________________________________________

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

____________________________________________________

.. _tagjumpsCOI:

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

____________________________________________________

.. _taxAssignCOI:

Taxonomy assignment
~~~~~~~~~~~~~~~~~~~

| Taxonomy assignment with the RDP-classifier against `CO1Classifier v5.1.0 database. <https://github.com/terrimporter/CO1Classifier>`_ 
| **---** `Download the CO1Classifier v5.1.0 for RDP here (click) <https://github.com/terrimporter/CO1Classifier/releases/download/RDP-COI-v5.1.0/RDP_COIv5.1.0.zip>`_ **---**

.. code-block:: bash
   :caption: assign taxonomy with RDP-classifier
   :linenos:

    #!/bin/bash

    # download the CO1Classifier reference databse
    wget \
      "https://github.com/terrimporter/CO1Classifier/releases/download/RDP-COI-v5.1.0/RDP_COIv5.1.0.zip"
    # unzip the database and edit name
    unzip RDP_COIv5.1.0.zip && mv mydata CO1Classifier_v5.1.0_RDP
    
    # specify reference database for RDP
    reference_database="CO1Classifier_v5.1.0_RDP/rRNAClassifier.properties"
    reference_database=$(realpath $reference_database) # get database names with full path

    # specify input fasta file
    cd ASV_table
    ASV_fasta="ASVs_TagJumpFiltered.fasta"

    # Run RDP-classifier
    time rdp_classifier \
            -Xmx12g \
            classify \
            -t $reference_database \
            -f allrank \
            -o RDP.taxonomy.txt \
            -q $ASV_fasta

____________________________________________________

.. _sorttaxaCOI:

Get target taxa
~~~~~~~~~~~~~~~

| This part filters the ASV dataset to include only target taxonomic group for the following analyses. 
| For example, if you are interested in Hymenoptera, then discard all ASVs that do not match to the target taxon based on the user defined threshold (default = 0.8). 
|
| :ref:`See below if no pre-selection is preferred <donotsorttaxaCOI>`

.. code-block:: R
   :caption: get only target taxon annotations
   :linenos:

    #!/usr/bin/env Rscript
    ### Filter dataset based on RDP classifier results to include target taxa 

    # specify taxon and threshold
    taxon="Metazoa"  # target taxonomic group(s); 
                         # multiple groups should be from the same taxonomic level
                         # separator is "," (e.g., "Hymenoptera, Lepidoptera")
    tax_level="kingdom"  # allowed levels: kingdom | phylum | class | order | family | genus
    threshold="0.8"      # threshold for considering an ASV as a target taxon

    # specify the ASV table and ASVs.fasta file that would be filtered to include only target taxa 
    ASV_fasta = "ASVs_TagJumpFiltered.fasta"
    ASV_table = "ASV_table_TagJumpFiltered.txt"

    # specify the RDP-classifier output file (taxonomy file)
    taxtab="RDP.taxonomy.txt"
    
    #--------------------------------------#
    library(stringr)
    library(dplyr)

    # read ASV table
    table = read.table(ASV_table, sep = "\t", check.names = F, header = T, row.names = 1)
    
    # read taxonomy table
    tax = read.table(taxtab, sep = "\t", check.names = F, row.names = 1)
    cat("\n Input =", nrow(tax), "features.\n")
    # remove not needed columns from tax dataframe
    tax = tax[, -c(1, 2, 3, 4, 6, 9, 12, 15, 18, 21, 24, 27)]
    # assign colnames for tax
    colnames(tax) = c("superkingdom", "superkingdom_BootS",
                    "kingdom", "kingdom_BootS",
                    "phylum","phylum_BootS",
                    "class", "class_BootS",
                    "order", "order_BootS",
                    "family", "family_BootS",
                    "genus", "genus_BootS",
                    "species", "species_BootS")

    # taxon list
    taxon_list = strsplit(taxon, ", ")[[1]]
 
    ### extract only target-taxon ASVs from the 'raw' RDP results
    tax_filtered = tax %>%
        filter(.data[[tax_level]] %in% taxon_list)

    ### change all tax ranks to "unclassified_*" when 
        # the bootstrap values is less than the specified threshold
    #kingdom
    tax_filtered = tax_filtered %>% mutate(kingdom = ifelse(kingdom_BootS < 
        threshold, paste0("unclassified_", superkingdom), as.character(kingdom)))
    #phylum
    tax_filtered = tax_filtered %>% mutate(phylum = ifelse(phylum_BootS < 
        threshold, paste0("unclassified_", kingdom), as.character(phylum)))
    #replace potential "unclassified_unclassified_" with "unclassified_"
    tax_filtered$class = stringr::str_replace(tax_filtered$class, "unclassified_unclassified_", 
                                                                            "unclassified_")
    #class
    tax_filtered = tax_filtered %>% mutate(class = ifelse(class_BootS < 
        threshold, paste0("unclassified_", phylum), as.character(class)))
    #replace potential "unclassified_unclassified_" with "unclassified_"
    tax_filtered$class = stringr::str_replace(tax_filtered$class, "unclassified_unclassified_", 
                                                                            "unclassified_")
    #order
    tax_filtered = tax_filtered %>% mutate(order = ifelse(order_BootS < 
        threshold, paste0("unclassified_", class), as.character(order)))
    #replace potential "unclassified_unclassified_" with "unclassified_"
    tax_filtered$order = stringr::str_replace(tax_filtered$order, "unclassified_unclassified_", 
                                                                            "unclassified_")
    #family
    tax_filtered = tax_filtered %>% mutate(family = ifelse(family_BootS < 
        threshold, paste0("unclassified_", order), as.character(family)))
    #replace potential "unclassified_unclassified_" with "unclassified_"
    tax_filtered$family = stringr::str_replace(tax_filtered$family, "unclassified_unclassified_", 
                                                                            "unclassified_")
    #genus
    tax_filtered = tax_filtered %>% mutate(genus = ifelse(genus_BootS < 
        threshold, paste0("unclassified_", family), as.character(genus)))
    #replace potential "unclassified_unclassified_" with "unclassified_"
    tax_filtered$genus = stringr::str_replace(tax_filtered$genus, "unclassified_unclassified_", 
                                                                            "unclassified_")

    # species to genus_sp when the bootstrap values is < 0.9
    tax_filtered = tax_filtered %>% mutate(species = ifelse(species_BootS < 0.9, 
                                                        paste0(genus, "_sp"), species))
   
    ### count occurrences of each taxon in df (RDP results)
    count_taxa = function(df, taxa) {
    sapply(taxa, function(taxon) sum(apply(df, 1, function(row) any(row == taxon))))
    }
    taxon_counts = count_taxa(tax_filtered, taxon_list)

    # Check the counts
    if (all(taxon_counts == 0)) {
        print("ERROR: None of the specified taxa are present in the RDP results.")
    } else {
        if (any(taxon_counts == 0)) {
        warning("One or more of the specified taxa are not present in the RDP results.")
        }
        print(taxon_counts)
    }

    ### extract only target-taxon ASVs from the 'threshold filtered' RDP results
    tax_filtered_thresh = tax_filtered %>%
        filter(.data[[tax_level]] %in% taxon_list)
    # write filtered RDP taxonomy table
    tax_filtered_thresh = cbind(ASV = rownames(tax_filtered_thresh), tax_filtered_thresh)
    write.table(tax_filtered_thresh, 
                file = "RDP.taxonomy.filt.txt",  
                quote = F, 
    	        row.names = F,
                sep = "\t")
    
    ### filter the ASV table to match ASVs that were kept in the tax_filtered table
    table_filt = table[rownames(table) %in% rownames(tax_filtered_thresh), ]

    ### check ASV table; if 1st col is sequence, then remove it for metaMATE
    if (colnames(table_filt)[1] == "Sequence") {
        cat(";; 2nd column was 'Sequence', removing this ... \n")
        table_filt = table_filt[, -1]
    }

    # write filtered table
    table_filt = cbind(ASV = rownames(table_filt), table_filt)
    write.table(table_filt, 
                file = paste0(sub("\\.[^.]*$", "_tax_filt.txt", ASV_table)),  
                quote = F, 
    	        row.names = F,
                sep = "\t")

    # filter ASV_fasta
    library(Biostrings)
    fasta = readDNAStringSet(ASV_fasta)
    fasta.tax_filt = fasta[names(fasta) %in% rownames(table_filt)]
    # write filtered ASV_fasta
    writeXStringSet(fasta.tax_filt, 
                    paste0(sub("\\.[^.]*$", "_tax_filt.fasta", ASV_fasta)), 
                    width = max(width(fasta.tax_filt)))

    
____________________________________________________

.. _donotsorttaxaCOI:

**If no pre-selection is preferred, then just remove "Sequence" column from the ASV table**

.. code-block:: R
   :caption: remove "Sequence" column from the ASV table
   :linenos:

    # read ASV table
    ASV_table = "ASV_table.txt"
    table = read.table(ASV_table, sep = "\t", check.names = F, header = T, row.names = 1)

    # check ASV table; if 1st col is sequence, then remove it for metaMATE
    if (colnames(table)[1] == "Sequence") {
        cat("## removing 'Sequence' column ... \n")
        table = table[, -1]

        # write filtered table
        table_filt = cbind(ASV = rownames(table), table)
        write.table(table_filt, 
                file = paste0(sub("\\.[^.]*$", ".noSeq.txt", ASV_table)),  
                quote = F, 
    	        row.names = F,
                sep = "\t")

    } else {
        cat("## there was no 'Sequence' column; proceed with the current table ... \n")
    }

____________________________________________________    
  
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
    # specify specifications file for metaMATE
    specifications="specifications.txt"
    specifications=$(realpath $specifications) # get full directory path


    # download the CO1Classifier reference databse
    wget \
     "https://github.com/terrimporter/CO1Classifier/releases/download/SINTAX-COI-v5.1.0-ref/SINTAX_COIv5.1.0_ref.zip"
    # unzip the database and edit name
    unzip SINTAX_COIv5.1.0_ref.zip && mv training CO1Classifier_v5.1.0 
    mv CO1Classifier_v5.1.0/sintax.fasta CO1Classifier_v5.1.0/CO1Classifier_v5.1.0.fasta
    
    # specify reference database for metaMATE
    reference_database="CO1Classifier_v5.1.0/CO1Classifier_v5.1.0.fasta"
    reference_database=$(realpath $reference_database) # get full directory path


.. code-block:: bash
   :caption: run metaMATE-find
   :linenos:

    #!/bin/bash
    ## run metaMATE-find
  
    ## go to the directory that hosts your ASVs.fasta and ASV table files.
  
    # specify input ASVs table and fasta
    ASV_table="ASV_table_TagJumpFiltered_tax_filt.txt" # make sure that the 2nd col is not "Sequence"
    ASV_fasta="ASVs_TagJumpFiltered_tax_filt.fasta"    # specify ASVs fasta file 

    # specify variables
    genetic_code="5"        # the standard genetic code. 5 is invertebrate mitochondrial code
    length="313"            # the expected length of an amplicon
    basesvariation="9"      # allowed length variation (bp) from the expected length of an amplicon
    taxgroups="undefined"   # (optional); if sequence binning is to be performed on 
                               # a per-taxon basis (as in specifications file) 
                               # then specify the taxon grouping file
    NA_abund_thresh="0.05"  # nonauthentic_retained_estimate_p < 0.05 (value from mateMATE results);
                               # the allowed abundance threshold of 
                               # non-validated OTUs/ASVs in the filtered dataset.
    abundance_filt="FALSE"  # TRUE/FALSE; if FALSE, then NA_abund_thresh is ineffective, 
                               # and no filtering is done based on the ASV abundances,
                               # i.e., filter only based on length, basesvariation and genetic_code.
                               # FALSE may be used when the seq-depth for the target taxa is low.
                               # If TRUE, then NA_abund_thresh will be applied. 
                               
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
    # remove old $output_dir if exists
    if [[ -d $output_dir ]]; then
        rm -rf $output_dir
    fi

    # if perfoming clade binning, then WARNING when processing more than 65,536 ASVs
    ASVcount=$(grep -c "^>" $ASV_fasta)
    if (( $ASVcount > 65536 )); then
        printf '%s\n' "WARNING]: clade binning NOT performed, 
         because the input ASVs limit is 65,536 for that.
         Current input has $ASVcount ASVs."
    fi

    # check abundance_filt; 
     # if FALSE then make new specifications file, that excludes abundance filtering
    if [[ $abundance_filt == "FALSE" ]]; then
        printf '%s\n' "[library; n; 0-1/2]" > specifications0.txt
        specifications=$(realpath specifications0.txt)
    fi

    # quick check of the specifications file, has to contain "library" | "total" | "clade" | "taxon"
    if ! grep -q -e "library" -e "total" -e "clade" -e "taxon" $specifications; then
        printf '%s\n' "ERROR]: specifications file seems to be wrong. 
         Does not contain any of the terms (library, total, clade, taxon)."
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
        --table $genetic_code \
        --threads 8 \
        --output $output_dir \
        --overwrite $taxgroups

    # check for the presence of "metamate_out" dir and "resultcache" file (did metaMATE-find finish)
    if [[ -d $output_dir ]] && [[ -e $output_dir/resultcache ]] && [[ -e $output_dir/results.csv ]]; then
        printf '\n%s\n\n' "metaMATE-find finished, proceed"
        # export variables for below script (Rscript)
        if [[ $abundance_filt != "FALSE" ]]; then
            printf '%s\n' "exporting NA_abund_thresh of $NA_abund_thresh for metaMATE-dump"
            export NA_abund_thresh
        else
            export output_dir
            export abundance_filt
        fi
    else 
        printf '%s\n' "ERROR]: cannot find the $output_dir (metaMATE-find output) 
         to start metaMATE-dump OR no authentic ASVs found??"
    fi


.. code-block:: bash
   :caption: get the results_index from the metamate_out/results.csv file
   :linenos:

    #!/usr/bin/env Rscript
  
    ## read results.csv
    output_dir = Sys.getenv('output_dir') # = "metamate_out" as specified above
    find_results = read.csv(file.path(output_dir, "results.csv"))

    # get variables
    abundance_filt = Sys.getenv('abundance_filt')

    ## filter results if abundance_filt is FALSE
    if (abundance_filt == "FALSE"){
        result_index = "0" # get first result_index (library_n = 0) 
        write(result_index, file.path(output_dir, "selected_result_index.txt"))
    }

    ## filter results based on NA_abund_thresh 
    if (abundance_filt != "FALSE"){
        NA_abund_thresh = as.numeric(Sys.getenv('NA_abund_thresh'))
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
    }



.. code-block:: bash
   :caption: run metaMATE-dump to discard putative artefact ASVs
   :linenos:

    #!/bin/bash

    ## metaMATE-dump 
    ASV_fasta=$(basename $ASV_fasta)
    
    # read result_index
    read -r result_index < $output_dir/selected_result_index.txt
    printf '%s\n' " - selected result_index = $result_index"

    # run metaMATE-dump
    printf '%s\n' "# Running metaMATE-dump"
    metamate dump \
    --asvs $ASV_fasta \
    --resultcache $output_dir/resultcache \
    --output $output_dir/${ASV_fasta%.*}_metaMATE.filt \
    --overwrite \
    --resultindex $result_index

    # generate a list of ASV IDs 
    seqkit seq -n $output_dir/${ASV_fasta%.*}_metaMATE.filt.fasta > \
                        $output_dir/${ASV_fasta%.*}_metaMATE.filt.list

    # filter the ASV table; include only the ASVs that are in ${ASV_fasta%.*}_metaMATE.filt.list
    awk -v var="$output_dir/${ASV_fasta%.*}" 'NR==1; NR>1 {print $0 | \
            "grep -Fwf "var"_metaMATE.filt.list"}' $ASV_table > \
             $output_dir/${ASV_table%.*}_metaMATE.filt.txt

    # filter the RDP.taxonomy.filt.txt file to include only ASVs retained by metaMATE
    awk -v var="$output_dir/${ASV_fasta%.*}" 'NR==1; NR>1 {print $0 | \
            "grep -Fwf "var"_metaMATE.filt.list"}' RDP.taxonomy.filt.txt > \
            $output_dir/RDP.taxonomy.metaMATE.filt.txt
                                                                              
                                                                              
    # write discarded ASVs list
    comm -3 <(sort <(seqkit seq -n $ASV_fasta)) \
            <(sort $output_dir/${ASV_fasta%.*}_metaMATE.filt.list) \
            > $output_dir/metaMATE.discarded.list
    
    
.. code-block:: bash
   :caption: optionally rescue discarded ASVs that have GENUS level bootstrap value > 0.9
   :linenos:

    #!/bin/bash
    
    # get discarded ASVs (RDP taxonomy list)
    grep -Fwf $output_dir/metaMATE.discarded.list RDP.taxonomy.filt.txt \
    > $output_dir/metaMATE.discarded.RDP.taxonomy.txt

    # get the rescued ASVs that have GENUS level bootstrap value >= 0.9
    awk -F'\t' '$26 >= 0.9' $output_dir/metaMATE.discarded.RDP.taxonomy.txt \
    > $output_dir/rescued.txt

    # check if rescued.txt exists and is not empty
    if [[ -s $output_dir/rescued.txt ]]; then
        # add the rescued ASVs to $output_dir/RDP.taxonomy.metaMATE.filt.txt
        cat $output_dir/rescued.txt >> $output_dir/RDP.taxonomy.metaMATE.filt.txt

        # add the rescued ASVs to $output_dir/${ASV_fasta%.*}_metaMATE.filt.fasta
        seqkit grep -w 0 -f <(awk -F'\t' '{print $1}' $output_dir/rescued.txt) $ASV_fasta \
            >> $output_dir/${ASV_fasta%.*}_metaMATE.filt.fasta

        # add the rescued ASVs to $output_dir/${ASV_table%.*}_metaMATE.filt.txt
        grep -wf <(awk -F'\t' '{print $1}' $output_dir/rescued.txt) $ASV_table \
            >> $output_dir/${ASV_table%.*}_metaMATE.filt.txt

        printf '%s\n' "Rescued $(wc -l < $output_dir/rescued.txt) ASVs"
    else
        printf '%s\n' "No ASVs to rescue"
    fi


.. note:: 

    Herein case, the final filtered data is ``ASV_table_tax_filt_metaMATE.filt.txt`` and ``ASVs_tax_filt_metaMATE.filt.fasta`` in the ``metamate_out`` directory.
    The filtered RDP-classifier results (matching the ASVs in the latter files) is ``RDP.taxonomy.metaMATE.filt.txt`` in the ``metamate_out`` dir.
    
    If deemed relevant, then you may proceed with the below workflow below that includes clustering ASVs to OTUs. 

____________________________________________________

.. _clusteringCOI:

Clustering ASVs to OTUs
~~~~~~~~~~~~~~~~~~~~~~~

| Clustering ASVs to OTUs with vsearch. 

.. code-block:: R
   :caption: get the size of ASVs
   :linenos:

    #!/usr/bin/env Rscript

    # specify input ASVs table and fasta
    ASV_table="ASV_table_tax_filt_metaMATE.filt.txt" # specify ASV table file  
    ASV_fasta="ASVs_tax_filt_metaMATE.filt.fasta"    # specify ASVs fasta file  

    ################################
    library(Biostrings)
    # Read the ASV table
    ASV_table = read.table(ASV_table, sep = "\t", check.names = F, 
                                header = T, row.names = 1)

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
    ASV_table="ASV_table_tax_filt_metaMATE.filt.txt"
    
    # read output dir
    output_dir = Sys.getenv('output_dir')

    # read output from vsearch clustering (-uc OTU.uc)
    inp_UC = file.path(output_dir, "OTUs.uc") 
    ################################
    library(data.table)
    # load input data - ASV table
    ASV_table = fread(file = ASV_table, header = TRUE, sep = "\t")

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


.. _postclusteringCOI:

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

    # get matching RDP taxonomy results
    head -n 1 ../RDP.taxonomy.metaMATE.filt.txt > RDP.taxonomy.txt
    cat ../RDP.taxonomy.metaMATE.filt.txt | \
      grep -wf OTUs_LULU.list >> RDP.taxonomy.txt

    # remove unnecessary files
    rm OTUs.fasta.n*

    # move OTU_table two directories down
    cd ..
    mv $output_dir ../..

    
.. note:: 

    The final OTUs data is ``OTU_table_LULU.txt`` and ``OTUs_LULU.fasta`` in the ``OTU_table`` directory.

    The matching RDP taxonomy files are ``RDP.taxonomy.txt`` in the ``OTU_table`` directory.

____________________________________________________

|eufund| |chfund| |ukrifund|
