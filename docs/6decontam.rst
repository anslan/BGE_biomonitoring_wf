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

.. _remove_contaminants:

Identify and remove contaminants
********************************

Here, we apply the `decontam <https://github.com/benjjneb/decontam>`_ R package to identify and remove contaminants from the metabarcoding feature table
(based on the `decontam vignette <https://benjjneb.github.io/decontam/vignettes/decontam_intro.html>`_).

Alongside with the biological samples, the metabarcoding workflow needs to include also **control samples to assess the potential contaminations** 
due to the sensitivity of the DNA-based methods. 
The control sample can include field controls (e.g., a sample that should contain no target DNA and is exposed to the same field conditions and handling as real samples),
blank DNA extractions (i.e., no sample substrate added to DNA extraction tubes) and 
blank PCR reactions (i.e., no reaction with template DNA).

In many cases, even if the control samples appear to be clean, they may still contain some sequences in the final feature table.
For smaller projects, including only few samples, the deconamintion process may be easy via visual inspections. 
However, **for larger projects, including many samples, the decontamination process may be more challenging, and here tools such as decontam can be very helpful**.

___________________________________________________

Input data
==========

The input data for the decontamination process include first **metabarcoding feature table** (with control samples) and **sample metadata**.
Then, later, the identified contaminanted features can be removed also from the :ref:`taxonomy table and fasta file <taxonomy_table>`.

- The metabarcoding feature table (per sequencing run/per batch processed simultaneously) with sequence abundances (example):

+------------+---------+---------+---------+--------------------+---------------------+-------------------+
|            | sample1 | sample2 | sample3 | :red:`s_ExtrBlank` | :red:`s_FiledBlank` | :red:`s_pcrBlank` |
+------------+---------+---------+---------+--------------------+---------------------+-------------------+
| **OTU_01** | 4290    | 40550   | 3902    | 0                  | 0                   | 0                 |
+------------+---------+---------+---------+--------------------+---------------------+-------------------+
| **OTU_02** | 550     | 34501   | 2       | 0                  | 0                   | 0                 |
+------------+---------+---------+---------+--------------------+---------------------+-------------------+
| **OTU_03** | 0       | 0       | 0       | 136                | 7                   | 201               |
+------------+---------+---------+---------+--------------------+---------------------+-------------------+
| **OTU_04** | 3061    | 0       | 3465    | 1                  | 0                   | 0                 |
+------------+---------+---------+---------+--------------------+---------------------+-------------------+
| ...        | ...     | ...     | ...     | ...                | ...                 | ...               |
+------------+---------+---------+---------+--------------------+---------------------+-------------------+


.. admonition:: control samples

    Samples **s_ExtrBlank**, **s_FiledBlank** and **s_pcrBlank** are control samples for DNA extraction, field blank and PCR blank, respectively.
    They are used to assess the potential contaminations. 

In the below script, this table is represented as **tab-delimited** ``OTU_table.txt`` file.


- The sample metadata (example): 

+------------------+-------------------+----------------+------------+
|                  | DNA_concentration | sample_type    | seq_counts |
+------------------+-------------------+----------------+------------+
| **sample1**      | 19.3              | true_sample    | 133902     |
+------------------+-------------------+----------------+------------+
| **sample2**      | 20.1              | true_sample    | 155423     |
+------------------+-------------------+----------------+------------+
| **sample3**      | 20                | true_sample    | 250432     |
+------------------+-------------------+----------------+------------+
| **s_ExtrBlank**  | 0.12              | control_sample | 247        |
+------------------+-------------------+----------------+------------+
| **s_FiledBlank** | 1.3               | control_sample | 34         |
+------------------+-------------------+----------------+------------+
| **s_pcrBlank**   | 0.19              | control_sample | 325        |
+------------------+-------------------+----------------+------------+
| ...              | ...               | ...            | ...        |
+------------------+-------------------+----------------+------------+

The ``DNA_concentration`` is **optional** and can be included if available. This represents the DNA concentration of the sample 
prior pooling the samples into a single library for sequencing.
The ``sample_type`` indicates the type of the sample: true_sample or control_sample. 
The ``seq_counts`` represents the number of sequences assigned to the sample. 


.. note::

    ``sample_type`` and ``seq_counts`` columns can be generated from the metabarcoding feature table (see below). 
    **So, essentially, no specific metadata file is required prior** *decontam* **if DNA concentration is not available.**


___________________________________________________

Decontamination process
=======================

The decontamination process is performed using the *decontam* R package.
Here, we do not use DNA concentration data, but only sequence counts and sample type (true_sample or control_sample). 
**The procedure follows "prevalence" method:**
the presence/absence of each feature in true samples is compared to the presence/absence in control samples to identify contaminants.
More information about that here: `decontam vignette <https://benjjneb.github.io/decontam/vignettes/decontam_intro.html#identify-contaminants---prevalence>`_ 
and in the `decontam paper <https://doi.org/10.1186/s40168-018-0605-2>`_. 


.. code-block:: R
   :caption: decontamination
   :linenos:

   #!/usr/bin/Rscript

   # install decontam package if not already installed 
   if (!requireNamespace("decontam", quietly = TRUE)) {
   BiocManager::install("decontam")
   } else {
   cat("decontam is already installed, version", 
        as.character(packageVersion("decontam")), "\n")
   }

   # load decontam package
   library(decontam)
   set.seed(1)

   # read OTU table 
   # tab-delimted txt file, samples as columns, OTUs in rows
   OTU_table = read.table("OTU_table.txt", 
                        header = TRUE,        
                        sep = "\t",
                        row.names = 1)

   #------------------------#
   ### Make metadata file ###
   #------------------------#
   # add seq_counts columns per sample
   seq_counts <- colSums(OTU_table)

   metadata <- data.frame(
   sample_id = names(seq_counts),
   seq_counts = seq_counts
   )

   # add sample_type column
   metadata$sample_type <- ifelse(
   grepl("^s_(ExtrBlank|FiledBlank|pcrBlank)", metadata$sample_id),
   "control_sample",
   "true_sample"
   )

   # sort by seq_counts
   metadata <- metadata[order(metadata$seq_counts),]
   # add Index column
   metadata$Index <- seq(nrow(metadata))

   # check metadata
   head(metadata)


   ### check summary statistics
   summary_stats <- metadata %>%
   group_by(sample_type) %>%
   summarise(
        n_samples = n(),
        total_sequences = sum(seq_counts),
        mean_sequences = mean(seq_counts),
        .groups = 'drop'
   )
   print(summary_stats)

   #-----------------------#
   ### find contaminants ###
   #-----------------------#
   # ensure sample order matches
   metadata <- metadata[match(colnames(OTU_table), metadata$sample_id), ]

   # define control samples for decontam 
   is.neg <- metadata$sample_type == "control_sample"

   # run decontam with default settings
   # note that t(OTU_table) is used to transpose the OTU table for decontam
   contam_result <- isContaminant(t(OTU_table), 
                                method = "prevalence", 
                                threshold = 0.1,
                                neg = is.neg)

   # check how many OTUs are considered as contaminants
   table(contam_result$contaminant)

   #------------------------#
   ### check contaminants ###
   #------------------------#
   # transform OTU table presence/absence data
   otu_pa <- OTU_table
   otu_pa[otu_pa > 0] <- 1

   # make dataframe for plotting
   pos_sample_ids <- na.omit(metadata$sample_id[metadata$sample_type == "true_sample"])
   neg_sample_ids <- na.omit(metadata$sample_id[metadata$sample_type == "control_sample"])
   otu_pa_pos <- otu_pa[, pos_sample_ids]
   otu_pa_neg <- otu_pa[, neg_sample_ids]
   pa_neg <- rowSums(otu_pa_neg)
   pa_pos <- rowSums(otu_pa_pos)
   df.pa <- data.frame(
   pa.neg = pa_neg,
   pa.pos = pa_pos,
   contaminant = contam_result$contaminant
   )

   # visualize contaminant prevalence
   ggplot(data = df.pa, aes(x = pa.neg, 
                            y = pa.pos, 
                            color = contaminant)) +
   geom_point() +
   xlab("Prevalence of negative controls") +
   ylab("Prevalence of true samples") +
   theme_minimal()

   #-------------------------#
   ### remove contaminants ###
   #-------------------------#
   # get contaminant OTU names
   contaminants <- rownames(contam_result)[contam_result$contaminant == TRUE]

   # filter out contaminant OTUs
   OTU_table_decontam <- OTU_table[!rownames(OTU_table) %in% contaminants, ]

   #----------------------------#
   ### remove control samples ###
   #----------------------------#
   # Get true sample IDs
   true_samples <- metadata %>%
   filter(sample_type == "true_sample") %>%
   pull(sample_id)

   # filter OTU table to include only true_samples
   OTU_table_decontam <- OTU_table_decontam[, colnames(OTU_table_decontam) %in% true_samples]

   # Remove OTUs with zero sequences (if any)
   OTU_table_decontam <- OTU_table_decontam[rowSums(OTU_table_decontam) > 0, ]

   # Summary 
   initial_number_of_OTUs = nrow(OTU_table)
   number_of_OTUs_after_excluding_controls = nrow(OTU_table_decontam)
   cat("Removed", initial_number_of_OTUs-number_of_OTUs_after_excluding_controls, 
        "contaminant OTUs\n")
   cat("OTUs left in the table:", number_of_OTUs_after_excluding_controls)

   ### Write decontaminanted OTU table to Excel file
   library(openxlsx)
   write.xlsx(OTU_table_decontam, 
            "OTU_table_decontam.xlsx",
            rowNames = TRUE)

   # write the list of contaminant OTUs to a file;
   # this will be used to filter also the taxonomy table
   write.xlsx(as.data.frame(contaminants), 
            "contaminants.xlsx",
            rowNames = FALSE)

**The decontaminated OTU table is written to an Excel file (in your output directory).**

.. _taxonomy_table:

___________________________________________________

Now, the OTUs that were considered as **contaminants should be removed also from the taxonomy table** and **fasta file**.


- Example of an input taxonomy table:

+------------+----------------+---------------+------------------+---------------------+
|            | Order          | Family        | Genus            | species             |
+------------+----------------+---------------+------------------+---------------------+
| **OTU_01** | Lepidoptera    | Nymphalidae   | Aglais           | Aglais urticae      |
+------------+----------------+---------------+------------------+---------------------+
| **OTU_02** | Lepidoptera    | Nymphalidae   | Genus_0022       | Genus_0022_sp       |
+------------+----------------+---------------+------------------+---------------------+
| **OTU_03** | Carnivora      | Felidae       | Felis            | Felis catus         |
+------------+----------------+---------------+------------------+---------------------+
| **OTU_04** | Sarcoptiformes | Pyroglyphidae | Dermatophagoides | Dermatophagoides_sp |
+------------+----------------+---------------+------------------+---------------------+
| ...        | ...            | ...           | ...              | ...                 |
+------------+----------------+---------------+------------------+---------------------+

In the below script, this table is represented as **tab-delimited** ``taxonomy_table.txt`` file.

- Example of an input fasta file:

.. admonition:: fasta file 

    .. code-block:: text
       :class: small-font

       >OTU_01
       ACTTTATATTTTATTTTTGGAATTTGAGCAGGAATAGTAGGAACTTCTCTTAGTTTATTAATTCG...
       >OTU_02
       ATAGTAGGAACATCTCTTAGTTTATTAATTCGAACTGAACTAGGAAATCCAGGTTCACTTATTGG...
       >OTU_03
       TCTTTACCTTTTATTCGGTGCCTGAGCTGGCATGGTGGGGACTGCTCTTAGTCTTCTAATCCGGG...
       >OTU_04
       TACTTTGTATTTTGTTTTTGGGGTGTGATCTGGTATGTTGGGGACTAGGTTCAGAAGACTAATTC...
       ...

In the below script, this file is represented as ``OTUs.fasta``.


.. code-block:: R 
   :caption: remove contaminants from the taxonomy table
   :linenos:

   #!/usr/bin/Rscript

   # load taxonomy table 
   taxonomy = read.table("taxonomy_table.txt", 
                        header = TRUE,
                        sep = "\t",
                        row.names = 1)

   # load list of contaminant OTUs 
   library(openxlsx)
   contaminants = read.xlsx("contaminants.xlsx",
                            rowNames = FALSE,
                            colNames = TRUE)

   # extract contaminant OTU names (assuming first column is header)
   contaminants = contaminants[, 1]

   #---------------------------#
   ### Filter taxonomy table ###
   #---------------------------#
   # filter out contaminant OTUs from taxonomy table
   taxonomy_filtered = taxonomy[!rownames(taxonomy) %in% contaminants, ]

   cat("Taxonomy table - Initial OTUs:", nrow(taxonomy), "\n")
   cat("Taxonomy table - After filtering:", nrow(taxonomy_filtered), "\n")
   cat("Removed", nrow(taxonomy) - nrow(taxonomy_filtered), 
        "contaminant OTUs from taxonomy table\n")

   # write filtered taxonomy table
   write.table(taxonomy_filtered,
                "taxonomy_table_decontam.txt",
                sep = "\t",
                quote = FALSE,
                row.names = TRUE,
                col.names = NA)

   #-----------------------#
   ### Filter fasta file ###
   #-----------------------#
   library(seqinr)

   # load fasta file
   OTUs_fasta = seqinr::read.fasta("OTUs.fasta")

   # filter out contaminant OTUs from fasta file
   OTUs_fasta_filtered = OTUs_fasta[!names(OTUs_fasta) %in% contaminants]

   cat("\nFasta file - Initial sequences:", length(OTUs_fasta), "\n")
   cat("Fasta file - After filtering:", length(OTUs_fasta_filtered), "\n")
   cat("Removed", length(OTUs_fasta) - length(OTUs_fasta_filtered), 
                  "contaminant sequences from fasta file\n")

   # write filtered fasta file
   write.fasta(sequences = OTUs_fasta_filtered, 
                names = names(OTUs_fasta_filtered),
                file.out = "OTUs_decontam.fasta",
                nbchar = 999)


.. admonition:: output

    The final output data is ``OTU_table_decontam.txt``, 
    ``taxonomy_table_decontam.txt`` and ``OTUs_decontam.fasta``. 
    Additionally, the list of contaminant OTUs is written to ``contaminants.xlsx``.


____________________________________________________

|logo_BGE_small| |eufund| |chfund| |ukrifund|
