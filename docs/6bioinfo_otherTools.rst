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
            font-size: 20px; /* text size */
            padding: 5px;   /* add  padding */
        }
    </style>

.. role:: yellow-background


|logo_BGE_alpha|


Other tools
***********


BLAST
~~~~~

| Assign taxonomy with `BLAST <https://pubmed.ncbi.nlm.nih.gov/2231712/>`_. 
| Reference database is a fasta file; no special formatting is required.

.. code-block:: bash
   :caption: BLAST
   :linenos:

    #!/bin/bash

    # specify the query fasta file
    fasta=$"ASVs.fasta"

    # specify reference database for BLAST 
    reference_database="../CO1Classifier_v5.1.0/CO1Classifier_v5.1.0.fasta"
    reference_database=$(realpath $reference_database) # get full directory path

    ## if the database is just in fasta format, then convert it to BLAST format

    ### Check and assign BLAST database
    d1=$(echo $reference_database | awk 'BEGIN{FS=OFS="."}{print $NF}')
    #make blast database if db is not formatted for BLAST
    db_dir=$(dirname $reference_database)
    check_db_presence=$(ls -1 $db_dir/*.nhr 2>/dev/null | wc -l)
    if (( $check_db_presence != 0 )); then
        if [[ $d1 == "fasta" ]] || [[ $d1 == "fa" ]] || \
            [[ $d1 == "fas" ]] || [[ $d1 == "fna" ]] || \
            [[ $d1 == "ffn" ]]; then
            database=$"-db $reference_database"
        elif [[ $d1 == "ndb" ]] || [[ $d1 == "nhr" ]] || \
             [[ $d1 == "nin" ]] || [[ $d1 == "not" ]] || \
             [[ $d1 == "nsq" ]] || [[ $d1 == "ntf" ]] || \
             [[ $d1 == "nto" ]]; then
            reference_database=$(echo $reference_database | awk 'BEGIN{FS=OFS="."}NF{NF-=1};1')
            database=$"-db $reference_database"
        fi
    elif [[ $d1 == "fasta" ]] || [[ $d1 == "fa" ]] || \
         [[ $d1 == "fas" ]] || [[ $d1 == "fna" ]] || \
         [[ $d1 == "ffn" ]]; then
            printf '%s\n' "Note: converting fasta formatted database for BLAST"
            makeblastdb -in $reference_database -input_type fasta -dbtype nucl
            database=$"-db $reference_database"
    fi

    #BLAST
    printf '%s\n' "# Running BLAST for $(grep -c "^>" $fasta) sequences"
    blastn -strand plus \
                -num_threads 20 \
                -query $fasta \
                $database \
                -out 10BestHits.txt -task blastn \
                -max_target_seqs 10 -evalue=0.001 \
                -word_size=7 -reward=1 \
                -penalty=-1 -gapopen=1 -gapextend=2 \
    -outfmt "6 qseqid stitle qlen slen qstart qend sstart send evalue length nident mismatch gapopen gaps sstrand qcovs pident"
    
    #qseqid = Query Seq-id
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
    #qcovs = Query Coverage Per Subject

    ### parse BLAST 1st hit 
    awk 'BEGIN{FS="\t"}''!seen[$1]++' 10BestHits.txt > BLAST_1st_hit.txt
    #check which seqs got a hit
    gawk 'BEGIN{FS="\t"}{print $1}' < BLAST_1st_hit.txt | \
        uniq > gothits.names
    #add no_hits flag
    seqkit seq -n $fasta > $fasta.names
    grep -v -w -F -f gothits.names $fasta.names | \
        sed -e 's/$/\tNo_significant_similarity_found/' >> BLAST_1st_hit.txt
    #add header
    sed -i '1 i\
    qseqid\t1st_hit\tqlen\tslen\tqstart\tqend\tsstart\tsend\tevalue\tlength\tnident\tmismatch\tgapopen\tgaps\tsstrand\tqcovs\tpident' \
    BLAST_1st_hit.txt

    #remove unnecessary files 
    rm *.names


Taxonomy assignment with the RDP-classifier

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
    ASV_fasta_tmp="ASVs_TagJumpFiltered_minmax.fasta"

    # select by size to only retain ASVs that are within the range of expected variation. In this case we set it to 400 up to 430 bps
    vsearch --fastx_filter $ASV_fasta \
            --fastq_minlen 400 \
            --fastq_maxlen 430 \
            --fastaout $ASV_fasta_tmp

    mv $ASV_fasta_tmp $ASV_fasta

    # Run RDP-classifier
    time rdp_classifier \
            -Xmx12g \
            classify \
            -t $reference_database \
            -f allrank \
            -o RDP.taxonomy.txt \
            -q $ASV_fasta

____________________________________________________


Get target taxa
~~~~~~~~~~~~~~~

| This part filters the ASV dataset to include only target taxonomic group for the following analyses. 
| For example, if you are interested in Hymenoptera, then discard all ASVs that do not match to the target taxon based on the user defined threshold (default = 0.8). 


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

|logo_BGE_small| |eufund| |chfund| |ukrifund|
