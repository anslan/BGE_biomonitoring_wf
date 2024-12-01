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

____________________________________________________

|logo_BGE_small| |eufund| |chfund| |ukrifund|
