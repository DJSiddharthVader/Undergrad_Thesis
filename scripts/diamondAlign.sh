#!/bin/bash
#$1 is dir with all_nucleotide.fna and all_proteins.faa files
cd "$1"
diamond makedb --in all_proteins.faa -d all_proteins
diamond blastp --more-sensitive --matrix BLOSUM62 -e 5e-2 --query-cover 85 --no-self-hits -f 6 qseqid qlen length sseqid slen qcovhsp pident score bitscore evalue -q all_proteins.faa -d all_proteins --out allvsallproteins.dmnd.out
/home/sid/thesis_SidReed/scripts/addcolsdmnd.sh allvsallproteins.dmnd.out
cd -
