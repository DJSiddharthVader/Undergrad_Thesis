#!/bin/bash
#$1 is the all proteins fasta file
#$2 is the genusname
fasta="$1"
db="mmseq_db_$2"
out="clusters_$2"
mmseqs createdb "$fasta" "$db"
mmseqs cluster -e 0.05 -c 0.8 --min-seq-id 0.85 --alignment-mode 3 --cluster-mode 1 "$db" "$out" mmseq_tmp
mmseqs createtsv "$db" "$db" "$out" "$out".tsv
\rm -rf mmseqs_tmp