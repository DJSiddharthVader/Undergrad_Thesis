#!/bin/bash

fasta="../$1" #the all proteins fasta file
genus="$2" #is the genusname
mmseq_dir="$3"

db="mmseq_db_$genus"
out="clusters_$genus"
#run mmseq clustering
mkdir -p "$mmseq_dir" && cd "$mmseq_dir"
mmseqs createdb "$fasta" "$db" -v 3 >| mmseq.log 2>| mmseq.errors
mmseqs cluster -e 0.05 -c 0.8 --min-seq-id 0.85 --alignment-mode 3 --cluster-mode 1 "$db" "$out" mmseq_tmp >> mmseq.log 2>> mmseq.errors
mmseqs createtsv "$db" "$db" "$out" "$out".tsv >| mmseq.log 2>| mmseq.errors
\rm -rf mmseqs_tmp
cd -


#DEPRECIATED
##$1 is the all proteins fasta file
##$2 is the genusname
#fasta="$1"
#db="mmseq_db_$2"
#out="clusters_$2"
#mmseqs createdb "$fasta" "$db"
#mmseqs cluster -e 0.05 -c 0.8 --min-seq-id 0.85 --alignment-mode 3 --cluster-mode 1 "$db" "$out" mmseq_tmp >| mmseq.log 2>| mmseq.errors
#mmseqs createtsv "$db" "$db" "$out" "$out".tsv
#\rm -rf mmseqs_tmp
