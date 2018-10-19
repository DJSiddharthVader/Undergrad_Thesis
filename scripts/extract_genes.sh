#!/bin/bash

#will get all the genes with descriptive headers from a gbff file, format them as fasta and write to stdout,

gbffdir=$1

faadir="bacterialGeneFastas"
mkdir -p $faadir
for gfile in $gbffdir;
do
    home/sid/thesis_SidReed/scripts/modified_get_features.pl <(zcat $gfile) > "$faadir"/`basename "$gfile"`.faa: 2>> "$faadir"/modified_get_features_errorfile.txt
done

allfasta="$faadir"/all_bacterial_genes.faa
cat "$faadir"/*.faa >  "$allfasta"

blastdbdir="bacterial_blastdbs"
mkdir -p "$blastdbdir"
blastdbfile="$blastdbdir"/all_bacterial_genes.db
makeblastdb -in "$allfasta" -out "$blastdbfile" -dbtype prot -parse_seqids

blastp -db "$blastdbfile"
       -query
       -out
       -outfmt '7 qseqid qlen qstart qend length sseqid slen qcovs score bitscore evalue'
       -soft_masking true
       -use_sw_tback
