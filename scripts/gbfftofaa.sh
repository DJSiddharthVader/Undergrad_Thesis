#!/bin/bash

outdir=$(basename "$1" | cut -d'.' -f1)
mkdir -p $outdir
while read file;
do
     ~/thesis_SidReed/scripts/modified_get_features.pl -a -t "$file" >| "$outdir""/"`basename $file`".faa" 2>> "$outdir""/errors.txt"
done < $1
