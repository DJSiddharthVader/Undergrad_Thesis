#!/bin/bash

outdir="$2"
pdir="prot_$2"
mkdir -p "$pdir"
ndir="nuc_$2"
mkdir -p "$ndir"

while read file; do
    ~/thesis_SidReed/scripts/modified_get_features.pl -a -t "$file" >| "$pdir""/"`basename $file`".faa" 2>> "$pdir""/errors.txt"
    ~/thesis_SidReed/scripts/modified_get_features.pl -c -t "$file" >| "$ndir""/"`basename $file`".fna" 2>> "$ndir""/errors.txt"
done < $1

#if [ "$3" = 'prot' ]; then
#    while read file; do
#         ~/thesis_SidReed/scripts/modified_get_features.pl -a -t "$file" >| "$outdir""/"`basename $file`".faa" 2>> "$outdir""/errors.txt"
#    done < $1
#elif [ "$3" = 'nuc' ]; then
#    while read file; do
#         ~/thesis_SidReed/scripts/modified_get_features.pl -c -t "$file" >| "$outdir""/"`basename $file`".fna" 2>> "$outdir""/errors.txt"
#    done < $1
#else
#    echo "Error: arg 3 must be prot or nuc"
#fi
