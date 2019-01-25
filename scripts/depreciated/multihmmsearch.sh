#!/bin/bash

# $1 is dir of hmms to run against
# $2 is dir of query fasta files

outdir=`basename "$1"`"_"`basename "$2"`"_hmmresults"
mkdir -p "$outdir"
echo "$outdir"
for fasta in $2/*.faa
do
    fd="$outdir""/"`basename "$fasta"`
    echo "$fd"
    mkdir -p "$fd"
    for hmm in $1/*.hmm;
    do
        /home/sid/thesis_SidReed/apps/bin/hmmsearch "$hmm" "$fasta" > "$fd"/`basename "$hmm"`"_result.txt"
    done
done
