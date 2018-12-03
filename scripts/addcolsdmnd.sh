#!/bin/bash

echo "qseqid qlen length sseqid slen qcovhsp pident score bitscore evalue" | tr [:blank:] \\t >| tmp_tabs
cat "$1" >> tmp_tabs
mv tmp_tabs "$1"

