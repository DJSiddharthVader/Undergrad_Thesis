#!/bin/bash
fam="$1"
fdir="gene_tree_files/trees/$fam"
mkdir -p "$fdir"
mv "$fam"_* \
   "$fam".mcmc \
   "$fam".parts \
   "$fam".con.tre \
   "$fam".run1.p \
   "$fam".run1.t \
   "$fam".run2.p \
   "$fam".run2.t \
   "$fam".trprobs \
   "$fam".tstat \
   "$fam".vstat \
   "$fdir" 2> /dev/null
