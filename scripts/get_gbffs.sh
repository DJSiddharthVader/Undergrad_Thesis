#!/bin/bash

#file to download gbffs from ncbi summary txt file

# 1 assembly_accession
# 2 bioproject
# 3 biosample
# 4 wgs_master
# 5 refseq_category
# 6 taxid
# 7 species_taxid
# 8 organism_name
# 9 infraspecific_name
# 10 isolate
# 11 version_status
# 12 assembly_level
# 13 release_type
# 14 genome_rep
# 15 seq_rel_date
# 16 asm_name
# 17 submitter
# 18 gbrs_paired_asm
# 19 paired_asm_comp
# 20 ftp_path
# 21 excluded_from_refseq
# 22 relation_to_type_material
ftpurl="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt"
rawsummary="assembly_summary_110918SR.txt"
filesuffix="_genomic.gbff.gz"
outfile="all_refseq_bacteria_ftp_paths_110918.txt"

wget -qO- $ftpurl | tee "$rawsummary" | sed 's/\t\{1,\}/,/g' | awk -F "," '{print $0 ~ "latest" && $0 ~ "Complete Genome"?$(NF-1):""}' | sed '/^$/d' | perl -pe 's/(\/[^\/]+?)$/\1\1/g' | sed "s/$/$filesuffix/g" >| "$outfile"
echo "Got all paths in $outfile"

for url in $(cat "$outfile");
do
    wget $url
done
