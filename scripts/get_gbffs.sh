#!/bin/bash

#file to download gbffs from ncbi summary txt file

ftpurl="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" #tsv-ish of bacteria, accs, names and ftp urls for everything in refseq
rawsummary="assembly_summary_110918SR.txt"
filesuffix="_genomic.gbff.gz" #file suffix for every organism
outfile="all_refseq_bacteria_ftp_paths_110918.txt" #outputfile with all ftp links

wget -qO- $ftpurl | tee "$rawsummary" | sed 's/\t\{1,\}/,/g' | awk -F "," '{print $0 ~ "latest" && $0 ~ "Complete Genome"?$(NF-1):""}' | sed '/^$/d' | perl -pe 's/(\/[^\/]+?)$/\1\1/g' | sed "s/$/$filesuffix/g" >| "$outfile"
echo "Got all paths in $outfile"

for url in $(cat "$outfile");
do
    wget $url
done


#wget -qO- $ftpurl #get file from ncbi, quiet, direct to stdout
#tee "$rawsummary" #write file to file, continues to pass through stdout
#sed 's/\t\{1,\}/,/g' #sub tabs, multitabs for commas
#awk -F "," '{print $0 ~ "latest" && $0 ~ "Complete Genome"?$(NF-1):""}' #for rows where the genome is complete and the latest, print the 2nd last column (ftp link for organism) else print nothing
#sed '/^$/d' # remove blank lines
#perl -pe 's/(\/[^\/]+?)$/\1\1/g' #fix links so they point to specific file name without suffix
#sed "s/$/$filesuffix/g" #add suffix to the end of each ftp link (line)
## >| "$outfile" # write to output file specified above
