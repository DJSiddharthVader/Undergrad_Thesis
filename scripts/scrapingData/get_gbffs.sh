#!/bin/bash

#file to download gbffs from ncbi summary txt file

function downloadSuffix(){
    acc_source="$1"
    suffix="$2"
    urlfile="$3"
    mkdir -p "$(basename $urlfile | cut -d'.' -f1)"
#echo "$suffix $acc_source $urlfile $(basename $urlfile | cut -d'.' -f1)"
    cd "$(basename $urlfile | cut -d'.' -f1)"
    tail -n +3 "$acc_source" | perl -pe "s/^.*(ftp:\/\/[^ \t]*)\t.*$/\1/g" | perl -pe "s/(^.*)(\/[^\/\t\n]*)/\1\2\2/g" | sed "s/$/$suffix/g" >| "$urlfile"
    numurls="$(wc -l $urlfile | cut -d' ' -f1)"
    parallel --eta --bar --eta -a "$urlfile" wget -q
    cd -
}

summary="/home/sid/thesis_SidReed/data/bacterialGBFFs/assembly_summary.txt"
#filesuffix="_genomic.gbff.gz" #file suffix for every organism
nuc_filesuffix="_cds_from_genomic.fna.gz" #file suffix for genomic dna fasta
prot_filesuffix="_translated_cds.faa.gz" #file suffix for protien fasta, all headers match genomic fasta
outfile="$1"
downloadSuffix "$summary" "$nuc_filesuffix"  "nuc_$outfile"
downloadSuffix "$summary" "$prot_filesuffix"  "prot_$outfile"

#DEPRECIATED
#GETTING ASSEMBLYSUMMARY
#ftpurl="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" #tsv-ish of bacteria, accs, names and ftp urls for everything in refseq
#filesuffix="_genomic.gbff.gz" #file suffix for every organism
#outfile="all_refseq_bacteria_ftp_paths_110918.txt" #outputfile with all ftp links

#WORKED EARILIER
#filesuffix="_genomic.gbk.gz" #file suffix for every organism
#cat "$summary" | perl -pe 's/(\/[^\/]+?)$/\1\1/g' | sed "s/$/$prot_filesuffix/g" >| "prot_$outfile"
#cat "$summary" | perl -pe 's/(\/[^\/]+?)$/\1\1/g' | sed "s/$/$nuc_filesuffix/g" >| "nuc_$outfile"
#echo "Got all paths in nuc_$outfile , prot_$outfile"
#for url in $(cat "$outfile");
#do
#    wget $url
#done

#COMMENTS ON SUBSTITUTIONS
#wget -qO- $ftpurl | tee "$rawsummary" | sed 's/\t\{1,\}/,/g' | awk -F "," '{print $0 ~ "latest" && $0 ~ "Complete Genome"?$(NF-1):""}' | sed '/^$/d' | perl -pe 's/(\/[^\/]+?)$/\1\1/g' | sed "s/$/$filesuffix/g" >| "$outfile"
#wget -qO- $ftpurl #get file from ncbi, quiet, direct to stdout
#tee "$rawsummary" #write file to file, continues to pass through stdout
#sed 's/\t\{1,\}/,/g' #sub tabs, multitabs for commas
#awk -F "," '{print $0 ~ "latest" && $0 ~ "Complete Genome"?$(NF-1):""}' #for rows where the genome is complete and the latest, print the 2nd last column (ftp link for organism) else print nothing
#sed '/^$/d' # remove blank lines
#perl -pe 's/(\/[^\/]+?)$/\1\1/g' #fix links so they point to specific file name without suffix
#sed "s/$/$filesuffix/g" #add suffix to the end of each ftp link (line)
## >| "$outfile" # write to output file specified above
