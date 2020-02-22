#!/bin/bash

#file to download gene sequences from urllust
function downloadList(){
    urllist="$1"
    output="$2"
    mode="$(basename $urllist | cut -d'_' -f1)"
    case $mode in
        'nucleotide')
            outdir="$output/nucleotide"
            ;;
        'protein')
            outdir="$output/protein"
            ;;
        'genome')
            outdir="$output/genome"
            ;;
        *)
            echo "Invalid mode, pick from {nuc|prot}"
            exit 1
            ;;
    esac
    mkdir -p $outdir && cd $outdir #directory
    parallel --eta --bar -a "$urllist" wget -q -a '../download_errors.txt'
    cd -
}
function main() {
    filelist_dir="$1"
    outdir="$2"
    for file in $filelist_dir/*FTPURLs.txt;
    do
        downloadList "$(readlink -e $file)" $outdir
    done
}
filelist_dir="$(readlink -e $1)"
output="$2"
main $filelist_dir $output

#DEPRECIATED
#GETTING ASSEMBLYSUMMARY
#ftpurl="ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt" #tsv-ish of bacteria, accs, names and ftp urls for everything in refseq
#filesuffix="_genomic.gbff.gz" #file suffix for every organism
#outfile="all_refseq_bacteria_ftp_paths_110918.txt" #outputfile with all ftp links

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
#function downloadSuffix(){
#    assembly_summary="$1"
#    output="$2"
#    mode="$3"
#    mkdir -p $output && cd $output #directory
#    case $mode in
#        'nuc')
#            urllist="./nucleotide_data_urls.txt"
#            downloaddir="./nucleotide"
#            suffix="_cds_from_genomic.fna.gz" #file suffix for genomic dna fasta
#            ;;
#        'prot')
#            urllist="./protein_data_urls.txt"
#            downloaddir="./protein"
#            suffix="_translated_cds.faa.gz" #file suffix for protien fasta, all headers match genomic fasta
#            ;;
#        *)
#            echo "Invalid mode, pick from {nuc|prot}"
#            exit 1
#            ;;
#    esac
#    #create file of fasta ftp urls for download, nucleotide or protein
#    tail -n +3 "$assembly_summary" |\
#               grep 'Complete Genome' |\
#               perl -pe "s/^.*(ftp:\/\/[^ \t]*)\t.*$/\1/g" |\
#               perl -pe "s/(^.*)(\/[^\/\t\n]*)/\1\2\2/g" |\
#               sed "s/$/$suffix/g" >| "$urllist"
#    mkdir -p $downloaddir && cd $downloaddir
#    parallel --eta --bar -a "../$urllist" wget -q
#}
#assembly_summary="$(readlink -e $1)"
#output="$2"
#downloadSuffix $assembly_summary $output 'nuc'
#downloadSuffix $assembly_summary $output 'prot'

