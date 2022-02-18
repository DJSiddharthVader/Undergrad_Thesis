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
