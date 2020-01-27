#!/bin/bash

#use ./get_features_csv.sh dir_of_gbffs/ >| outfile.csv
#takes ~50minutes on the 10,614 gbffs I have, probably speed it up by grepping single column for every file at once, paste cols?

declare -a gcounts=("Genes (total)" "CDS (total)" "Genes (coding)" "CDS (coding)" "Genes (RNA)" "tRNAs" "ncRNAs" "Pseudo Genes (total)")
headers="Filepath~GCF~Accession Number~Organism~Genus~Taxonomy~IsComplete~NCBI_CRISPR_Arrays~Genes (total)~CDS (total)~Genes (coding)~CDS (coding)~Genes (RNA)~tRNAs~ncRNAs~Pseudo Genes (total)"
singlefileinfo() {  #get info for a single gbfffile, , delimited
    #identifying info
    local GCF=$(grep -m 1 'Assembly' "$1" | cut -d':' -f2 | tr -d ' ') #gca for linking back to crispr annoattions
    local accnum=$(head -1 "$1" | sed 's/^\w\+\s\+\(\w\+\)\s.*$/\1/g') #accession
    local orgname=$(grep -m 1 ORGANISM "$1" | sed 's/^\s\+ORGANISM\s\+//') #name
    local genus=$(echo $orgname | perl -pe 's/^(\w.*?)\s.*$/\1/') #genus only
    local taxonomy=$(sed -n '/ORGANISM/,/REFERENCE/{/(ORGANISM|REFERENCE|\n)/d;p}' "$1" | head -n -1 | sed 's/ORGANISM//' | tail -n+2 | tr -d '\n' | tr -d ' ') #gets taxonomy until genus
    local complet=$(grep -m1 COMPLETENESS "$1" | cut -d':' -f2)
    local crisprs=$(grep -m1 CRISPR.*:: "$1" | cut -d':' -f3 | tr -d ' ')
    if [[ $crisprs == '' ]];then
        crisprs=0
    fi
    info=("$GCF" "$accnum" "$orgname" "$genus" "$taxonomy" "$complet" "$crisprs" )
    #gene counts
    for gcount in "${gcounts[@]}" #get gene type counts for bacteria
    do
        pattern="^\s*""$gcount"
        local tmp=$(grep -m 1 "$pattern" "$1" | tr -d ',' | perl -pe 's/^.*:: (\d+)[^\d]*$/\1/')
        info+=("$tmp")
    done
    #join with ~ and print
    outvar="$1" #absolute file path
    for i in "${info[@]}" #for each piece of info
    do
        outvar="$outvar~$i" #join with ~
    done
    echo "$outvar"
}
function main() {
    gbff_dir="$1"
    echo "$headers"
    for file in $gbff_dir/* #dir of gbff files
    do
        singlefileinfo $(readlink -e "$file")
    done
}

main "$1"

