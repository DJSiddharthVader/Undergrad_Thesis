#!/bin/bash

dedupFasta() {
    sed -e '/^>/s/$/@/' -e 's/^>/#/' "$1" |\
        tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
        sort -u -t ' ' -f -k1,1 |\
        sed -e 's/^/>/' -e 's/\t/\n/' |\
        tail -n+2 > tmp_dedup
    mv tmp_dedup "$1"
}

rmEmptyHeads (){
    awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' "$1" > tmp_rmempt
    mv tmp_rmempt "$1"
}

extractFastas() {
    mkdir -p "$2"
    mkdir -p "$3"
    while read file; do
        ~/thesis_SidReed/scripts/modified_get_features.pl -a -t "$file" >| "$2""/"`basename $file`".faa" 2>> "$2""/errors.txt"
        ~/thesis_SidReed/scripts/modified_get_features.pl -c -t "$file" >| "$3""/"`basename $file`".fna" 2>> "$3""/errors.txt"
    done < "$1"
}

main() {
    mkdir -p "$2"
    pdir="protein_fastas"
    ndir="nucleotide_fastas"
    extractFastas "$1" "$2/$pdir" "$2/$ndir"
    cd "$2"
    cat "$pdir"/*.faa >| all_proteins.faa
    #dedupFasta all_proteins.faa
    rmEmptyHeads all_proteins.faa
    cat "$ndir"/*.fna >| all_nucleotides.fna
    #dedupFasta all_nucleotides.fna
    rmEmptyHeads all_nucleotides.fna
    cd -
}

extractMixedFastas() {
    mkdir -p "$2" #mixed_fastas
    while read file; do
        ~/thesis_SidReed/scripts/modified_get_features.pl -b -t "$file" >| "$2/`basename $file | cut -d'.' -f-2`.fnaa" 2>> "$2/errors.txt"
    done < "$1"
}

seperateMixedFastas() {
    mkdir -p "$2" #protein_fastas
    mkdir -p "$3" #nucleotide_fastas
    for file in "$1"/*.fnaa; do
        python ~/thesis_SidReed/scripts/seperateFasta.py "$file" >| tmp_fnames
        base=`basename $file | cut -d'.' -f-2`
        protname=`tail -1 tmp_fnames`
        mv "$protname" "./$2/$base".faa
        nucname=`head -1 tmp_fnames`
        mv "$nucname" "./$3/$base".fna
    done
}

main2() {
    mkdir -p "$2"
    cd "$2"
    mdir="mixed_fastas"
    extractMixedFastas ../"$1" "$mdir"
    pdir="protein_fastas"
    ndir="nucleotide_fastas"
    seperateMixedFastas "$mdir" "$pdir" "$ndir"
    cat "$pdir"/*.faa >| all_proteins.faa
    dedupFasta all_proteins.faa
    rmEmptyHeads all_proteins.faa
    cat "$ndir"/*.fna >| all_nucleotides.fna
    dedupFasta all_nucleotides.fna
    rmEmptyHeads all_nucleotides.fna

}

#need conda environment activated, cant do it in script for some reason so just make sure its activated first
shopt -s extglob #allows for globbing in scripts
main2 "$1" "$2" #$1 is file with paths to gbff files, $2 is genus name

#DEPRECIATED
