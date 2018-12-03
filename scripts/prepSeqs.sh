#!/bin/bash

extractFastas() {
    pdir="prot_$2"
    mkdir -p "$pdir"
    ndir="nuc_$2"
    mkdir -p "$ndir"
    while read file; do
        ~/thesis_SidReed/scripts/modified_get_features.pl -a -t "$file" >| "$pdir""/"`basename $file`".faa" 2>> "$pdir""/errors.txt"
        ~/thesis_SidReed/scripts/modified_get_features.pl -c -t "$file" >| "$ndir""/"`basename $file`".fna" 2>> "$ndir""/errors.txt"
    done < "$1"
}

extractMixedFastas() {
    mdir="mix_$2"
    mkdir -p "$mdir"
    while read file; do
        ~/thesis_SidReed/scripts/modified_get_features.pl -b -t "$file" >| "$mdir/`basename $file | cut -d'.' -f-2`.fnaa" 2>> "$mdir/errors.txt"
    done < "$1"
}

seperateMixedFastas() {
    pdir="prot_$1"
    mkdir -p "$pdir"
    ndir="nuc_$1"
    mkdir -p "$ndir"
    for file in "mix_$1"/*.fnaa; do
        python ~/thesis_SidReed/scripts/seperateFasta.py "$file" >| tmp_fnames
        base=`basename $file | cut -d'.' -f-2`
        protname=`tail -1 tmp_fnames`
        mv "$protname" "./$pdir/$base".faa
        nucname=`head -1 tmp_fnames`
        mv "$nucname" "./$ndir/$base".fna
    done
}

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

main() {
    extractFastas "$1" "$2"
    cat "prot_$2"/*.faa >| all_proteins.faa
    dedupFasta all_proteins.faa
    rmEmptyHeads all_proteins.faa
    cat "nuc_$2"/*.fna >| all_nucleotides.fna
    dedupFasta all_nucleotides.fna
    #rmEmptyHeads all_nucleotides.fna
}

main2() {
    extractMixedFastas "$1" "$2"
    seperateMixedFastas "$2"
    cat "prot_$2"/*.faa >| all_proteins.faa
    dedupFasta all_proteins.faa
    rmEmptyHeads all_proteins.faa
    cat "nuc_$2"/*.fna >| all_nucleotides.fna
    dedupFasta all_nucleotides.fna
    rmEmptyHeads all_nucleotides.fna

}

shopt -s extglob #allows for globbing in scripts
#main "$1" "$2"
#need conda environment activated, cant do it in script for some reason so just make sure its activated first
main2 "$1" "$2"


