#!/bin/bash
#$1 is genus name, capitalized
function dedupFasta() {
    sed -e '/^>/s/$/@/' -e 's/^>/#/' "$1" |\
        tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
        sort -u -t ' ' -f -k1,1 |\
        sed -e 's/^/>/' -e 's/\t/\n/' |\
        tail -n+2 > tmp_dedup
    mv tmp_dedup "$1"
}
function rmEmptyHeads(){
    awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' "$1" > tmp_rmempt
    mv tmp_rmempt "$1"
}
function main(){
    genus="$1"
}
main "$1"
