#!/bin/bash

sed -e '/^>/s/$/@/' -e 's/^>/#/' "$1" |\
    tr -d '\n' | tr "#" "\n" | tr "@" "\t" |\
    sort -u -t ' ' -f -k1,1 |\
    sed -e 's/^/>/' -e 's/\t/\n/' |\
    tail -n+2 > "$2"
