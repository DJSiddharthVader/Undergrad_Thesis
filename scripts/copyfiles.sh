#!/bin/bash

counter=1

for file in "$1"*.newick;
do
    cp $file "$2"
    counter=$((counter+1))
    if (( $counter > $3 ))
    then
        break
    fi
done
cp "$1"species.newick "$2"
