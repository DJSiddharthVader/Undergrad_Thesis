#!/bin/bash
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} $2 {print ">"$0}' "$1" > tmp_rmempt
mv tmp_rmempt "$1"
