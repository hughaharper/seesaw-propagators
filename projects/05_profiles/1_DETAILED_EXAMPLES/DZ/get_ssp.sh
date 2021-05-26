#!/usr/bin/env bash

string=$1
output_name=$2

grep -A50 ${string}  output.xy |\
    awk '{ if ($1 != ">") print $0; else print ">" }' |\
        awk '{ if (NF <=2) print $0 }' | sed '1,3d' |\
            sed '/>/q' | sed '$d' | awk '{ print $1, $2}' > ${output_name}
