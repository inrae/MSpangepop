#!/bin/bash

in="$1"
out_dir="$2"
out=""

while read -r line; do
    if [[ "$line" == \>* ]]; then
        [ -n "$out" ] && exec 3>&-
        out="${line##>}.fasta"
        exec 3> "$out_dir/$out"
        echo "$line" >&3
    else
        echo "$line" >&3
    fi
done < "$in"

[ -n "$out" ] && exec 3>&-
