#!/usr/bin/bash

set -e -o pipefail

file1=$@

if [[ "$#" -ne 1 || ! -f "$file1" ]]; then
    echo "Need an existing input file"
    exit 1;
fi    

sed  -e 's/Sebastes_dalli/Sebastes_dalli #1/'  -e 's/Sebastes_hopkinsi/Sebastes_hopkinsi #1/'  -e 's/Sebastes_inermis/Sebastes_inermis #1/'  -e 's/Sebastes_minor/Sebastes_minor #1/'  -e 's/Sebastes_semicinctus/Sebastes_semicinctus #1/'  -e 's/Sebastes_steindachneri/Sebastes_steindachneri #1/'  $file1

