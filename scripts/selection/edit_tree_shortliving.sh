#!/usr/bin/bash

set -e -o pipefail

file1=$@

if [[ "$#" -ne 1 || ! -f "$file1" ]]; then
    echo "Need an existing input file"
    exit 1;
fi    

sed -e 's/Helicolenus_lengerichi/Helicolenus_lengerichi{Foreground}/' -e 's/Sebastes_atrovirens/Sebastes_atrovirens{Foreground}/'  -e 's/Sebastes_cheni/Sebastes_cheni{Foreground}/'  -e 's/Sebastes_dalli/Sebastes_dalli{Foreground}/'  -e 's/Sebastes_emphaeus/Sebastes_emphaeus{Foreground}/'  -e 's/Sebastes_exsul/Sebastes_exsul{Foreground}/'  -e 's/Sebastes_hopkinsi/Sebastes_hopkinsi{Foreground}/'  -e 's/Sebastes_inermis/Sebastes_inermis{Foreground}/'  -e 's/Sebastes_lentiginosus/Sebastes_lentiginosus{Foreground}/'  -e 's/Sebastes_minor/Sebastes_minor{Foreground}/'  -e 's/Sebastes_rastrelliger/Sebastes_rastrelliger{Foreground}/' -e 's/Sebastes_schlegelii/Sebastes_schlegelii{Foreground}/'  -e 's/Sebastes_semicinctus/Sebastes_semicinctus{Foreground}/'  -e 's/Sebastes_serriceps/Sebastes_serriceps{Foreground}/'  -e 's/Sebastes_steindachneri/Sebastes_steindachneri{Foreground}/'  -e 's/Sebastes_ventricosus/Sebastes_ventricosus{Foreground}/' $file1

