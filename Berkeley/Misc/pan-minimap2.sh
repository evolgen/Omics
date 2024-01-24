#!/bin/bash

set -e

while [[ "$@" != "" ]];
do
    minimap2 -cx asm20 $@
    shift
done

