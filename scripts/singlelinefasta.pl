#!/usr/bin/bash

file1=$@

perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' $file1


