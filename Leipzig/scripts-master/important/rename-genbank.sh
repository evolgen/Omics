#!/usr/bin/bash

for f in *.gb; do
mv "$f" "$(basename "$f" .gb).gbk"
done


