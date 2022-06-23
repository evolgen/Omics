#! /bin/bash -c
exec < $1
while read line
do
$line2=$line+3
paste $line $line2 > $line.out
done
