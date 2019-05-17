#!/usr/bin/sh

set -e

if [ $# -ne 2 ]; then
    echo "Your command line contains $# arguments; only one ID-list needed";
    exit
fi

idlist=$1
threads=$2

if [[ $idlist != *.list ]]; then
    echo "File should end with .list"; exit
fi

number='^[0-9]+$'
if ! [[ $threads =~ $number ]] ; then
   printf "  #Threads: Not a valid number since, you gave:\t $threads\n"; exit;
fi

if [[ $threads == 0 ]] ; then
   printf "  #Threads: Cannot be 0\n"; exit;
fi

cat $idlist | while read line1; do
    if [ ${line1:0:3} != "[A-Z][A-Z]A-Z]" ]; then
      printf "The ID should be of SRR, yours is:\t $line1\n"; exit
    fi  

    initial1=$(echo ${line1:0:6})
    initial2=$(echo ${line1: -1})
    url="ftp://ftp.sra.ebi.ac.uk/vol1/srr/${initial1}/00${initial2}/$line1"
    printf "\t\tDownloading from $url\n"
    wget $url
    
done




: '

	SRR2048513	->	ftp://ftp.sra.ebi.ac.uk/vol1/srr/SRR204/003/SRR2048513
'
