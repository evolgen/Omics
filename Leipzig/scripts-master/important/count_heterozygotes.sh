#!/usr/bin/env bash

# Dobzhansky Center for Genome Bioinformatics, 2014
# Gaik Tamazian
# gaik (dot) tamazian (at) gmail (dot) com

# -------------------------------------------------

set -u

function usage {
	echo count_heterozygosity.sh - count the number of homozygous and 
	echo heterozygous mutations in a VCF file with a single invididual.
	echo
	echo Usage:
	echo	./count_heterozygosity.sh variants.vcf
	echo
	echo The script counts the numbers of homozygous and heterozygous
	echo mutations in the specified VCF file. The file must contain
	echo information on mutations for a single individual, otherwise an
	echo error is returned.
	
	exit 1
}

function check_single_individual {
	# The function checks if the specified VCF file contains mutations
	# of a single individual. It returns 1 if the VCF file satisfies it,
	# otherwise 0 is returned.

	local VCF_FILE=$1

	# get the number of columns for each string of the specified VCF
	# file
	COL_NUM=`awk '!/^($|#)/{ print NF }' "$VCF_FILE" | sort -nu` 
	MIN_COL_NUM=`echo -e "$COL_NUM" | head -n 1`
	MAX_COL_NUM=`echo -e "$COL_NUM" | tail -n 1`
	if [ $MAX_COL_NUM -ne $MIN_COL_NUM ]; then
		# if the least and the greates values of column numbers are not
		# equal to each other, then there are rows which contain
		# different number of columns
		return 0
	elif [ $MAX_COL_NUM -ne 10 ]; then
		# if all rows have the same name number of columns which differs
		# from 10, then the specified VCF file describes mutations of
		# more than one individual
		return 0
	else
		# all rows have 10 columns, so the specified VCF file contains
		# mutations of a single individual
		return 1
	fi
}

# We use an AWK script to output the numbers of homozygous and
# heterozygous mutations in the specified VCF file.

read -d '' GET_ZYGOSITY << 'EOF'
BEGIN {
	homozygotes = 0
	heterozygotes = 0
	OFS = "\\t"
}
!/^($|#)/{
	# scan the format column (#9) to determine which part of the 
	# individual info column corresponds to the genotype field
	n = split($9, format_fields, ":")
	genotype_field = 0
	for (i = 1; i <= n; i++) {
		if (format_fields[i] == "GT") {
			genotype_field = i
			break
		}
	}
	# get fields for an individual
	split($10, individual, ":")
	# check the right separator sign: "/" or "|"
	if (individual[genotype_field] ~ /\\//) {
		allele_separator = "/"
	} else {
		allele_separator = "|"
	}
	split(individual[genotype_field], alleles, allele_separator)
	# compare alleles to each other
	if (alleles[1] == alleles[2]) {
		homozygotes++
	} else {
		heterozygotes++
	}
}
END {
	print "# heterozygotes:", heterozygotes
	print "# homozygotes:", homozygotes
	print "# total:", homozygotes + heterozygotes
}
EOF

# End of the AWK script.

if [ $# -ne 1 ]; then
	usage
	exit 1
fi

# check if the specified VCF file exists
if [ ! -f "$1" ]; then
	>&2 echo Error: the specified file "$1" does not exist.
	exit 1
fi

# check if the specified VCF file contains mutations of the single
# individual 
check_single_individual "$1"
if [ $? -ne 1 ]; then
	>&2 echo Error: the specified VCF file contain mutations of multiple
	>&2 echo individuals
	exit 1
fi

# run the AWK script defined above
echo `basename "$1"`
awk "$GET_ZYGOSITY" "$1"


