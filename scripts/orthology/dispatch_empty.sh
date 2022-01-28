#!/bin/bash

#SBATCH --job-name=empty
#SBATCH -A co_genomicdata
#SBATCH --partition=savio,savio2,savio2_bigmem #,savio_bigmem
#SBATCH --time=100:00
#SBATCH -N 1
#SBATCH --qos=savio_lowprio
#SBATCH -o dispatch.log


