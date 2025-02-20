#!/bin/bash

set -euo pipefail

graph=$1 # .full.pg
# freq=$2 # .freq.tsv (verity simulation)

# Create a fasta ref
echo "ref" > paths.tsv
vg paths --extract-fasta -p paths.tsv --xg $graph > ref.fa
samtools faidx ref.fa

# Deconstruct / identify position/chr from the pangenome graph
vg deconstruct -p ref -a $graph > merged.deconstruct.vcf

# Left align and normalize using ref.fasta
bcftools norm merged.deconstruct.vcf -f ref.fa > merged.deconstruct.norm.vcf

# python3 add_freq.py merged.deconstruct.norm.vcf $freq

