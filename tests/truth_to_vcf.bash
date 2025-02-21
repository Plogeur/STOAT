#!/bin/bash

set -euo pipefail

graph=$1 # .full.pg
path=$2 # .path
freq=$3 # .freq.tsv (verity simulation)
output_vcf=$4 # .vcf

# Create a fasta ref
echo "ref" > paths.tsv
vg paths --extract-fasta -p paths.tsv --xg $graph > ref.fa
samtools faidx ref.fa

# Deconstruct / identify position/chr from the pangenome graph
vg deconstruct -p ref -a $graph > merged.deconstruct.vcf

# Left align and normalize using ref.fasta
bcftools norm merged.deconstruct.vcf -f ref.fa > merged.deconstruct.norm.vcf

# Add frequency verity info to the vcf
touch "$output_vcf"
python3 tests/add_vcf_freq.py merged.deconstruct.norm.vcf $path $freq $output_vcf

# bash tests/truth_to_vcf.bash tests/simulation/binary_data/pg.full.pg tests/simulation/binary_data/list_snarl_paths.tsv tests/simulation/binary_data/pg.snarls.freq.tsv tests/simulation/binary_data/truth.deconstruct.norm.vcf
