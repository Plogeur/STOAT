#!/bin/bash

set -euo pipefail

graph=$1 # .full.pg
path=$2 # .path
freq=$3 # .freq.tsv (verity simulation)
output_vcf=$4 # .vcf

# Create a fasta ref
echo "ref" > tests/paths.tsv
vg paths --extract-fasta -p tests/paths.tsv --xg $graph > tests/ref.fa
samtools faidx tests/ref.fa

# Add fake haplotype to the pangenome graph
vg gbwt -x $graph -o tests/pg.cover.gbwt -P -n32
vg deconstruct -a -p ref -g tests/pg.cover.gbwt $graph > tests/deconstruct.vcf

# Left align and normalize using ref.fasta
bcftools norm tests/deconstruct.vcf -f tests/ref.fa > tests/deconstruct.norm.vcf

# Add frequency verity info to the vcf
touch "$output_vcf"
python3 tests/add_vcf_freq.py tests/deconstruct.norm.vcf $path $freq $output_vcf

# bash tests/truth_to_vcf.bash tests/simulation/binary_data/pg.full.pg tests/simulation/binary_data/snarl_paths.tsv tests/simulation/binary_data/pg.snarls.freq.tsv tests/truth.deconstruct.norm.vcf

# bash tests/truth_to_vcf.bash tests/simulation/quantitative_data/pg.full.pg tests/simulation/quantitative_data/snarl_paths.tsv tests/simulation/quantitative_data/pg.snarls.freq.tsv tests/quantitative.truth.deconstruct.norm.vcf
