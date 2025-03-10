#!/bin/bash

INPUT_VCF=$1
OUTPUT_VCF=$2

# Check if the input file is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input.vcf output.vcf"
    exit 1
fi

# Filter SNPs from the VCF file
bcftools view -v snps "$INPUT_VCF" -o "$OUTPUT_VCF" -O v

echo "Filtered SNPs saved to $OUTPUT_VCF"

# bash tests/only_snp_vcf.bash tests/truth.deconstruct.norm.vcf tests/truth.deconstruct.snp.norm.vcf 