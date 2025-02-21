#!/bin/bash

set -euo pipefail

# Input files
gwas_file="$1"
vcf_file_ref="$2"
form_gwas="tests/formatted_gwas.tsv"
output_vcf="tests/updated_vcf.vcf"

# Step 1: Formatting plink gwas output
awk '{gsub(/^ +/, ""); gsub(/ +/, "\t")}1' "$gwas_file" > "$form_gwas"
echo "Formatted GWAS output saved to $form_gwas"

# Step 2: Add p_value info from gwas output to the vcf ref
touch "$output_vcf"
python3 tests/add_p_value.py $form_gwas $vcf_file_ref $output_vcf

# Step 3: Normalization
bcftools norm $output_vcf -m -any > $output_vcf.norm.vcf
echo "Normalized VCF saved to $output_vcf.norm.vcf"

# bash tests/plink_to_vcf.bash tests/plink_tests_output/binary_plink.natif.assoc tests/plink_tests_output/binary.merged.decomposed.normalized.svs.vcf
