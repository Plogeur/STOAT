#!/bin/bash

# Input files
gwas_file="$1"
vcf_file_ref="$2"
vcf_file_at="$3"
out_file="formatted_gwas.tsv"

# Step 1: Format the GWAS file by removing the leading tab in the header
awk '{gsub(/^ +/, ""); gsub(/ +/, "\t")}1' "$gwas_file" > "$out_file"
echo "Formatted GWAS output saved to $out_file."

# Step 2: Call the Python script add_num_sample.py
python3 tests/add_info_plink.py "$vcf_file_ref" "$out_file"
echo "Python script add_info_plink.py executed."

mv updated_formatted_gwas.tsv tests/plink_tests_output/$out_file
rm formatted_gwas.tsv

# bash tests/modify_plink.bash tests/plink_tests_output/binary_plink.assoc tests/plink_tests_output/binary.merged.decomposed.normalized.svs.vcf tests/plink_tests_output/binary.merged.vcf
