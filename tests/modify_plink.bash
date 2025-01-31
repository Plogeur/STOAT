#!/bin/bash

# Input files
gwas_file="$1"
vcf_file="$2"
out_file="formatted_gwas.tsv"

# Check if files exist
if [[ ! -f "$gwas_file" || ! -f "$vcf_file" ]]; then
    echo "Error: One or both input files do not exist."
    exit 1
fi

# Step 1: Format the GWAS file by removing the leading tab in the header
awk 'NR==1{sub(/^\t/, ""); print} NR>1{print}' "$gwas_file" > "$out_file"
echo "Formatted GWAS output saved to $out_file."

# Step 2: Replace the second column in GWAS with the SS ID from VCF
# Save the starting index of non-header rows in the VCF
vcf_data_start=$(awk 'BEGIN {start_line = 0} /^#/ {next} {start_line = NR; exit} END {print start_line}' "$vcf_file")

awk -v start_index="$vcf_data_start" 'BEGIN { FS = OFS = "\t" } \
    NR==FNR { \
        if ($0 ~ /^#/) { next } \
        if (NF >= 8) { \
            split($8, info, ";"); \
            for (i in info) if (info[i] ~ /^SS=/) { \
                split(info[i], ss, "="); ss_map[NR-start_index+1] = ss[2]; break \
            } \
        } \
        next \
    } \
    FNR==1 { print; next } \
    { \
        if ((FNR - 1) in ss_map) $2 = ss_map[FNR - 1]; \
        print \
    }' "$vcf_file" "$out_file" > tmp_gwas.txt

mv tmp_gwas.txt "$out_file"
echo "Second column in GWAS replaced with SS IDs."

# Step 3: Call the Python script add_num_sample.py
python3 tests/add_info_plink.py "$vcf_file" "$out_file"
echo "Python script add_info_plink.py executed."

mv updated_formatted_gwas.tsv tests/plink_tests_output/$out_file
rm formatted_gwas.tsv
