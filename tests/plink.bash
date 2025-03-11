#!/bin/bash

set -euo pipefail

# Input arguments
vcf_file=$1
pg_file=$2
phenotype=$3
snarl_paths=$4
snarl_freq=$5
snp_only_vcf=$6 # "True" or "False"

# Output file names
output_vcf_plink="plink.vcf"
output_vcf_truth="truth.vcf"

# Generate Plink VCF from GWAS output
bash plink_to_vcf.bash "$vcf_file" "$phenotype" "$output_vcf_plink"

# Generate truth VCF from pangenome graph
bash truth_to_vcf.bash "$pg_file" "$snarl_paths" "$snarl_freq" "$output_vcf_truth"

# Check if SNP filtering is required
if [[ "$snp_only_vcf" == "True" ]]; then
    echo "Filtering SNPs from VCF files..."
    bash only_snp_vcf.bash "$output_vcf_truth" "${output_vcf_truth}.snp.vcf"
    bash only_snp_vcf.bash "$output_vcf_plink" "${output_vcf_plink}.snp.vcf"

    output_vcf_truth="${output_vcf_truth}.snp.vcf"
    output_vcf_plink="${output_vcf_plink}.snp.vcf"
fi

# Compare the processed VCF files
echo "Comparing VCF files..."
python3 compare_truth_vcf.py "$output_vcf_truth" "$output_vcf_plink"

echo "Pipeline completed successfully."

# BINARY
# bash tests/plink.bash tests/simulation/binary_data/binary.decomposed.vcf ../data/binary_data/pg.full.pg \
# ../data/binary_data/binary.plink.phenotype.tsv ../data/binary_data/list_snarl_paths.tsv \
# ../data/binary_data/pg.snarls.freq.tsv True






# -------------------------------- GITLAB --------------------------------
# BINARY
# bash plink.bash ../data/binary_data/binary.decomposed.vcf ../data/binary_data/pg.full.pg \
# ../data/binary_data/binary.plink.phenotype.tsv ../data/binary_data/list_snarl_paths.tsv \
# ../data/binary_data/pg.snarls.freq.tsv ../data/binary_data/plink_vcf.vcf \
# ../data/binary_data/truth.deconstruct.norm.vcf

# QUANTITATIVE
# bash plink.bash ../data/binary_data/binary.decomposed.vcf ../data/binary_data/pg.full.pg \
# ../data/binary_data/binary.plink.phenotype.tsv ../data/binary_data/list_snarl_paths.tsv \
# ../data/binary_data/pg.snarls.freq.tsv ../data/binary_data/plink_vcf.vcf \
# ../data/binary_data/truth.deconstruct.norm.vcf
