#!/bin/bash

set -euo pipefail

# Input files
vcf_file="$1"
phenotype="$2"
form_gwas="formatted_gwas.tsv"
output_vcf=$3

# Step : reformated vcf file (delete all _ in the vcf file)
sed -i 's/_//g' $vcf_file

# Step : Genotyping plink
plink --vcf $vcf_file --make-bed --allow-extra-chr --chr-set 9 --out genotype

# Step : Make plink gwas output # change to --assoc or --linear depending to the phenotype
plink --bfile genotype --pheno $phenotype --pheno-name PHENO --linear --allow-no-sex \
    --allow-extra-chr --out plink

# Step : Formatting plink gwas output
awk '{gsub(/^ +/, ""); gsub(/ +/, "\t")}1' plink.assoc > "$form_gwas"
echo "Formatted GWAS output saved to $form_gwas"

# Step : Add p_value info from gwas output to the vcf ref
touch "$output_vcf"
python3 add_p_value.py $form_gwas $vcf_file $output_vcf

# # Step : Normalization
# bcftools norm $output_vcf -m -any > $output_vcf.norm.vcf
# echo "Normalized VCF saved to $output_vcf.norm.vcf"

# bash plink_to_vcf.bash simulation/binary_data/binary.decomposed.normalized.svs.vcf simulation/binary_data/binary.plink.phenotype.tsv binary.plink_vcf.vcf
# bash plink_to_vcf.bash simulation/quantitative_data/quantitative.decomposed.normalized.svs.vcf simulation/quantitative_data/quantitative.plink.phenotype.tsv quantitative.plink_vcf.vcf

