#!/bin/bash

set -euo pipefail

# Input files
vcf_file="$1"
phenotype="$2"
form_gwas="tests/formatted_gwas.tsv"
output_vcf="tests/plink_vcf.vcf"

# Step : reformated vcf file (delete all _ in the vcf file)
sed -i 's/_//g' $vcf_file

# Step : Genotyping plink
plink --vcf $vcf_file --make-bed --allow-extra-chr --chr-set 9 --out tests/genotype

# Step : Make plink gwas output
plink --bfile tests/genotype --pheno $phenotype --pheno-name PHENO --assoc --allow-no-sex \
    --allow-extra-chr --out tests/plink

# Step : Formatting plink gwas output
awk '{gsub(/^ +/, ""); gsub(/ +/, "\t")}1' tests/plink.assoc > "$form_gwas"
echo "Formatted GWAS output saved to $form_gwas"

# Step : Add p_value info from gwas output to the vcf ref
touch "$output_vcf"
python3 tests/add_p_value.py $form_gwas $vcf_file $output_vcf

# Step : Normalization
bcftools norm $output_vcf -m -any > $output_vcf.norm.vcf
echo "Normalized VCF saved to $output_vcf.norm.vcf"

# bash tests/plink_to_vcf.bash tests/simulation/binary_data/binary.decomposed.vcf tests/simulation/binary_data/binary.plink.phenotype.tsv

