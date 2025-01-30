#!/bin/bash

dir_vcf=$1

# Decompress all .vcf files in the directory
for file in $dir_vcf/*.vcf.gz; do
    [ -e "$file" ] || continue # Skip if no files match
    gunzip "$file"
done

# Re-compress all .vcf files in the directory
for file in $dir_vcf/*.vcf; do
    [ -e "$file" ] || continue # Skip if no files match
    bgzip "$file"
    # bcftools norm "$file"
    bcftools index "$file"
done

bcftools merge $dir_vcf/*.vcf.gz --threads 6 -0 -Oz -o $dir_vcf/merged_output_0.vcf
