import sys
from cyvcf2 import VCF

vcf_file = sys.argv[1]
gwas_file = sys.argv[2]
output_file = "updated_" + gwas_file

# Open the VCF file using cyvcf2
vcf = VCF(vcf_file)

with open(gwas_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Read the header from the GWAS file and append new column names
    header = infile.readline().strip()
    outfile.write(header + "\tNUM_SAMPLE\n")  # Add new columns to the header

    # Iterate through the GWAS file and VCF variants
    for line, variant in zip(infile, vcf):
        parts = line.strip().split("\t")

        # Get genotypes from the VCF for the current variant
        sample_genotypes = variant.genotypes  # List of [allele1, allele2, phased]

        # Count the number of samples with non-zero alleles
        num_sample = sum(sum(gt[:2]) > 0 for gt in sample_genotypes if gt[0] != -1)

        # Write the updated line to the output file
        outfile.write(line.strip() + f"\t{num_sample}\n")

print(f"Updated GWAS file with NUM_SAMPLE and AT columns saved to {output_file}")
