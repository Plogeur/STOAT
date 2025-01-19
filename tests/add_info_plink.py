import sys
from cyvcf2 import VCF

# Check if input files are provided
if len(sys.argv) < 3:
    print("Usage: python3 add_columns.py <vcf_file> <gwas_file>")
    sys.exit(1)

vcf_file = sys.argv[1]
gwas_file = sys.argv[2]
output_file = "updated_" + gwas_file

# Open the VCF file using cyvcf2
vcf = VCF(vcf_file)

try:
    with open(gwas_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Read the header from the GWAS file and append new column names
        header = infile.readline().strip()
        outfile.write(header + "\tNUM_SAMPLE\tAT\n")  # Add new columns to the header

        # Iterate through the GWAS file and VCF variants
        for line, variant in zip(infile, vcf):
            parts = line.strip().split("\t")

            # Get genotypes from the VCF for the current variant
            sample_genotypes = variant.genotypes  # List of [allele1, allele2, phased]

            # Count the number of samples with non-zero alleles
            num_sample = sum(sum(gt[:2]) > 0 for gt in sample_genotypes if gt[0] != -1)

            # Extract the 'AT' field from the VCF INFO column
            allele_type = variant.INFO.get("AT", "NA")  # Use "NA" if AT is not found

            # Write the updated line to the output file
            outfile.write(line.strip() + f"\t{num_sample}\t{allele_type}\n")

    print(f"Updated GWAS file with NUM_SAMPLE and AT columns saved to {output_file}")

except FileNotFoundError:
    print(f"Error: File {gwas_file} or {vcf_file} not found.")
    sys.exit(1)
except Exception as e:
    print(f"An error occurred: {e}")
    sys.exit(1)
