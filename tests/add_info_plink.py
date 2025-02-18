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

# def create_dict_AT(vcf_file_ref, vcf_file_at):
#     vcf_dict = {}

#     # Process reference VCF file
#     for variant in VCF(vcf_file_ref):
#         sample_genotypes = variant.genotypes  # List of genotype lists
#         num_allele = sum(1 for sample in sample_genotypes for allele in sample[:2] if allele != -1)  # Count all alleles, including 0
#         var = variant.ID.split("_")[0]
#         if var in vcf_dict :
#             vcf_dict[var][0][variant.ID] = num_allele

#         else :
#             vcf_dict[variant.ID] = [{variant.ID:num_allele}, "NA"]
            
#     # Process AT values from another VCF
#     for variant in VCF(vcf_file_at):
#         at_value = variant.INFO.get('AT', "NA")  # Ensure missing values default to "NA"
#         if variant.ID in vcf_dict:
#             vcf_dict[variant.ID][1] = at_value

#     return vcf_dict

# def fill_gwas(gwas_file, vcf_dict):
#     with open(gwas_file, 'r') as infile, open(output_file, 'w') as outfile:
#         # Read the header from the GWAS file and append new column names
#         header = infile.readline().strip()
#         outfile.write(header + "\tNUM_SAMPLE\tAT\n") # Add new columns to the header

#         # Iterate through the GWAS file and VCF variants
#         for line in infile:
#             parts = line.strip().split("\t")

#             snarl = parts[1]
#             parts[1] = parts[1].split('_')[0]
#             dict_sample, allele = vcf_dict.get(parts[1], ["NA", "NA"])  # Corrected lookup
#             if dict_sample != "NA" :
#                 num_sample = dict_sample.get(snarl, "NA")
#             else :
#                 num_sample = "NA"
            
#             # Write the updated line to the output file
#             outfile.write("\t".join(parts) + f"\t{num_sample}\t{allele}\n")

#     print(f"Updated GWAS file with NUM_SAMPLE and AT columns saved to {output_file}")
