import sys

def gwas_generator(gwas_file):
    with open(gwas_file, 'r') as gwas:
        gwas.readline()  # Skip header line
        for line in gwas:
            yield line.strip().split()

# Input files
plink_gwas_file = sys.argv[1]  # PLINK GWAS results file
vcf_file = sys.argv[2]  # Input VCF file
output_vcf = sys.argv[3]  # Output VCF file with updated INFO column

# Process the VCF file and append PV info
gwas_iter = gwas_generator(plink_gwas_file)

with open(vcf_file, 'r') as vcf, open(output_vcf, 'w') as out_vcf:
    for line in vcf:
        if line.startswith("##"):
            out_vcf.write(line)
        elif line.startswith("#CHROM"):
            out_vcf.write("##INFO=<ID=PV,Number=1,Type=Float,Description=\"P-value from GWAS\">\n")
            out_vcf.write(line)
        else:
            gwas_line = next(gwas_iter)
            pval = gwas_line[-2] # Last last column... is the p-value
            if pval == "NA":
                pval = "1"

            if float(pval) > 1 :
                raise ValueError("P-value is greater than 1. Please check the GWAS file.")
            
            fields = line.strip().split("\t")
            info_field = fields[7]
            updated_info = info_field + ";PV=" + pval if info_field != "." else "PV=" + pval
            fields[7] = updated_info
            out_vcf.write("\t".join(fields) + "\n")

print("Annotated VCF file saved as", output_vcf)
