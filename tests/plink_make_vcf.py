import sys
from cyvcf2 import VCF

vcf_file = sys.argv[1]
gwas_file = sys.argv[2]
output_file = "updated_" + gwas_file

# Open the VCF file using cyvcf2
vcf = VCF(vcf_file)
