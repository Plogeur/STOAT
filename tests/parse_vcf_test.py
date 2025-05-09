from cyvcf2 import VCF # type: ignore

vcf_path = "variants.vcf"
# Parse variant line by line
for variant in VCF(vcf_path):
    
    variant_path = variant.INFO.get('AT') # float
    #genotypes = variant.genotypes  # Extract genotypes once per variant
    print("variant_path : ", variant_path)