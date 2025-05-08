import random
import argparse

def generate_binary_pheno(num_samples, output_pheno):
    with open(output_pheno, "w") as f:
        f.write("FID\tIID\tPHENO\n")
        phenos = [1] * (num_samples // 2) + [2] * (num_samples - num_samples // 2)
        for i in range(num_samples):
            fid = iid = f"sample{i+1}"
            f.write(f"{fid}\t{iid}\t{phenos[i]}\n")

def generate_quantitative_pheno(num_samples, output_pheno):
    with open(output_pheno, "w") as f:
        f.write("FID\tIID\tPHENO\n")
        for i in range(num_samples):
            pheno = random.uniform(-1, 1)
            fid = iid = f"sample{i+1}"
            f.write(f"{fid}\t{iid}\t{pheno}\n")

def generate_qtl_gene_position(num_samples, output_qtl, output_gene_position, number_genes=100):

    chroms = [str(i+1) for i in range(10)]
    list_of_genes = [f"gene_{i}" for i in range(number_genes)]
    with open(output_gene_position, "w") as f:
        f.write("gene_name\tchr\tstart_pos\tend_pos\n")
        for i in range(0, number_genes):
            start_pos = (i + 1) * 100
            end_pos = start_pos + 1000
            f.write(f"{list_of_genes[i]}\t{chroms[i//10]}\t{start_pos}\t{end_pos}\n")

    with open(output_qtl, "w") as f:
        list_of_samples = "\t".join([f"sample{i+1}" for i in range(num_samples)])
        f.write(f"gene_name\t{list_of_samples}\n")
        for i in range(0, number_genes):
            gene_expression = "\t".join([random.uniform(0, 10) for _ in range(num_samples)])
            f.write(f"{list_of_genes[i]}\t{gene_expression}\n")

def generate_covar(num_samples, output_covar):
    age_list = [age for age in range(0, 100)]
    with open(output_covar, "w") as f:
        f.write("FID\tIID\tSEX\tAGE\tPC1\tPC2\tPC3\n")
        for i in range(num_samples):
            fid = iid = f"sample{i+1}"
            age = random.choice(age_list)
            sex = 1 if i % 2 == 0 else 2 
            pc1 = random.uniform(-1, 1)
            pc2 = random.uniform(-1, 1)
            pc3 = random.uniform(-1, 1)
            f.write(f"{fid}\t{iid}\t{sex}\t{age}\t{pc1}\t{pc2}\t{pc3}\n")

def random_genotype():
    genotype = ['0/0', '1/0', '0/1', '1/1']
    return random.choice(genotype)

def generate_paths(start_node):
    # Each path is 3 nodes long: ref and alt share first and last node, differ in the middle
    ref_path = f">{str(start_node)}>{str(start_node + 2)}>{str(start_node + 3)}"
    alt_path = f">{str(start_node)}>{str(start_node + 1)}>{str(start_node + 3)}"
    paths = f"{ref_path},{alt_path}"
    return paths, start_node + 4

def generate_vcf_and_paths(num_samples, num_variants, vcf_file, paths_file):
    with open(vcf_file, "w") as vcf, open(paths_file, "w") as paths:
        # VCF header
        vcf.write("##fileformat=VCFv4.2\n")
        vcf.write("##INFO=<ID=AT,Number=1,Type=String,Description=\"Graph path with ref and alt\">\n")
        vcf.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
        
        # Add contig headers for chromosomes 1–10
        for chrom_id in range(1, 11):
            vcf.write(f"##contig=<ID={chrom_id}>\n")
        
        # Column header line
        vcf.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t")
        vcf.write("\t".join([f"sample{i+1}" for i in range(num_samples)]) + "\n")

        # TSV header
        paths.write("chr\tpos\tsnarl\tpaths\ttype\n")

        current_node = 2
        variants_per_chrom = num_variants // 10
        chroms = [str(i+1) for i in range(10)]
        ref = "A"
        alt = "T"
        fmt = "GT"

        for chrom in chroms:
            for i in range(variants_per_chrom):
                pos = (i + 1) * 10
                all_paths, current_node = generate_paths(current_node)
                at_info = f"AT={all_paths}"

                genotypes = '\t'.join([random_genotype() for _ in range(num_samples)])
                var_id = f"{chrom}_{i+1}"
                vcf.write(f"{chrom}\t{pos}\trs{var_id}\t{ref}\t{alt}\t.\tPASS\t{at_info}\t{fmt}\t{genotypes}\n")
                paths.write(f"{chrom}\t{pos}\trs{var_id}\t{all_paths}\t{ref},{alt}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate phenotype, VCF, and path info for simulated pangenome variants.")
    parser.add_argument("num_samples", type=int, default=100, help="Number of samples")
    parser.add_argument("num_variants", type=int, default=100000, help="Total number of variants (must be divisible by 10)")
    parser.add_argument("--vcf_file", default="variants.vcf", help="Output VCF file")
    parser.add_argument("--paths_file", default="paths_snarl.tsv", help="Output paths/snarl file")
    parser.add_argument("--covar_file", default="covar.tsv", help="Output covariate file")
    parser.add_argument("--qtl", default="qtl.tsv", help="Output qtl file")
    parser.add_argument("--gene_position", default="gene_position.tsv", help="Output gene position file")
    
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--binary", default="binary_phenotype.tsv", help="Binary phenotype")
    group.add_argument("--quantitative", default="quantitative_phenotype.tsv", help="Quantitative phenotype")

    args = parser.parse_args()

    if args.num_variants % 10 != 0:
        parser.error("Number of variants must be divisible by 10 (for 10 chromosomes).")

    if args.binary:
        generate_binary_pheno(args.num_samples, args.binary)

    elif args.quantitative:
        generate_quantitative_pheno(args.num_samples, args.quantitative)
    
    generate_qtl_gene_position(args.num_samples, args.num_variants, args.qtl, args.gene_position)
    generate_covar(args.num_samples, args.covar_file)
    generate_vcf_and_paths(args.num_samples, args.num_variants, args.vcf_file, args.paths_file)

    print(f"Files generated:\n- {args.pheno_file}\n- {args.vcf_file}\n- {args.paths_file}\n- {args.covar_file}")

# python3 tests/create_dataset.py 200 1000 --pheno_file data/simu/phenotypes.txt --vcf_file data/simu/variants.vcf --paths_file data/simu/paths_snarl.tsv --covar_file data/simu/covar.tsv
