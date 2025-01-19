from cyvcf2 import VCF
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def count_alleles_in_paths(vcf_file, paths):
    """
    Count the number of alleles passing through specified paths in a VCF file.

    Args:
        vcf_file (str): Path to the input VCF file.
        paths (list of str): List of paths to check, e.g., [">4242>4243>4245", ">4242>4244>4245"].

    Returns:
        dict: A dictionary with paths as keys and counts of valid alleles as values.
    """
    # Validate and decompose paths
    decomposed_paths = {
        path: [f">{'>'.join(path.split('>')[i:i + 2])}" for i in range(1, len(path.split(">")) - 1)]
        for path in paths
    }
    logger.info("Decomposed paths: %s", decomposed_paths)

    # Open the VCF file
    vcf = VCF(vcf_file)
    logger.info("Loaded VCF file: %s, Number of samples: %d", vcf_file, len(vcf.samples))

    # Initialize path counts using numpy for performance
    path_counts = {sub_path: np.zeros(len(vcf.samples) * 2, dtype=bool) for sub_paths in decomposed_paths.values() for sub_path in sub_paths}

    # Iterate over each variant in the VCF
    for variant in vcf:
        allele_paths = variant.INFO.get("AT", "").split(",")

        for sub_path in path_counts.keys():
            for path_index, allele_path in enumerate(allele_paths):
                if sub_path in allele_path:
                    for sample_index, genotype in enumerate(variant.genotypes):
                        path_counts[sub_path][sample_index * 2] |= (genotype[0] == path_index)
                        path_counts[sub_path][sample_index * 2 + 1] |= (genotype[1] == path_index)

    # Combine sub-path counts to compute full path counts
    result_counts = {}
    for path, sub_paths in decomposed_paths.items():
        sub_path1, sub_path2 = sub_paths
        valid_samples = path_counts[sub_path1] & path_counts[sub_path2]
        result_counts[path] = np.sum(valid_samples)

    logger.info("Final path counts: %s", result_counts)
    return result_counts

# Example usage
vcf_file = "tests/simulation/binary_data/norm_merged.vcf"
paths_to_check = [">4242>4243>4245", ">4242>4244>4245"]

path_counts = count_alleles_in_paths(vcf_file, paths_to_check)

for path, count in path_counts.items():
    print(f"Path {path}: {count} valid alleles")
