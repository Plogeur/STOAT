import argparse
from stoat import list_snarl_paths
from stoat import snarl_analyser
from stoat import utils
from stoat import p_value_analysis
from stoat import gaf_creator
import time
import logging
import os
import sys
from datetime import datetime

def main() : 

    # Argument Parsing
    parser = argparse.ArgumentParser(description="Run the Stoat GWAS analysis pipeline")
    parser.add_argument("-v", "--vcf",type=utils.check_format_vcf_file, help="Path to the merged VCF file (.vcf or .vcf.gz)", required=True)
    parser.add_argument("-l", "--listpath", type=str, help="Path to the list paths", required=False)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=utils.check_format_pheno, help="Path to the binary group file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=utils.check_format_pheno, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-o", "--output", type=str, required=False, help="Base path for the output directory")
    args = parser.parse_args()

    if not args.listpath :
        parser.error("When --listpath (-l) is not provided, both -p and -d must be specified to compute it.")

    # Generate unique output directory based on timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = os.path.join(args.output or "output", f"run_{timestamp}")
    os.makedirs(output_dir, exist_ok=True)

    # Configure logging with both console and file handlers
    log_file = os.path.join(output_dir, "run.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )

    logger = logging.getLogger(__name__)

    # Log the exact command used to launch the script
    command_line = " ".join(sys.argv)
    logger.info(f"Command: {command_line}")
    start_time = time.time()

    # Check vcf samples matching other files (pheno, covar)
    list_samples = utils.parsing_samples_vcf(args.vcf)

    if args.binary:
        logger.info("Parsing binary phenotype...")
        pheno = utils.parse_pheno_binary_file(args.binary)

    elif args.quantitative:
        logger.info("Parsing quantitative phenotype...")
        pheno = utils.parse_pheno_quantitatif_file(args.quantitative)
    
    utils.check_matching(pheno, list_samples, args.quantitative)

    input_snarl_path = args.listpath
    snarl_paths = utils.parse_snarl_path_file_dict(input_snarl_path)

    vcf_object = snarl_analyser.SnarlProcessor(args.vcf, list_samples)
    logger.info("Starting fill matrix...")
    vcf_object.fill_matrix()

    # Handle Binary Analysis
    if args.binary:
        gaf = False
        output_snarl = os.path.join(output_dir, "binary_analysis.tsv")
        logger.info("Binary table creation...")
        vcf_object.binary_table(snarl_paths, pheno, gaf, output_snarl)
        logger.info("Writing position...")

        output_manh = os.path.join(output_dir, "manhattan_plot_binary.png")
        output_qq = os.path.join(output_dir, "qq_plot_binary.png")
        output_significative = os.path.join(output_dir, "top_variant_binary.tsv")
        logger.info("Binary p-value analysis...")
        p_value_analysis.significative_snarl_binary(output_snarl, output_significative)
        p_value_analysis.qq_plot_binary(output_snarl, output_qq)
        p_value_analysis.plot_manhattan_binary(output_snarl, output_manh)

    # Handle Quantitative Analysis
    elif args.quantitative:
        output_file = os.path.join(output_dir, "quantitative_analysis.tsv")
        logger.info("Quantitative table creation...")
        vcf_object.quantitative_table(snarl_paths, pheno, output_file)
        logger.info("Writing position...")

        output_manh = os.path.join(output_dir, "manhattan_plot_quantitative.png")
        output_qq = os.path.join(output_dir, "qq_plot_quantitative.png")
        output_significative = os.path.join(output_dir, "top_variant_quantitative.tsv")
        logger.info("Quantitative p-value analysis...")
        p_value_analysis.significative_snarl_quantitatif(output_file, output_significative)
        p_value_analysis.qq_plot_quantitatif(output_file, output_qq)
        p_value_analysis.plot_manhattan_quantitatif(output_file, output_manh)

    logger.info(f"GWAS analysis completed in {time.time() - start_time:.2f} seconds.")
    logger.info(f"Output directory: {output_dir}")

if __name__ == "__main__":
    main()

"""
Usage test:
    stoat -l tests/simulation/binary_data/snarl_paths.tsv -v tests/simulation/binary_data/merged_output.vcf \
    -b tests/simulation/binary_data/phenotype.tsv -o output

    stoat -l tests/simulation/quantitative_data/snarl_paths.tsv -v tests/simulation/quantitative_data/merged_output.vcf \
    -q tests/simulation/quantitative_data/phenotype.tsv -o output
"""
