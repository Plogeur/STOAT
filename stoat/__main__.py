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
    parser.add_argument("-p", "--pg",type=utils.check_file, help='The input pangenome .pg file', required=False)
    parser.add_argument("-d", "--dist",type=utils.check_file, help='The input distance index .dist file', required=False)
    parser.add_argument("-n", "--name",type=utils.check_file, help='The input chromosome prefix reference file', required=False)
    parser.add_argument("-t", "--threshold",type=list_snarl_paths.check_threshold, help='Children threshold', required=False)
    parser.add_argument("-v", "--vcf",type=utils.check_format_vcf_file, help="Path to the merged VCF file (.vcf or .vcf.gz)", required=True)
    parser.add_argument("-l", "--listpath", type=utils.check_format_list_path, help="Path to the list paths", required=False)

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-b", "--binary", type=utils.check_format_pheno, help="Path to the binary group file (.txt or .tsv)")
    group.add_argument("-q", "--quantitative", type=utils.check_format_pheno, help="Path to the quantitative phenotype file (.txt or .tsv)")
    parser.add_argument("-c", "--covariate", type=utils.check_covariate_file, required=False, help="Path to the covariate file (.txt or .tsv)")
    parser.add_argument("-k", "--kinship", type=utils.check_kinship_prefix, required=False, help="Kinship prefix files (.grm.bin, .grm.id)")
    parser.add_argument("-g", "--gaf", action="store_true", required=False, help="Prepare binary gwas output to do gaf file + make gaf for the 10th significant paths")
    parser.add_argument("-o", "--output", type=str, required=False, help="Base path for the output directory")
    args = parser.parse_args()

    if args.quantitative and args.gaf:
        parser.error("The '--gaf' argument cannot be used with the '--quantitative' ('-q') argument.")

    if not args.listpath and (not args.pg or not args.dist) :
        parser.error("When --listpath (-l) is not provided, both -p and -d must be specified to compute it.")

    if (args.covariate and not args.kinship) or (args.kinship and not args.covariate):
        parser.error("Both --covariate (-c) and --kinship (-k) must be provided together.")

    if args.gaf and not args.pg:
        parser.error("The '--gaf' argument cannot be used without the '--pg' ('-p') argument.")

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

    if args.covariate :
        covar = utils.parse_covariate_file(args.covariate)
        utils.check_matching(covar, list_samples, args.covariate)
        kinship_matrix, kinship_ind = utils.parse_plink_grm(args.kinship)
        utils.check_matching(kinship_ind, list_samples, args.kinship)
    else :
        covar = None
        kinship_matrix = None

    if args.binary:
        logger.info("Parsing binary phenotype...")
        pheno = utils.parse_pheno_binary_file(args.binary)

    elif args.quantitative:
        logger.info("Parsing quantitative phenotype...")
        pheno = utils.parse_pheno_quantitatif_file(args.quantitative)
    
    utils.check_matching(pheno, list_samples, args.quantitative)

    if not args.listpath :
        
        if args.name :
            reference_chr = utils.parse_chr_reference(args.chr)
        else :
            reference_chr = {"ref":0}
            
        logger.info("Starting snarl path decomposition...")
        stree, pg, root, pp_overlay = list_snarl_paths.parse_graph_tree(args.pg, args.dist)
        snarls = list_snarl_paths.save_snarls(stree, root, pg, reference_chr, pp_overlay)
        logger.info(f"Total of snarls found : {len(snarls)}")
        logger.info("Saving snarl path decomposition...")
        output_snarl_path_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")
        output_snarl_path = os.path.join(output_dir, "snarl_paths.tsv")
        threshold = int(args.threshold) if args.threshold else 10 
        snarl_paths, paths_number_analysis = list_snarl_paths.loop_over_snarls_write(stree, snarls, pg, output_snarl_path, output_snarl_path_not_analyse, threshold)
        logger.info(f"Total of paths analyse : {paths_number_analysis}")
    else :
        if args.pg or args.dist : 
            logger.info("list snarls path are provided, .pg and .dist will be not use to make another list paths snarl")
        if args.gaf :
            pg = list_snarl_paths.parse_pg(args.pg)
        input_snarl_path = args.listpath
        snarl_paths, paths_number_analysis = utils.parse_snarl_path_file(input_snarl_path)
        logger.info(f"Total of snarls found : {paths_number_analysis}")

    vcf_object = snarl_analyser.SnarlProcessor(args.vcf, list_samples)
    logger.info("Starting fill matrix...")
    vcf_object.fill_matrix()

    # Handle Binary Analysis
    if args.binary:
        gaf = True if args.gaf else False
        output_snarl = os.path.join(output_dir, "binary_analysis.tsv")
        logger.info("Binary table creation...")
        vcf_object.binary_table(snarl_paths, pheno, kinship_matrix, covar, gaf, output_snarl)

        output_manh = os.path.join(output_dir, "manhattan_plot_binary.png")
        output_qq = os.path.join(output_dir, "qq_plot_binary.png")
        output_significative = os.path.join(output_dir, "top_variant_binary.tsv")
        logger.info("Binary p-value analysis...")
        p_value_analysis.significative_snarl_binary(output_snarl, output_significative)
        p_value_analysis.qq_plot_binary(output_snarl, output_qq)
        p_value_analysis.plot_manhattan_binary(output_snarl, output_manh)
        if gaf :
            output_gaf = os.path.join(output_dir, "group_paths.gaf")
            logger.info("GAF creation...")
            gaf_creator.parse_input_file(output_significative, snarl_paths, pg, output_gaf)
 
    # Handle Quantitative Analysis
    elif args.quantitative:
        output_file = os.path.join(output_dir, "quantitative_analysis.tsv")
        logger.info("Quantitative table creation...")
        vcf_object.quantitative_table(snarl_paths, pheno, kinship_matrix, covar, output_file)

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
Usage example:
    stoat -p ../droso_data/fly/fly.pg -d ../droso_data/fly/fly.dist -v ../droso_data/pangenome.dm6.vcf \
    -r ../droso_data/fly/fly.deconstruct.vcf -q ../droso_data/pangenome_phenotype.tsv -o output
 
Usage test:
    stoat -p tests/simulation/binary_data/pg.full.pg -d tests/simulation/binary_data/pg.dist -v tests/simulation/binary_data/merged_output.vcf \
    -b tests/simulation/binary_data/phenotype.tsv -o output

    stoat -p tests/simulation/quantitative_data/pg.full.pg -d tests/simulation/quantitative_data/pg.dist -v tests/simulation/quantitative_data/merged_output.vcf \
    -q tests/simulation/quantitative_data/phenotype.tsv -o output
"""
