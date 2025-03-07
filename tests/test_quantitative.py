import os
from pathlib import Path
from PIL import Image, ImageChops

# Import necessary modules
from stoat_node import snarl_analyser
from stoat_node import list_snarl_paths
from stoat_node import utils
from stoat_node import p_value_analysis

def test_snarl_simulation_analyser():
    vcf_file = "tests/simulation/quantitative_data/merged_output.vcf"
    phenotype_file = "tests/simulation/quantitative_data/phenotype.tsv"
    dist_file = "tests/simulation/quantitative_data/pg.dist"
    pg_file = "tests/simulation/quantitative_data/pg.full.pg"
    output_dir = Path("tests/quantitative_tests_output")

    expected_output_snarl_path = "tests/simulation/quantitative_data/snarl_paths.tsv"
    expected_output = "tests/simulation/quantitative_data/quantitative_analysis.tsv"
    expected_significative = "tests/simulation/quantitative_data/top_variant_quantitative.tsv"
    expected_manh = "tests/simulation/quantitative_data/manhattan_plot_quantitative.png"
    expected_qq = "tests/simulation/quantitative_data/qq_plot_quantitative.png"

    # Perform test logic
    list_samples = utils.parsing_samples_vcf(vcf_file)
    vcf_object = snarl_analyser.SnarlProcessor(vcf_file, list_samples)
    vcf_object.fill_matrix()
    reference_chr = {"ref":0}
    stree, pg, root, pp_overlay = list_snarl_paths.parse_graph_tree(pg_file, dist_file)
    snarls = list_snarl_paths.save_snarls(stree, root, pg, reference_chr, pp_overlay)

    output_snarl_path_not_analyse = os.path.join(output_dir, "snarl_not_analyse.tsv")
    output_snarl_path = os.path.join(output_dir, "snarl_paths.tsv")    
    snarl_paths = list_snarl_paths.loop_over_snarls_write(stree, snarls, pg, output_snarl_path, output_snarl_path_not_analyse, children_treshold=10)[0]

    assert os.path.exists(output_snarl_path), f"Output file {output_snarl_path} was not created."

    with open(expected_output_snarl_path, 'r') as expected_snarl_path:
        expected_content_snarl_path = expected_snarl_path.read()

    with open(output_snarl_path, 'r') as output_snarl_path:
        output_content_output_snarl_path = output_snarl_path.read()

    assert output_content_output_snarl_path == expected_content_snarl_path, "The output content does not match the expected content."

    os.makedirs(output_dir, exist_ok=True)
    output = os.path.join(output_dir, "quantitative_test.assoc.tsv")

    quantitative_pheno = utils.parse_pheno_quantitatif_file(phenotype_file)
    vcf_object.quantitative_table(snarl_paths, quantitative_pheno, output=output)
    # write_position.write_pos_snarl(vcf_file, output, "quantitatif")

    assert os.path.exists(output), f"Output file {output} was not created."

    with open(expected_output, 'r') as expected_file:
        expected_content = expected_file.read()

    with open(output, 'r') as output_file:
        output_content = output_file.read()

    assert output_content == expected_content, "The output content does not match the expected content."

    output_significative = os.path.join(output_dir, "top_variant_quantitative.tsv")
    p_value_analysis.significative_snarl_quantitatif(output, output_significative)

    with open(expected_significative, 'r') as output_file:
        expected_signif = output_file.read()

    with open(output_significative, 'r') as output_file:
        output_content_signif = output_file.read()

    assert output_content_signif == expected_signif, "The output content does not match the expected content."

    output_manh = os.path.join(output_dir, "manhattan_plot_quantitative.png")
    output_qq = os.path.join(output_dir, "qq_plot_quantitative.png")
    p_value_analysis.qq_plot_quantitatif(output, output_qq)
    p_value_analysis.plot_manhattan_quantitatif(output, output_manh)

    # Open the images
    img_manh = Image.open(output_manh)
    img_qq = Image.open(output_qq)

    expected_img_manh = Image.open(expected_manh)
    expected_img_qq = Image.open(expected_qq)

    # Compare the images
    diff_manh = ImageChops.difference(img_manh, expected_img_manh)
    diff_qq = ImageChops.difference(img_qq, expected_img_qq)

    assert diff_manh.getbbox() is None, "The manhanttan expected and manhanttan content is differente."
    assert diff_qq.getbbox() is None, "The qq expected and manhanttan content is differente."

# pytest tests/test_quantitative.py
