import argparse

def extract_vcf_header(input_vcf):
    header_lines = [
        '##fileformat=VCFv4.2\n',
        '##FILTER=<ID=PASS,Description="All filters passed">\n',
        '##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">\n',
        '##INFO=<ID=P,Number=1,Type=String,Description="P-value">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    ]
    column_header = None
    sample_names = []
    chr_name = []

    # Open the VCF file and process line by line
    with open(input_vcf, 'r') as vcf:
        for line in vcf:
            if line.startswith("##contig=<ID="):
                chr_name.append(line.strip())
            
            elif line.startswith("#CHROM"):
                column_header = line.strip()  # Get the column header line
                columns = column_header.split('\t')
                sample_names = columns[9:]  # Extract sample names from column 10 onward
                break  # No need to read further; the header ends here

    # Add sample names to the header lines
    header_lines.extend([f"{chr}\n" for chr in chr_name])
    header_lines.extend([f"##SAMPLE=<ID={sample}>\n" for sample in sample_names])

    return header_lines, column_header

def create_vcf_from_gwas(gwas_file, input_vcf, info_vcf, output_vcf):
    # Extract the header lines, column header, and number of samples
    header_lines, column_header = extract_vcf_header(input_vcf)

    # Open the output VCF for writing
    with open(output_vcf, 'w') as out_vcf:
        out_vcf.writelines(header_lines)
        out_vcf.write(column_header + '\n')

        # Process GWAS file and generate VCF body
        with open(gwas_file, 'r') as gwas, open(info_vcf, 'r') as info:
            next(gwas)
            next(info)
            for line_gwas, line_info in zip(gwas, info):
                fields = line_gwas.strip().split('\t')
                information = line_info.strip().split('\t')
                #CHR POS SNARL TYPE	REF	ALT	P_FISHER
                chrom, pos_str, _, _, ref_str_brut, alt_str_brut = fields[:6]

                placeholder = str(information[2]).replace(',', '\t')
                qual = "."
                filter_field = "PASS"
                snarl_id = information[0]
                info_field = information[1]
                format_field = "GT"
                pos = pos_str.split(',')[0] # take only the first position
                list_ref = ref_str_brut.split(':')
                list_list_ref = [ref_str.split(',') for ref_str in list_ref]
                list_alt = alt_str_brut.split(':')
                list_list_alt = [alt_str.split(',') for alt_str in list_alt]

                dict_pos = {}
                dict_pos[pos] = [list_list_ref[0], list_list_alt[0]]

                for pos, [list_ref, list_alt] in dict_pos.items() :
                    for ref, alt in zip(list_ref, list_alt) :
                        # Create and write the VCF line
                        vcf_line = f"{chrom}\t{pos}\t{snarl_id}\t{ref}\t{alt}\t{qual}\t{filter_field}\t{info_field}\t{format_field}\t{placeholder}\n"
                        out_vcf.write(vcf_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GWAS output to a VCF file using an input VCF header.")
    parser.add_argument('-g', '--gwas', type=str, required=True, help="Path to the GWAS output file (TSV format).")
    parser.add_argument('-v', '--input_vcf', type=str, required=True, help="Path to the input VCF file (to extract header information).")
    parser.add_argument('-i', '--info', type=str, required=True, help="Path to the input VCF file (to extract header information).")
    parser.add_argument('-o', '--output_vcf', type=str, required=True, help="Path to the output VCF file.")

    args = parser.parse_args()

    # Run the function with the parsed arguments
    create_vcf_from_gwas(args.gwas, args.input_vcf, args.info, args.output_vcf)

# python3 stoat/vcf_maker.py -g output/run_20250117_105311/quantitative_analysis.tsv -v tests/simulation/quantitative_data/merged_output.vcf -i output/vcf_from_stoat.vcf -o output/quantitative_test_vcf.vcf
