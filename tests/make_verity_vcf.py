import argparse
import re
import pandas as pd

def split_snarl(input_str):
    # Split the string and filter empty elements, then convert to integers
    return [str(num) for num in re.split(r'[><]', input_str) if num]

def modify_path(path:str) -> str:
    match = re.findall(r"[<>]\d+", path)
    return match[0] + match[-1]

def extract_vcf_header():
    header_lines = [
        '##fileformat=VCFv4.2\n',
        '##INFO=<ID=F,Number=1,Type=String,Description="Absolute value diff Frequency verity">\n',
        '##FILTER=<ID=PASS,Description="All filters passed">\n',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n',
        '##SAMPLE=<ID=SAMP0>\n'
    ]
    
    column_header = 'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMP0'
    return header_lines, column_header

def parse_snarl_path(path_file:str, freq_file:str) -> dict :
    
    snarl_paths = {}
    df = pd.read_csv(path_file, sep='\t', dtype=str)
    df['paths'] = df['paths'].str.split(',')
    df['type'] = df['type'].str.split(',')
    for chr, pos, paths, seqs in zip(df['chr'], df['pos'], df['paths'], df['type']):
        id_node = int(split_snarl(paths[0])[0]) # get the first node id from a paths
        path_modify = modify_path(paths[0])
        # all case : SNP, INS, DEL, COMPLEX
        snarl_paths[id_node] = [chr, path_modify, seqs, pos, -1] # sequence, position, abs(freq group 0 - freq group 1)

    df = pd.read_csv(freq_file, sep='\t')
    for i in range(0, len(df) - 1, 2):
        row_1 = df.iloc[i]
        row_2 = df.iloc[i + 1]
        freq_diff = abs(row_1['freq'] - row_2['freq'])
        id_node = row_1['start_node']
        snarl_paths[id_node][4] = f"{freq_diff:.4}"
            
    return snarl_paths

def create_vcf_verity(snarl_paths_file, freq_file, output_vcf):
    # Extract the header lines, column header, and number of samples
    header_lines, column_header = extract_vcf_header()
    dict_node = parse_snarl_path(snarl_paths_file, freq_file)

    # Open the output VCF for writing
    with open(output_vcf, 'w') as out_vcf:
        out_vcf.writelines(header_lines)
        out_vcf.write(column_header + '\n')

        for chr, snarl_id, seq, pos, freq_diff in dict_node.values():
            ref, alt = seq
            qual = filter = '.'
            info = f'F={freq_diff}'
            gt = './.'
            # remove the last \t
            vcf_line = f"{chr}\t{pos}\t{snarl_id}\t{ref}\t{alt}\t{qual}\t{filter}\t{info}\t{gt}\n"
            out_vcf.write(vcf_line)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GWAS output to a VCF file using an input VCF header.")
    parser.add_argument('-l', '--snarl_paths_file', type=str, required=True, help="Path to the input snarl paths file (from stoat).")
    parser.add_argument('-f', '--freq_file', type=str, required=True, help="Path to the input frequency verity simulation file (from simulation)")
    parser.add_argument('-o', '--output_vcf', type=str, required=True, help="Path to the output VCF file.")

    args = parser.parse_args()

    # Run the function with the parsed arguments
    create_vcf_verity(args.snarl_paths_file, args.freq_file, args.output_vcf)

# python3 tests/make_verity_vcf.py -l output/test/list_snarl_paths.tsv -f tests/simulation/binary_data/pg.snarls.freq.tsv -o tests/binary_tests_output/test_verity.vcf
