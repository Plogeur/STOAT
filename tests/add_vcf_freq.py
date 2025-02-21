import pandas as pd
import re
import sys
from cyvcf2 import VCF # type: ignore

def split_snarl(input_str):
    # Split the string and filter empty elements, then convert to integers
    return [str(num) for num in re.split(r'[><]', input_str) if num]

def parse_snarl_path_file_dict(path_file:str, freq_file:str) -> dict :
    
    snarl_paths = {}
    df = pd.read_csv(path_file, sep='\t', dtype=str)
    df['paths'] = df['paths'].str.split(',')
    for paths in df['paths']:
        first_node = int(split_snarl(paths[0])[0])
        # all case : SNP, INS, DEL, COMPLEX
        snarl_paths[first_node] = [paths, (-1, -1)]

    df = pd.read_csv(freq_file, sep='\t')

    # Iterate in pairs and calculate the differences
    for i in range(0, len(df), 4):
        row_1 = df.iloc[i]
        row_2 = df.iloc[i + 1]
        row_3 = df.iloc[i + 2]
        row_4 = df.iloc[i + 3]
       
        freq_diff_1 = abs(row_1['freq'] - row_2['freq'])
        freq_diff_2 = abs(row_3['freq'] - row_4['freq'])

        paths = snarl_paths[int(row_1['start_node'])][0]
        second_node_path1 = int(split_snarl(paths[0])[1])
        second_node_path2 = int(split_snarl(paths[1])[1])

        # Define the order of freq by paths
        if second_node_path1 == int(row_1["next_node"]) :
            snarl_paths[int(row_1['start_node'])][1] = (f"{freq_diff_1:.4}",f"{freq_diff_2:.4}")

        elif second_node_path2 == int(row_1["next_node"]) :
            snarl_paths[int(row_1['start_node'])][1] = (f"{freq_diff_2:.4}",f"{freq_diff_1:.4}")

        else :
            raise ValueError("second_node_path1 and second_node_path2 not matched")

    return snarl_paths

def modify_vcf(input_vcf:str, snarl_paths:dict, output_vcf="tests/modified_vcf.vcf"):
    with open(input_vcf, 'r') as infile, open(output_vcf, 'w') as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)  # Write header lines unchanged
            else:
                fields = line.strip().split("\t")
                first_node_snarl_id = int(split_snarl(fields[2])[0])
                freq_1, freq_2 = snarl_paths[first_node_snarl_id][1]
                info_field = fields[7]
                info_field += f";FR={freq_1},{freq_2}"  # Append new INFO field
                fields[7] = info_field
                outfile.write("\t".join(fields) + "\n")

vcf_file = sys.argv[1]
path_file = sys.argv[2]
freq_file = sys.argv[3]
output_vcf = sys.argv[4]

snarl_paths = parse_snarl_path_file_dict(path_file, freq_file)
modify_vcf(vcf_file, snarl_paths, output_vcf)

# python3 tests/add_vcf_freq.py tests/simulation/binary_data/pg.snarls.freq.tsv tests/simulation/binary_data/list_snarl_paths.tsv tests/simulation/binary_data/merged_output.vcf