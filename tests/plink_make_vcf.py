import sys
from cyvcf2 import VCF
import pandas as pd 
import re

vcf_file = sys.argv[1]
path_file = sys.argv[2]
freq_file = sys.argv[3]

def split_snarl(input_str):
    # Split the string and filter empty elements, then convert to integers
    return [str(num) for num in re.split(r'[><]', input_str) if num]

def parse_snarl_path(path_file:str) -> dict :
    
    snarl_paths = {}
    df = pd.read_csv(path_file, sep='\t', dtype=str)
    for seq, pos in zip(df['type'], df['pos']):
        snarl_paths[pos] = seq.split(',')
    return snarl_paths

dict_snarl_pos = parse_snarl_path(path_file)

for variant in VCF(vcf_file):
    genotypes = variant.genotypes  # Extract genotypes once per variant
    