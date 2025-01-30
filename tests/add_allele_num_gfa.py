#!/usr/bin/env python3

import sys

# Check for correct number of arguments
if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <file1> <file2>")
    sys.exit(1)

file1 = sys.argv[1]
file2 = sys.argv[2]
output_file = "gfa2bin_plink_allele.tsv"

# Create a dictionary to store NODE -> ALLELE_NUM mapping
allele_dict = {}

# Read file2 and populate the dictionary
with open(file2, "r") as f2:
    header2 = f2.readline().strip().split("\t")  # Read header
    try:
        node_idx = header2.index("NODE")
        allele_idx = header2.index("ALLELE_NUM")
    except ValueError:
        print("Error: 'NODE' or 'ALLELE_NUM' column not found in file2.")
        sys.exit(1)

    for line in f2:
        parts = line.strip().split("\t")
        allele_dict[parts[node_idx]] = parts[allele_idx]  # Store NODE -> ALLELE_NUM mapping

# Process file1 and add ALLELE_NUM column
with open(file1, "r") as f1, open(output_file, "w") as out:
    header1 = f1.readline().strip().split("\t")  # Read header
    if "SNP" not in header1:
        print("Error: 'SNP' column not found in file1.")
        sys.exit(1)

    # Write new header with ALLELE_NUM
    out.write("\t".join(header1) + "\tALLELE_NUM\n")

    snp_idx = header1.index("SNP")

    for line in f1:
        parts = line.strip().split("\t")
        snp_value = parts[snp_idx]

        # Get ALLELE_NUM if SNP is in dictionary, otherwise use "NA"
        allele_num = allele_dict.get(snp_value, "NA")
        out.write("\t".join(parts) + f"\t{allele_num}\n")

print(f"Merged file saved as: {output_file}")
