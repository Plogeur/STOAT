#!/bin/bash

set -euo pipefail

mkdir -p output
sort -t$'\t' -k3,3n data/binary/snarl_paths.tsv > data/binary/snarl_paths.modify.tsv

plink --bfile output/genotype --pheno data/binary/binary.plink.phenotype.tsv \
    --pheno-name PHENO --assoc --allow-no-sex \
    --allow-extra-chr --out output/plink_binary

sed -i '' 's/[[:space:]]\{1,\}$//' output/plink_binary.assoc
awk '{{gsub(/^ +/, ""); gsub(/ +/, "\t")}}1' output/plink_binary.assoc > output/plink_binary.modify
sed -i '' -e 's/\tSNP\t/\tSNARL\t/' -e 's/\tBP\t/\tPOS\t/' -e 's/\tP\t/\tP_FISHER\t/' output/plink_binary.modify

sort -t$'\t' -k2,2n output/plink_binary.modify > output/plink_binary.modify2
mv output/plink_binary.modify2 output/plink_binary.modify

# Remove from snarl_paths the snarl that is remove with plink maf
# by removing the snarl that is not in plink_binary.modify in bash analysing line by line
bash filter_snarl.bash data/binary/snarl_paths.modify.tsv output/plink_binary.modify

python3 verify_truth.py --freq data/binary/pg.snarls.freq.tsv \
    --p_value output/plink_binary.modify --paths data/binary/snarl_paths.modify.tsv -b
