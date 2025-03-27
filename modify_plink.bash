#!/bin/bash

set -euo pipefail

plink_file=$1

sed -i 's/[[:space:]]\{1,\}$//' $plink_file
awk '{{gsub(/^ +/, ""); gsub(/ +/, "\t")}}1' $plink_file > $plink_file.modify
sed -i -e 's/\tSNP\t/\tSNARL\t/' -e 's/\tBP\t/\tPOS\t/' -e 's/\tP\t/\tP_FISHER\t/' $plink_file.modify

sort -t$'\t' -k2,2n $plink_file.modify > $plink_file.modify2
mv $plink_file.modify2 $plink_file.modify

echo "Done modifying plink file"
