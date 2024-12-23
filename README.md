# Slink

Slink is a C++ implementation of the STOAT tool, designed to create PLINK formats file using pangenome paths for advanced Genome-Wide Association Studies (GWAS) focusing on snarl structures within pangenome graphs.

## Installation

````bash
git clone https://github.com/Plogeur/STOAT.git
cd STOAT
pip install -r requirements.txt

# install bdsg version > 3.0.0 
# DO NOT use pip install bdsg cause version 3.0.0 instead do :
git clone --recursive https://github.com/vgteam/libbdsg.git
cd libbdsg
pip install .

# see more installation information on bdsg github
````

## Input format file

Required files :
- pg_file : Pangenome graph file, formats accepted: .pg or .xg.
- dist_file : Distance file generated with vg dist, format: .dist.
- vcf : Merged VCF file, created using bcftools merge, formats: .vcf or .vcf.gz.
- pheno : Three-column file with FID (family/sample name), IID (sample name), and PHENO (integer/float). Format: .txt or .tsv (tab-separated).

Optional file : 
- list_paths : Two-column file containing snarl names and the list of paths through the snarl's netgraph, separated by tabs. Format: .txt or .tsv. (avoid to compute list_paths file)

## Usage

```bash
# Without list_paths file
./slink -p <pg_file> -g <dist_file> --vcf <vcf_path> --pheno <pheno> -o <output_dir>

# With list_paths file
./slink --snarl <list_paths> --vcf <vcf_path> --pheno <pheno> -o <output_dir>
```

## Output file 
The PLINK format is used in genetic studies for managing and analyzing genotype data. Below is an example of a basic PLINK output section, formatted as a `.bed`, `.bim`, and `.fam` file for binary data.

### BED File
The `.bed` file is a binary file containing genotype information. It does not have a readable text format, but its structure includes:
- A magic number header to identify the file.
- A genotype matrix where each individual's genotype data is stored in binary format.

### BIM File Example
The `.bim` file contains information about the SNARL in the SNP column.

```
CHR SNP          CM      BP      A1 A2
1   >123>2423    0.0     123456  A  G
1   >124>2425    0.0     123789  T  C
2   >125>2426    0.0     456789  G  T
```

- **CHR**: Chromosome number
- **SNP**: SNP identifier
- **CM**: Genetic distance (centiMorgans)
- **BP**: Base-pair position
- **A1**: Reference allele
- **A2**: Alternate allele

### FAM File Example
The `.fam` file contains information about individuals. Each row represents an individual.

```
FID IID  PID  MID  SEX  PHENOTYPE
F001 I001 0    0    1    2
F001 I002 0    0    2    1
F002 I003 0    0    1    2
```

- **FID**: Family ID
- **IID**: Individual ID
- **PID**: Paternal ID (0 if unknown)
- **MID**: Maternal ID (0 if unknown)
- **SEX**: Sex (1=Male, 2=Female, 0=Unknown)
- **PHENOTYPE**: Disease status (1=Unaffected, 2=Affected, -9=Unknown)

---

### Example of GWAS command using slink output

To analyze this data using PLINK, run the following command:
```bash
# Quantitative phenotype
plink --bfile slink_prefix_output --linear --out results

# Binary phenotype
plink --bfile slink_prefix_output --assoc --out results
```

