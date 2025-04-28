# STOAT (Snarl Tree Orchestrated Association Test)

<p align="center">
    <a href="https://www.python.org/downloads/release/python-3100/"><img src="https://img.shields.io/badge/Python-3.10-blue.svg"></a>
    <a href="https://github.com/vgteam/libbdsg/releases/tag/v0.3"><img src="https://img.shields.io/badge/bdsg-0.3-green.svg"></a>
</p>

<!-- 
It will release one days, It will !!!
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/jmonlong/sveval)](https://github.com/jmonlong/sveval/releases/latest)
[![Docker Repository on Quay](https://quay.io/repository/jmonlong/sveval/status "Docker Repository on Quay")](https://quay.io/repository/jmonlong/sveval) 
-->

<img src="pictures/logo.png" width="150">

## Project Overview

<img src="pictures/stoat_rep.png">


## Dependency

Manual installation : 

- jansson 
- libbdsg
- htslib
- eigen3
- boost
- Catch2 v3

Or you can use the compiled unbuntu version provided.

## Building

```bash
git clone --branch stoat_cxx https://github.com/Plogeur/STOAT.git
cd stoat_cxx

mkdir build && cd build
cmake .. && make -j 4

# ./stoat
# ./unit_tests
```  

STOAT is a specialized tool developed for conducting Genome-Wide Association Studies (GWAS) with a unique focus on snarl structures within pangenome graphs. Unlike traditional GWAS tools that analyze linear genome variants, STOAT processes VCF files to extract and analyze snarl regions—complex structural variations that capture nested and overlapping variant patterns within a pangenome. This approach allows for a more nuanced understanding of genetic variations in diverse populations and complex traits.

STOAT supports both binary and quantitative phenotypes:

- For binary phenotypes (e.g., case vs control studies), it utilizes chi-squared tests and Fisher’s exact test to evaluate associations between phenotype groups and snarl variants, providing robust statistical validation even in cases of sparse data.

- For quantitative phenotypes (e.g., traits measured on a continuous scale), the tool employs linear regression models to assess the association between snarl structures and phenotype values, allowing for continuous trait mapping with greater precision.

## Input format file

Required files :
- pg : Pangenome graph file, formats accepted: .pg or .xg.
- dist : Distance file generated with vg dist, format: .dist.
- vcf pangenomique : Merged VCF file, created using `vg pipeline` and bcftools merge, formats: .vcf or .vcf.gz. (ex : `bcftools merge -m none -Oz -o test`)
- phenotype : phenotype file organise in three-column with FID (family/sample name), IID (sample name), and PHENO (integer/float). Format: .txt or .tsv (tab-separated).
- chromosome : Txt file that containt the reference chromosome haplotype name in the pangenome graph. Format: .txt or .tsv.

Optional file : 
- paths : Snarl decoposition stoat output, Two-column file containing snarl names and the list of paths through the snarl's netgraph, separated by tabs. Format: .txt or .tsv.
- kinship : Kinship matrix file use in LMM analysis
- covariate : Covariate file. Format: .txt or .tsv.
- position gene : Gene position. Format: .txt or .tsv.

The VCF pangenomique is a VCF merged from a pangenomique mapping+calling, we recommand to use the [vg pipeline](https://github.com/vgteam/vg_snakemake)

VCF file :
```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	samp_g0_0	samp_g0_10	samp_g0_11
Chr1	411	>1>9	CGATTATGGA	C,CGATTA,CGATT,CGA	396.121	PASS	LV=0;DP=1241;AT=>1>2>4>5>6>8>9,>1>9,>1>2>3>5>6>8>9,>1>2>3>5>7>8>9,>1>2>4>5>7>8>9	GT:DP:AD:GL:GQ:GP:XD:MAD	./1	0/2	1/3
```

Phenotype file :
```
FID	IID	PHENO
samp_g0_0	samp_g0_0	1
samp_g0_1	samp_g0_1	1
```

Chromosome reference file :
```
Chr1
Chr2
Chr3
```

Covariate file : 
```
IID	SEX	CP1	CP2	CP42
samp_g0_0	1	1.2	562.25	42.25
samp_g0_1	0	5.2	359.25	65.24
```

## Usage

Use `stoat tool` if you want to launch the full tool at once, starting from snarl path identification (identifying the multiple paths that can be taken by a sample based on the pangenome graph) and ending with the results plots (Manhattan plot and QQ plot).

- Run full tool :
```bash
# binary trait
stoat -p <pg.pg> -d <dist.dist> -v <vcf.vcf.gz> -b <phenotype.txt> -o output

# quantative trait
stoat -p <pg.pg> -d <dist.dist> -v <vcf.vcf.gz> -q <phenotype.txt> -o output
```

Explanation of all options:
```bash 
-v, --vcf <path>            Path to the VCF file (.vcf or .vcf.gz)
-s, --snarl <path>          Path to the snarl file (.txt or .tsv)
-p, --pg <path>             Path to the pg file (.pg)
-d, --dist <path>           Path to the dist file (.dist)
-r, --chr <path>            Path to the chromosome reference file (.txt)
--children <int>            Max number of children for a snarl in the snarl decomposition process (default = 50)
--path-length <int>         Max length for a path snarl in the snarl decomposition process (default = 10 000)
-b, --binary <path>         Path to the binary group file (.txt or .tsv)
-g, --gaf                   Make a GAF file from the GWAS analysis
-q, --quantitative <path>   Path to the quantitative phenotype file (.txt or .tsv)
--covariate <path>          Path to the covariate file (.txt or .tsv)
--covar-name <string>       Covariate column name used in the gwas analyse
-e, --eqtl <path>           Path to the Expression Quantitative Trait Loci file (.txt or .tsv)
--gene-position <path>      Path to the Gene position file (.txt or .tsv)
-k, --kinship <path>        Path to the kinship matrix file (.txt or .tsv)
--make-bed                  Create a plink format files (.bed, .bim, .fam)
--table-threshold <int>     The N p-value digits threshold to use for plotting regression data file (defauld : 5 <=> 10-5)
--maf                       Add a maf (Minimum allele frequency) thresold (defauld : 0.01)
-o, --output <name>         Output dir name
-t, --thread <int>          Number of threads
-h, --help                  Print this help message;

```

- Example of usage : 

```bash
# decompose pangenome
./stoat -p <pg.pg> -d <dist.dist> -o <paths.txt>

# binary trait with list_path already computed and gaf creation 
./stoat -l <paths.txt> -v <vcf.vcf.gz> -r <ref.vcf.gz> -b <phenotype.txt> -o output.tsv

# quantitative trait with list_path already computed
./stoat -l <paths.txt> -v <vcf.vcf.gz> -r <ref.vcf.gz> -q <phenotype.txt> -o output.tsv
```

## Output

| Column Name       | Description                                                                                   |
|-------------------|-----------------------------------------------------------------------------------------------|
| **CHR**           | Chromosome name where the variation occurs.                                                   |
| **POS**           | Position of the snarl within the chromosome.                                                  |
| **SNARL**         | Identifier for the variant, snarl name/id                                                     |
| **TYPE**          | Type of genetic variation (e.g., SNP, INS, DEL).                                              |
| **P_FISHER**      | P-value calculated using Fisher's exact test (binary analysis).                               |
| **P_CHI2**        | P-value calculated using the Chi-squared test (binary analysis).                              |
| **ALLELE_NUM**    | Total number of alleles that pass in this snarl.                                              |
| **MIN_ROW_INDEX** | Minimum group of samples that pass through one path of the snarl. (binary analysis).          |
| **NUM_COLUM**     | Number of paths in the snarl. (binary analysis).                                              |
| **INTER_GROUP**   | Sum of the minimum samples that pass through each path. (binary analysis).                    |
| **AVERAGE**       | Average number of total samples passing through this snarl, divided by the number of paths. (binary analysis).  |
| **P**             | P-value calculated using linear regression (quantitative analysis).                           |
| **RSQUARED**      | R-squared value, proportion of variance explained by the model (quantitative analysis).       |
| **SE**            | Mean Standard error, estimatation coefficients of all paths in a snarl (quantitative analysis). |
| **BETA**          | Mean Beta coefficients, estimatation effect sizes of the prediction of all paths in a snarl (quantitative analysis). |

### Example of Output:

Below is an example of the output for a binary phenotype analysis:

```bash
CHR POS SNARL           TYPE  REF ALT   P_FISHER  P_CHI2  ALLELE_NUM  MIN_ROW_INDEX NUM_COLUM INTER_GROUP AVERAGE
1   12  5262721_5262719 SNP   A   T     0.4635    0.5182  286        2             137       46          143.0
1   15  5262719_5262717 INS   A   ATT   0.8062    0.8747  286        2             141       34          143.0
1   18  5262717_5262714 DEL   AA  T     0.2120    0.2363  286        2             134       32          143.0
```

Below is an example of the output for a quantitative phenotype analysis (-q option) :

```bash
CHR	POS	SNARL	TYPE	REF	ALT	RSQUARED	BETA	SE	P
1	12	5262721_5262719	SNP	A	T	8.3697e-01	1.3878e+01	6.5108e+00	4.0376e-01
1	15	5262719_5262717	INS	A	ATT	4.4237e-01	1.3238e+01	6.5345e+00	4.6574e-01
1	18	5262717_5262714	DEL	AA	T	6.3237e-01	1.6458e+01	6.6453e+00	4.7484e-01
1	19	5262717_5262714	COMPLEX	C	NA	4.2342e-01	2.3242e+01	5.3251e+00	1.3245e-01
```

## Visualization

### Manhattan and QQ plots 

STOAT will generated a manhattan and a QQ plot for binary and quantitatif analysis.

### SequenceTube

Use `--gaf` to geneate a GAF file and [sequenceTubeMap](https://github.com/vgteam/sequenceTubeMap) tool to visualize your gwas binary region results.

<p align="center">
<img src="pictures/seqTube.png" width="600">
</p>

Description : Color represente the different paths group (red : group 1 & blue : group 0) and opacity represente the number of samples in that paths (number of samples passing trought each paths % 60).
