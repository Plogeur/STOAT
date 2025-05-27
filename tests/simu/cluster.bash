#!/bin/bash
#SBATCH -J simu
#SBATCH -o simu.out
#SBATCH -e simu.err
#SBATCH -t 2:00:00
#SBATCH -p unlimitq
#SBATCH -n 1
#SBATCH -c 16
#SBATCH --mem=10G
#SBATCH --mail-user=matis.alias-bagarre@inserm.fr
#SBATCH --mail-type=END,FAIL

module load bioinfo/Snakemake/7.32.4 bioinfo/vg/1.57.0

## run pipeline
snakemake --cores 16 --snakefile Snakefile_simulation --configfile config/config_binary.yaml
