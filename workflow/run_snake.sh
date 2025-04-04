#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --mem=500M

module load bioinfo/Snakemake/7.20.0
snakemake -s Snakefile --profile . --rerun-incomplete --keep-going --dry-run


