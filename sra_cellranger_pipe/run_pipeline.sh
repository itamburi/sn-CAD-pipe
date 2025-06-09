#!/bin/bash
#SBATCH --job-name=snakemake_pipeline
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=logs/snakemake_%j.out
#SBATCH --error=logs/snakemake_%j.err

# Load modules needed for Snakemake and your tools
module load snakemake/6.0.0   # or whichever version is available
module load sra-tools/3.0.0
module load cellranger/8.0.1

# Run snakemake on the Snakefile in the current directory
snakemake --cores 8 --use-conda --printshellcmds --rerun-incomplete
