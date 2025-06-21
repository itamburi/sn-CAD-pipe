#!/bin/bash
#SBATCH --job-name=snakemake_workflow
#SBATCH --output=snakemake_%j.out
#SBATCH --error=snakemake_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=12:00:00
#SBATCH --account=mseldin_lab  # Use your SLURM account if required

# Optional: activate conda or module env
# module load snakemake
# source activate your_snakemake_env
module load miniconda3/24.9.2
source /opt/apps/miniconda3/24.9.2/etc/profile.d/conda.sh
conda activate snakemake_env

# Move to the directory where the Snakefile is
cd /dfs6/pub/itamburi/CAD_snrnaseq

# Run Snakemake
snakemake \
    --snakefile Snakefile \
    --jobs 20 \
    --use-conda \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    --printshellcmds \
    --cores $SLURM_CPUS_PER_TASK
