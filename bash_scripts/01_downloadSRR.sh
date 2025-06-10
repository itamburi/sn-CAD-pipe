#!/bin/sh
#SBATCH --job-name=fastq-dump    # Job name                  
#SBATCH --cpus-per-task 4
#SBATCH --array=1-18
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.err    # error based on jobname

#downloading GEO accession GSE131780
# CAD snRNAseq from Wirka 2019
# SRR9130237-SRR9130254, 18 samples


module load sra-tools/3.0.0  

START_NUM=9130236 #n-1
ACCESSION_NUM=$((START_NUM + SLURM_ARRAY_TASK_ID))
SRA_ID="SRR${ACCESSION_NUM}"

#had to run vdb-config -i and set repository location to /dfs6/pub/itamburi/CAD_snrnaseq/sra
fastq-dump --split-files --gzip $SRA_ID
