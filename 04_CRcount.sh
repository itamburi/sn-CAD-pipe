#!/bin/sh
#SBATCH --job-name=cellranger_cnt    # Job name                  
#SBATCH -A mseldin_lab
#SBATCH --array=1-18
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.err    # error based on jobname

module load cellranger/8.0.1  

# guide: https://bioinformatics-core-shared-training.github.io/UnivCambridge_ScRnaSeq_Nov2021/Markdowns/03_CellRanger.html
wd="/dfs6/pub/itamburi/CAD_snrnaseq"
ref="${wd}/cellranger_HS_ref"
fq="${wd}/fq"

# create prefixes.txt with: ls | sed 's/_[0-9]\.fastq\.gz//' | sort | uniq > prefixes.txt
prefix=`cat ${wd}/prefixes.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n 1`

# need to rename the fq files from SRR to follow the 10X format...
cellranger count --id=$prefix \
                 --transcriptome=$ref \
                 --fastqs=$fq \
                 --sample=$prefix \
                 --localcores=8 \
                 --localmem=24 \
                 --create-bam false
