#!/bin/sh
#SBATCH --job-name=rename_sra    # Job name                  
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.err    # error based on jobname



# fastq-dump downloaded three fqs for each sample. e.g. SRR##_1, SRR##_2, SRR##_3
# by zcat SRR##_1.fq.gz | head, and checking with the followng ref we learn:
# https://davetang.org/muse/2018/06/06/10x-single-cell-bam-files/
# reads in SRR##_1 are 26 bases = cell barcode + UMI = R1
# reads in SRR##_2 are 98 bases = cDNA = R2
# reads in SRR##_3 are 8 bases = i7 sample index = I1

# need to rename to be compatible with cellranger's expected names format

wd="/dfs6/pub/itamburi/CAD_snrnaseq/fq"
cd $wd

# before and after:
# SRR#######_1.fastq.gz 
# <SampleName>_S<SampleNumber>_L00<Lane>_<Read>_001.fastq.gz
# seen here: https://bioinformatics-core-shared-training.github.io/UnivCambridge_ScRnaSeq_Nov2021/Markdowns/03_CellRanger.html

for fq in SRR*_*.fastq.gz; do
  base=$(echo $fq | cut -d'_' -f1)
  read=$(echo $fq | cut -d'_' -f2 | cut -d'.' -f1)
  case $read in
    1) suffix="R1";;
    2) suffix="R2";;
    3) suffix="I1";;
    *) echo "Unknown read: $read"; continue;;
  esac
  cp ${base}_${read}.fastq.gz ${base}_S1_L001_${suffix}_001.fastq.gz
done


# ... just delete the old names when complete
