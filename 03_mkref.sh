#!/bin/sh
#SBATCH --job-name=mkref    # Job name                  
#SBATCH --cpus-per-task 8
#SBATCH --mem-per-cpu=3G
#SBATCH --output=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.err    # error based on jobname



#curl -o /dfs6/pub/itamburi/CAD_snrnaseq/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
#	https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#curl -o /dfs6/pub/itamburi/CAD_snrnaseq/ensembl/Homo_sapiens.GRCh38.114.gtf.gz \
#	https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/Homo_sapiens.GRCh38.114.gtf.gz

module load cellranger/8.0.1  

fasta="/dfs6/pub/itamburi/CAD_snrnaseq/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf="/dfs6/pub/itamburi/CAD_snrnaseq/ensembl/Homo_sapiens.GRCh38.114.gtf.gz"

cellranger mkref \
  --nthreads=8 \
  --genome="cellranger_HS_ref" \
  --fasta=$fasta \
  --genes=$gtf
