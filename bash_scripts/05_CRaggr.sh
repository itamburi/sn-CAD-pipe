#!/bin/sh
#SBATCH --job-name=cellranger_aggr    # Job name                  
#SBATCH -A mseldin_lab
#SBATCH --cpus-per-task 8
#SBATCH --mem=24GB
#SBATCH --output=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.out   # log based on jobname
#SBATCH --error=/dfs6/pub/itamburi/CAD_snrnaseq/%x_%j.err    # error based on jobname

module load cellranger/8.0.1  

wd="/dfs6/pub/itamburi/CAD_snrnaseq"

cd $wd

# make the csv file for aggr
#echo "sample_id,molecule_h5" > GSE131780_libraries.csv
#find . -type f -name "molecule_info.h5" -path "./SRR*" | while read filepath; do
#  dir=$(basename "$(dirname "$(dirname "$filepath")")")
#  fullpath=$(realpath "$filepath")
#  echo "$dir,$fullpath" >> GSE131780_libraries.csv
#done

cellranger aggr --id=GSE131780 --csv=GSE131780_libraries.csv --normalize=none
