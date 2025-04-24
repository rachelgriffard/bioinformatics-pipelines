#!/bin/bash
#SBATCH --job-name=bw    # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=e617n596@kumc.edu     # Where to send mail	
#SBATCH --ntasks=8                  # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=8gb                    # Job memory request
#SBATCH --time=5-00:00:00             # Time limit days-hrs:min:sec
#SBATCH --array=1-4
#SBATCH --output=bw_%A_%a.log         # Standard output and error log 

pwd; hostname; date

ml conda
conda activate deeptools
ml samtools
module load R/4.2

echo "Running R script"

echo "$SLURM_ARRAY_TASK_ID"
 
Rscript 3_Get_BigWig.R $SLURM_ARRAY_TASK_ID
 
date
