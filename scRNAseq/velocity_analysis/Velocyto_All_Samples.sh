#!/bin/bash
#SBATCH --job-name=velocyto      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=    # Where to send mail	
#SBATCH --ntasks=8                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=100gb                    # Job memory request
#SBATCH --time=5-00:00:00             # Time limit days-hrs:min:sec
#SBATCH --array=1-5
#SBATCH --output=velocyto_%.log          # Standard output and error log 

pwd; hostname; date

ml R/4.2
ml samtools
ml conda
conda activate velocity

export LC_ALL=en_US.utf-8 && export LANG=en_US.utf-8

Rscript Velocyto_All_Samples.R $SLURM_ARRAY_TASK_ID

date
