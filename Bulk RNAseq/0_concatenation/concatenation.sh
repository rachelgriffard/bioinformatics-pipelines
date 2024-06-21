#!/bin/bash
#SBATCH --job-name=concatenation      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rgriffard@kumc.edu     # Where to send mail	
#SBATCH --ntasks=16                   # Run on a single CPU
#SBATCH --mem=64gb                    # Job memory request
#SBATCH --time=1-0:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=concatenation_%j.log          # Standard output and error log 

pwd; hostname; date


module load R

Rscript concatenation.R

 
date