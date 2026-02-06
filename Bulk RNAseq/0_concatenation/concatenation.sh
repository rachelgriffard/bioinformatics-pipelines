#!/bin/bash
#SBATCH --job-name=concatenation      # Job name
#SBATCH --partition=sixhour           # Partition Name (Required)
#SBATCH --mail-type=ALL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=     # Where to send mail	
#SBATCH --ntasks=16                   # Run on a single CPU
#SBATCH --mem=64gb                    # Job memory request
#SBATCH --time=2:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=concatenation_%j.log          # Standard output and error log 

pwd; hostname; date


module load R

Rscript concatenation.R

 
date
 
date