#!/bin/bash
#SBATCH --job-name=visCellType      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=    # Where to send mail	
#SBATCH --ntasks=1                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=125gb                    # Job memory request
#SBATCH --time=2-0:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=cell_type_%a.log          # Standard output and error log 

pwd; hostname; date

source ~/.bashrc

ml R
ml hdf5

Rscript 3_cell-id.R

date
