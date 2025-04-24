#!/bin/bash
#SBATCH --job-name=mkfastq      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=     # Where to send mail	
#SBATCH --ntasks=16                   # Run on a single CPU
#SBATCH --mem=64gb                    # Job memory request
#SBATCH --time=1-0:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=mkfastq_%j.log          # Standard output and error log 

pwd; hostname; date

ml bcl2fastq2
module load R

path/to/cellranger mkfastq \
--id=fastq --run=/path/to/fastq/folder \
--samplesheet=path/to/SampleData.csv
 
date