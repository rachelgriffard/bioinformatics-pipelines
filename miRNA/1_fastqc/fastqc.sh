#!/bin/bash
#SBATCH --job-name=FastQC      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rgriffard@kumc.edu     # Where to send mail	
#SBATCH --ntasks=8                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=16gb                    # Job memory request
#SBATCH --time=140:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=RSEM_%j.log          # Standard output and error log 

pwd; hostname; date

module load bowtie2
module load java 
module load R
module load fastqc

 
echo "Running R script"
 
Rscript fastqc.R
date
