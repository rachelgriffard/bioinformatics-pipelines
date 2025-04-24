#!/bin/bash
#SBATCH --job-name=spaceRanger      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=ALL       # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=     # Where to send mail	
#SBATCH --ntasks=1                   # Run on a single CPU
#SBATCH --cpus-per-task=12
#SBATCH --nodes=1
#SBATCH --mem=64gb                    # Job memory request -- set to appropriate parameters
#SBATCH --time=03:00:00             # Time limit days-hrs:min:sec -- set to appropriate parameters
#SBATCH --array=1-6               # Set to number of arrays you want to split
#SBATCH --output=spaceranger_%A_%a.log          # Standard output and error log 

pwd; hostname; date

#source ~/.bashrc

ml R/4.4

IMAGEFOLDER=/path/to/image/folder

RUNFOLDER=/path/to/fastq/folder

TOOLS=/path/to/spaceranger
REFERENCES=/path/to/reference/sequences

# Execute the R script (uncomment when ready)
Rscript 0_Space_Ranger.R $SLURM_ARRAY_TASK_ID $SLURM_ARRAY_TASK_MAX $RUNFOLDER $IMAGEFOLDER $TOOLS $REFERENCES

date