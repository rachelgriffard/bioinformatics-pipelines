#!/bin/bash
#SBATCH --job-name=parallel_gzip    # Job name
#SBATCH --output=parallel_gzip_%A_%a.out  # Output file for each task
#SBATCH --error=parallel_gzip_%A_%a.err   # Error file for each task
#SBATCH --array=0-11                # Array range (0 to 11 for 12 tasks)
#SBATCH --time=4:00:00             # Time limit (HH:MM:SS)
#SBATCH --mem=100gb                    # Memory per task
#SBATCH --partition=biostat        # Partition name


# Create an array of FASTQ files
FILES=(*.fastq)

# Get the file corresponding to this SLURM task ID
FILE=${FILES[$SLURM_ARRAY_TASK_ID]}

# Compress the file
gzip "$FILE"