#!/bin/bash
#SBATCH --job-name=tapestriPipeline      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=START,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rgriffard@kumc.edu     # Where to send mail	
#SBATCH --ntasks=1                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=120gb                    # Job memory request
#SBATCH --cpus-per-task=24
#SBATCH --time=2-0:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=tapestri_%j.log          # Standard output and error log 

pwd; hostname; date

eval "$($WORK/conda/envs/TapestriPipeline/bin/conda shell.bash hook)"

for sample in $(cat samples.txt); do
  cmd="tapestri dna run \
  --output-folder $sample \
  --genome /kuhpc/work/biostat/r816g589/reference/sscrofa11.1/susscrofa11.fa \
  --genome-version Sscrofa11.1 \
  --panel panel \
  --r1 concat/${sample}_L001_R1_001.fastq.gz \
  --r2 concat/${sample}_L001_R2_001.fastq.gz \
  --output-prefix $sample --n-cores 32"

  # Print the command
  echo $cmd
  
  # Run the command
  eval $cmd
done
 
date