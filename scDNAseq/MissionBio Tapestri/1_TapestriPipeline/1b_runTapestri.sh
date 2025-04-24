#!/bin/bash
#SBATCH --job-name=tapestriPipeline      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=START,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=     # Where to send mail	
#SBATCH --ntasks=1                   # Run on a single CPU
#SBATCH --nodes=1
#SBATCH --mem=120gb                    # Job memory request
#SBATCH --cpus-per-task=32
#SBATCH --time=2-0:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=tapestri_%j.log          # Standard output and error log 

pwd; hostname; date

eval "$(path/to/TapestriPipeline/bin/conda shell.bash hook)" # load tapestri pipeline

for sample in $(cat samples.txt | tr -d '\r'); do # make sure dos2unix for samples.txt
  cmd="tapestri dna run \
  --output-folder $sample \
  --config ${sample}_config.yaml \
  --genome path/to/fasta/file.fa \
  --genome-version genome_name \
  --panel path/to/panel/ \ # match yaml
  --r1 fastq/${sample}_L001_R1_001.fastq.gz \
  --r2 fastq/${sample}_L001_R2_001.fastq.gz \
  --output-prefix $sample --n-cores 32"
  echo $cmd

  
  # Run the command
  eval $cmd
done
 
date