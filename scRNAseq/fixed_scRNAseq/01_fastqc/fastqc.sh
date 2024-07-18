#!/bin/bash
#SBATCH --job-name=mkfastq      # Job name
#SBATCH --partition=biostat           # Partition Name (Required)
#SBATCH --mail-type=BEGIN,END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rgriffard@kumc.edu     # Where to send mail	
#SBATCH --ntasks=16                   # Run on a single CPU
#SBATCH --mem=64gb                    # Job memory request
#SBATCH --time=1-0:00:00             # Time limit days-hrs:min:sec
#SBATCH --output=mkfastq_%j.log          # Standard output and error log 

pwd; hostname; date

ml bcl2fastq2
module load R

/kuhpc/work/biostat/r816g589/tools/cellranger-8.0.1/cellranger mkfastq \
--id=fastq --run=/kuhpc/work/biostat/r816g589/collaborations/sundar/20240712-sundar-scRNAseq/240531_A00484_0483_AH7GC3DRX5 \
--samplesheet=/kuhpc/work/biostat/r816g589/collaborations/sundar/20240712-sundar-scRNAseq/240531_A00484_0483_AH7GC3DRX5/SampleData.csv
 
date