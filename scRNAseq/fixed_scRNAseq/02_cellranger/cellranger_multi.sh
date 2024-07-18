#!/bin/bash
#SBATCH --job-name=10x_multi_1                # Job name
#SBATCH --partition=biostat                # Partition Name (Required)
#SBATCH --mail-type=ALL                    # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=rgriffard@kumc.edu     # Where to send mail	
#SBATCH --nodes=1
#SBATCH --ntasks=16                      # Run a single task	
#SBATCH --mem=120gb                  # Job memory request
#SBATCH --time=3-00:00:00                  # Time limit days-hrs:min:sec
#SBATCH --output=cellranger_%j.log           # Standard output and error log

/kuhpc/work/biostat/r816g589/tools/cellranger-8.0.1/cellranger multi --id=multi_output --csv config.csv

echo Finished execution at `date`