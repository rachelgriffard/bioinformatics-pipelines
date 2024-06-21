# Quantification using RSEM
# Rachel Griffard
# Last updated: 061124

setwd()

path = "~/cat_files" # file path to concatenated files
samples = scan("samples.txt", what="", sep="\n")
thread.num = 8

for (sample in samples[1:6]){
  # mapper using RSEM -p for thread
  system(paste("mkdir ", sample, sep = ""))
  rsem = paste('rsem-calculate-expression -p ', thread.num, ' --bowtie2 --paired-end --output-genome-bam ', 
                path, sample, "_R1.fastq ",  
                path, sample, "_R2.fastq ~/hg38-rsem ", # file path to reference genome
                sample, "/", sample, sep = "")
  
  system(rsem)
}
