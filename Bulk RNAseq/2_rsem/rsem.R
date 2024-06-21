# Quantification using RSEM
# Rachel Griffard
# Last updated: 061124

setwd("S:\\Biostats\\BIO-STAT\\Koestler Devin\\Rachel Griffard\\Collaborations\\Pyaram_Kalyani\\061024_bulkRNA_kalyani\\mouse")

path = "/panfs/pfs.local/home/r816g589/work/collaborative/Kalyani/240509/cat_files"
samples = scan("samples.txt", what="", sep="\n")
thread.num = 8

for (sample in samples[1:6]){
  # mapper using RSEM -p for thread
  system(paste("mkdir ", sample, sep = ""))
  rsem = paste('rsem-calculate-expression -p ', thread.num, ' --bowtie2 --paired-end --output-genome-bam ', 
                path, sample, "_R1.fastq ",  
                path, sample, "_R2.fastq /panfs/pfs.local/work/biostat/d324p169/reference/fridleyRef/ucsc-hg38-rsem/hg38-rsem ", 
                sample, "/", sample, sep = "")
  
  system(rsem)
}
