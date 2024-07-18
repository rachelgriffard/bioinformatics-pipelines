# Quality check using FastQC
# Rachel Griffard
# Last updated: 061124

setwd()

path = "~/cat_files" # file path to concatenated files
samples = scan("samples.txt", what="", sep="\n")
thread.num = 8

system(paste0('mkdir fastqc'))
system(paste0('fastqc -o fastqc ', path, '/*.fastq.gz'))

system(paste0('cd fastqc'))
system(paste0('multiqc .'))