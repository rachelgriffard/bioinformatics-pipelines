# Quality check using FastQC
# Rachel Griffard
# Last updated: 061124

setwd()

path = "~/cat_files" # file path to concatenated files
samples = scan("samples.txt", what="", sep="\n")
thread.num = 8

system(paste('mkdir fastqc'))
system(paste('fastqc -o fastqc *fastq.gz'))