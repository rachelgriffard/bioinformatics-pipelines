# Quality check using FastQC
# Rachel Griffard
# Last updated: 061124

setwd("S:\\Biostats\\BIO-STAT\\Koestler Devin\\Rachel Griffard\\Collaborations\\Pyaram_Kalyani\\061024_bulkRNA_kalyani\\mouse")

path = "/panfs/pfs.local/home/r816g589/work/collaborative/Kalyani/240509/cat_files"
samples = scan("samples.txt", what="", sep="\n")
thread.num = 8

system(paste('mkdir fastqc'))
system(paste('fastqc -o fastqc *fastq.gz'))