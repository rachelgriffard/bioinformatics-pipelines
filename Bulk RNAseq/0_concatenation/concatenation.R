# Concatenation process
# setwd("S:\\Biostats\\BIO-STAT\\Koestler Devin\\Rachel Griffard\\Collaborations\\Pyaram_Kalyani\\061024_bulkRNA_kalyani")


raw_path = "/panfs/pfs.local/home/r816g589/work/collaborative/Kalyani/240509/raw_files"
samples = scan("samples.txt", what="", sep="\n")
cat_path = "/panfs/pfs.local/home/r816g589/work/collaborative/Kalyani/240509/cat_files"

for (sample in samples[1:length(samples)]){
  # concatenation
  system(paste("zcat ", raw_path,
               "/*R1*.fastq.gz >> ", cat_path, samples[sample], "_R1.fastq", sep = ""))
  system(paste("zcat ", raw_path, 
               "/*R2*.fastq.gz >> ", cat_path ,samples[sample], "_R2.fastq", sep = ""))
}

