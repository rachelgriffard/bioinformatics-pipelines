# Concatenation process
# setwd()


raw_path = "~/raw_files" # file path to raw files
samples = scan("samples.txt", what="", sep="\n") # text file with sample names on each line
cat_path = "~/cat_files" # file path to new concatenated files

for (sample in samples[1:length(samples)]){
  # concatenation
  system(paste("zcat ", raw_path,
               "/*R1*.fastq.gz >> ", cat_path, samples[sample], "_R1.fastq", sep = ""))
  system(paste("zcat ", raw_path, 
               "/*R2*.fastq.gz >> ", cat_path ,samples[sample], "_R2.fastq", sep = ""))
}

