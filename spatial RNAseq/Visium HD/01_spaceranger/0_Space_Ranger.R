################################################
# Title: Run Space Ranger
# Author: Emily Schueddig
# Adapted by: Rachel Griffard-Smith
# Date: 02/14/2025
# Last modified: 4/24/25
################################################
library(rio)

samples.fastq = list.files("path/to/samples")

find_matches <- function(input, candidates) {
  # Convert both input and candidates to lowercase for case-insensitive comparison
  input_lower <- tolower(input)
  candidates_lower <- tolower(candidates)
  
  # Find candidates that match the input string or contain it as a substring
  matches <- candidates[grepl(input_lower, candidates_lower)]
  return(matches)
}

args <- commandArgs(trailingOnly=TRUE)
group <- as.numeric(args[1])
n <- as.numeric(args[2])
print(group)

runfolder <- args[3]
imagefolder1 <- args[4]

for (i in list.files(imagefolder1)) {  
  imagefolder = paste0(imagefolder1, i)
  cytafolder <- paste0(imagefolder,"/",list.files(imagefolder)[grep("assay",list.files(imagefolder))])
  # read in sample sheet
  csv = read.csv(paste0(cytafolder,"/",list.files(cytafolder)[grep("csv",list.files(cytafolder))]))
  # get sample names
  samples = csv[grep("Sample",rownames(csv)),1]
  # get slide ID
  slide = csv[grep("Visium Slide ID",rownames(csv)),1]
  
  # split samples
  l = split(samples, cut(seq_along(samples), n, labels=F))
  samples.run = l[[group]]
  
  print(samples.run)
  
  sample.fastq = find_matches(samples.run, samples.fastq)
  
  for (sample in samples.run){
    area=strsplit(sample,"_")[[1]][1]
    cytaimage=list.files(cytafolder)[grep(sample,list.files(cytafolder))]
    
    spaceranger.count <- paste0(args[5],"/spaceranger-3.1.2/spaceranger count --id=", sample, 
                                " --transcriptome=",args[6],"/spaceranger/refdata-gex-mm10-2020-A --fastqs=",runfolder,"/",sample.fastq,
                                " --probe-set=",args[5],"/path/to/prob-set.csv --slide=",slide,
                                " --area=",area," --cytaimage='",cytafolder,"/",cytaimage,"' --image='",imagefolder,"/",sample,"_",slide,".tif' --create-bam=false --custom-bin-size=4")
    print(spaceranger.count)
    system(spaceranger.count)
  }
}
