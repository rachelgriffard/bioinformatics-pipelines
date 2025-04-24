################################################
# Title: velocyto on scRNAseq
# Author: Rachel Griffard (adapted from Emily Nissen)
# Last updated: 070624
################################################

args <- commandArgs(trailingOnly=TRUE)
group <- as.numeric(args[1])

sample.ids <- c("HDM673_WT", "HDM671_Mut")

# if(group == 1){
#   sample.ids = sample.ids[1:5]
# }else if(group == 2){
#   sample.ids = sample.ids[6:10]
# }else if(group == 3){
#   sample.ids = sample.ids[11:15]
# }else if(group == 4){
#   sample.ids = sample.ids[16:20]
# }else if(group == 5){
#   sample.ids = sample.ids[21:25]
# }


for (sample.id in sample.ids){
  velo = paste0("velocyto run10x --samtools-threads 8 /path/to/sample/", sample.id, 
                "/count path/to/reference/genes/genes.gtf")

  print(velo)
  system(velo)
}

