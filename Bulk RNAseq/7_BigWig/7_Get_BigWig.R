###########################################
# Title: Get bigWig files
# Author: Emily Nissen
# Date Created: 4/8/2024
# Last Modified: 
##########################################

args <- commandArgs(trailingOnly=TRUE)
group <- as.numeric(args[1])

path = "/kuhpc/work/biostat/r816g589/collaborations/kalyani/250311"
data.path = "/kuhpc/work/biostat/r816g589/collaborations/kalyani/250311"
out.path = "/kuhpc/work/biostat/r816g589/collaborations/kalyani/250311/bigWigs"

thread.num = 8

samples = c('D2_WT3',
                'D2_WT2', 
                'D2_WT1',
                'D2_KKO3',
                'D2_KKO2',
                'D2_KKO1',
                'D1_WT3',
                'D1_WT2',
                'D1_WT1',
                'D1_KKO3',
                'D1_KKO2',
                'D1_KKO1')

# if(group == 1){
#   samples = samples.all[1:5]
# }else if(group == 2){
#   samples = samples.all[6:10]
# }else if(group == 3){
#   samples = samples.all[11:15]
# }else if(group == 4){
#   samples = samples.all[16:20]
# }

for(sample in samples){
  sort = paste0("samtools sort ", data.path, "/", sample, "/", sample, ".genome.bam -o ",
                data.path, "/", sample, "/", sample, ".genome.sorted.bam")
  print(sort)
  system(sort)

  idx = paste0("samtools index ", data.path, "/", sample, "/", sample, ".genome.sorted.bam")
  print(idx)
  system(idx)
  
  bw = paste0("bamCoverage --bam ", data.path, "/", sample, "/", sample, ".genome.sorted.bam -o ",
              out.path, "/", sample, "_v2.bw")
  print(bw)
  system(bw)
}