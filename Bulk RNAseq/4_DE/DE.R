# Differential Expression with DESeq2
# Rachel Griffard
# Last updated: 061324

library(tidyverse)

samples = scan('samples.txt', what='', sep='\n')

dat = read.delim('genes.count.matrix')

dat = dat %>% 
  column_to_rownames('X') # rownames as genes

colnames(dat) = samples # rename columns as samples (check)

# DE with DESeq2
library(DESeq2)
library(edgeR) # DGEList
library(DEFormats) # DGEList to DGESeqDataSet

group = c('WT', 'KO', 'KO', 'KO', 'WT', 'WT') # grouping ordered by columns

coldata = data.frame(row.names = samples,
                     condition = group)

all(rownames(coldata) %in% colnames(dat))

dds = DESeqDataSetFromMatrix(countData = round(dat), # round bc DESeq req integer
                             colData = coldata,
                             design = ~ condition)

sgs = 3 # smallest group size
keep = rowSums(counts(dds) >= 10) >= sgs
dds = dds[keep,]

dds$condition = relevel(dds$condition, ref = 'WT') # set wildtype to reference

dds = DESeq(dds)
res_list = results(dds)

res = data.frame(row.names = res_list@rownames, res_list@listData)

write.csv(res, '')