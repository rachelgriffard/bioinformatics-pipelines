# Upset plots for Enriched v Unenriched - Swerdlow Feb '25
# Rachel Griffard-Smith
#
# Created: 020425
# Last updated: 020425

library(ComplexHeatmap)
library(tidyverse)
library(openxlsx)

fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

overlapGroups <- function (listInput, sort = TRUE) {
  # listInput could look like this:
  # $one
  # [1] "a" "b" "c" "e" "g" "h" "k" "l" "m"
  # $two
  # [1] "a" "b" "d" "e" "j"
  # $three
  # [1] "a" "e" "f" "g" "h" "i" "j" "l" "m"
  listInputmat    <- fromList(listInput) == 1
  #     one   two three
  # a  TRUE  TRUE  TRUE
  # b  TRUE  TRUE FALSE
  #...
  # condensing matrix to unique combinations elements
  listInputunique <- unique(listInputmat)
  grouplist <- list()
  # going through all unique combinations and collect elements for each in a list
  for (i in 1:nrow(listInputunique)) {
    currentRow <- listInputunique[i,]
    myelements <- which(apply(listInputmat,1,function(x) all(x == currentRow)))
    attr(myelements, "groups") <- currentRow
    grouplist[[paste(colnames(listInputunique)[currentRow], collapse = ":")]] <- myelements
    myelements
    # attr(,"groups")
    #   one   two three 
    # FALSE FALSE  TRUE 
    #  f  i 
    # 12 13 
  }
  if (sort) {
    grouplist <- grouplist[order(sapply(grouplist, function(x) length(x)), decreasing = TRUE)]
  }
  attr(grouplist, "elements") <- unique(unlist(listInput))
  return(grouplist)
  # save element list to facilitate access using an index in case rownames are not named
}

# read in data
unenr = read.csv('Proteomics/DAresults_Proteomics.csv')
enr = read.csv('Phosphoproteomics/DAresults_Phosphoproteomics.csv')
head(unenr)
head(enr)

# pull only significant
nrow(unenr[unenr$FDR <= 0.05,])
nrow(enr[enr$FDR <= 0.05,])

unenr_prot = unenr$Accession[unenr$FDR <= 0.05]
enr_prot = enr$Accession[enr$FDR <= 0.05 & !is.na(enr$Accession)]

upset = list(Unenriched = unenr_prot,
             Enriched = enr_prot)

m = make_comb_mat(upset)

png('Overlap_Unenr_Enr.png', width = 600, height = 500, res = 150)
UpSet(m, top_annotation = upset_top_annotation(m, add_numbers = TRUE),
      right_annotation = upset_right_annotation(m, add_numbers = TRUE))#,
      #column_title = 'Overlaping Significant Proteins (FDR<0.05)')
dev.off()


comb = overlapGroups(upset)
c = list()
names = names(comb)
names[names == 'Unenriched:Enriched'] = 'Overlap'
for (i in 1:length(names)) {
  df = data.frame(names(comb[[i]]))
  colnames(df) = names[i]
  c[[i]] = df
  names(c)[[i]] = names[i]
}

wb = createWorkbook()

for (sheet_name in names(c)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, c[[sheet_name]])
}

saveWorkbook(wb, 'overlap_proteins.xlsx', overwrite = TRUE)
