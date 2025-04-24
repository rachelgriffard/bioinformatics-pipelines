# Variant table creation
# Dec 2024
# Created by: Rachel Griffard-Smith

library(biomaRt)
library(tidyverse)

mart = useEnsembl(biomart = 'snp', dataset = 'sscrofa_snp') # select right dataset

pig_vep_pkd1 = read.csv('reference/PKD1-pig.csv') # import desired reference from ensembl
pig_vep_pkd2 = read.csv('reference/PKD2-pig.csv') # import desired reference from ensembl
pig_var = rbind(pig_vep_pkd1, pig_vep_pkd2)
colnames(pig_var) = c('Variant ID',
                      'vf',
                      'Location',
                      'Chr: bp',
                      'vf_allele',
                      'Alleles',
                      'Class',
                      'Source',
                      'Evidence',
                      'Clin. Sig.',
                      'Conseq. Type',
                      'AA',
                      'AA coord',
                      'AA coord2',
                      'sift_sort',
                      'sift_class',
                      'SIFT',
                      'Transcript'
)

pig_var = pig_var[,c(1,3,6,16,17,18)]

write.csv(pig_var, 'reference/pig_var.csv')

pig_var = read.csv('reference/pig_var.csv', row.names = 1)

# pull in files with variants used
file_list = list.files(path = "tables/", pattern = "*.csv")

file_list = file_list[!file_list %in% c('variants_table.csv')]

setwd('./tables')

data_list = lapply(file_list, read.csv)

variants = unique(unlist(purrr::map(data_list, 1)))

variants_df = data.frame(do.call(rbind, strsplit(variants, ":")))

variants_call = variants_df %>%
  mutate('chromosome' = substr(X1, 4, nchar(X1)),
         'start' = X2,
         'reference_allele' = sapply(strsplit(X3, '/'), '[', 1),
         'alternate_allele' = sapply(strsplit(X3, '/'), '[', 2)) %>%
  select(-c('X1', 'X2', 'X3'))

variant_info = getBM(
  attributes = c('refsnp_id'),
  filters = c('chr_name', 'start', 'end'), # 'allele'),
  values = list(variants_call$chromosome, variants_call$start, variants_call$start),
                # paste0(variants_call$reference_allele, '/', variants_call$alternate_allele)),
  mart = mart
)

library(httr)
library(jsonlite)

unlist_single_element_list <- function(x) {
  if (length(x) == 1) {
    unlist(x)
  } else {
    x
  }
}

table(nchar(variants_df$X3))/nrow(variants_df)

variants_df %>%
  mutate(sapply(X3, nchar()))

r = fromJSON(toJSON(content(r)))

variant_ids = pig_var$Variant.ID

server = "https://rest.ensembl.org"
ext = "/vep/sus_scrofa/id/"

variant_list = list()

r = POST(
  paste0(server, ext),
  query = list('CADD' = 1),
  content_type_json(),
  accept_json(),
  body = toJSON(list('ids' = variant_ids))
)

r = fromJSON(toJSON(content(r)))

variant_df = r %>%
  as.data.frame() %>%
  unnest_wider(col = colocated_variants,
               names_sep = '-') %>%
  mutate(id = input) %>%
  unnest_longer(col = transcript_consequences) %>%
  unnest_wider(col = transcript_consequences,
               names_sep = "_") %>%
  group_by(input) %>%
  relocate(input, .before = transcript_consequences_gene_symbol_source)



variant_df = r %>%
  as.data.frame() %>%
  unnest_wider(col = colocated_variants,
               names_sep = '-') %>%
  mutate(id = input) %>%
  unnest_longer(col = transcript_consequences) %>%
  unnest_wider(col = transcript_consequences,
               names_sep = "_") %>%
  group_by(input) %>%
  relocate(input, .before = transcript_consequences_gene_symbol_source) %>%
  rename(variant_id = input,
         CADD_PHRED = transcript_consequences_cadd_phred,
         chromosome = seq_region_name,
         allele = allele_string) %>%
  select(variant_id, chromosome, start, end, allele, CADD_PHRED, transcript_consequences_gene_symbol) %>%
  mutate(CADD_PHRED = map(CADD_PHRED, unlist_single_element_list),
         transcript_consequences_gene_symbol = map(transcript_consequences_gene_symbol, unlist_single_element_list)) %>%
  unnest(c(variant_id, chromosome, start, end, allele, CADD_PHRED, transcript_consequences_gene_symbol)) %>%
  group_by(variant_id) %>%
  summarize(
    chromosome = unique(chromosome),
    start = unique(start),
    end = unique(end),
    allele = unique(allele),
    CADD_PHRED = unique(CADD_PHRED),
    transcript_consequences_gene_symbol = paste(unique(transcript_consequences_gene_symbol), collapse = ',')
  )
  
  
  # filter(transcript_consequences_gene_symbol %in% c('PKD1', 'PKD2'))

write.csv(variant_df, 'tables/variants_table.csv', row.names = FALSE)


# check all variants present in table
variants_table = read.csv('tables/variants_table.csv')


