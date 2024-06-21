# Generate count matrix with rsem
# Rachel Griffard
# Last updated: 061124

samples = scan('/panfs/pfs.local/home/r816g589/work/collaborations/kalyani/240509/samples.txt', what='', sep='\n')
path = '/panfs/pfs.local/home/r816g589/work/collaborations/kalyani/240509/'

summary.ret = function(path, samples){
  
  # combine rsem result
  all.samples.ret = paste(path, samples, '/',samples, '.genes.results ', sep = '')
  collapse.ret = paste(all.samples.ret, collapse = '')
  combine.gene = paste('rsem-generate-data-matrix ', 
                        collapse.ret,
                        '> ', path, 'genes.count.matrix', sep = '')
  system(combine.gene)
}


summary.ret(path, samples)
