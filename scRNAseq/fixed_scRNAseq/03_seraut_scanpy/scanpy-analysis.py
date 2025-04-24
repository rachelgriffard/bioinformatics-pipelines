# -*- coding: utf-8 -*-
"""
scanpy for scRNAseq analysis

@author: Rachel Griffard-Smith

last update: 042425

ref: https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html
"""

import scanpy as sc
import anndata as ad
# import plotnine as pn
import scanpy.external as sce
import pandas as pd
# import bbknn


# set plot settings
sc.settings.set_figure_params(dpi=50, facecolor="white")

### Read in data
# read in the count matrix to an AnnData object
# samples = open('../samples.txt').read().splitlines()

samples = {
    "S1":"..path/to/sample/S1/count/sample_filtered_feature_bc_matrix.h5",
    "S2":"..path/to/sample/S2/count/sample_filtered_feature_bc_matrix.h5",
    "S3":"..path/to/sample/S3/count/sample_filtered_feature_bc_matrix.h5",
    "S4":"..path/to/sample/S4/count/sample_filtered_feature_bc_matrix.h5",
}

adatas = {}

for sample,filename in samples.items():
    sample_adata = sc.read_10x_h5(filename)
    sample_adata.var_names_make_unique()
    adatas[sample] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
print(adata.obs["sample"].value_counts())
adata

### Quality control (**done separately for each sample)
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(P)]")

# calculate qc
sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

# violin plot to compare computed metrics from above
sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts")

# filter low quality
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

### Doublet detection - although perhaps not necessary, set upper range on filtering quality
sc.pp.scrublet(adata, batch_key="sample")

### Normalization
# Saving count data
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)


### Feature selection
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata)

var_select = adata.var.highly_variable_nbatches > 1
var_genes = var_select.index[var_select]

### Dimension reduction
# pca
sc.tl.pca(adata)

# var ratio
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# plot pcas against features like batch, qc metrics, etc
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
)

### Nearest neighbors
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

### Clustering
# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])

### QC/cell filtering using UMAP
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)

sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

### Integration
alldata = {}
for sample in samples:
    alldata[sample] = adata[adata.obs['sample'] == sample,]
alldata

sce.pp.harmony_integrate(adata, 'sample')
'X_pca_harmony' in adata.obsm

sc.pp.neighbors(adata, use_rep='X_pca_harmony')
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2
)

sc.pl.umap(
    adata,
    color="pct_counts_mt",
    # Setting a smaller point size to get prevent overlap
    size=2
)

sc.pl.umap(
    adata,
    color="leiden",
    # Setting a smaller point size to get prevent overlap
    size=2
)

"""
## BBKNN integration
sce.pp.bbknn(adata, batch_key="sample")  # running bbknn 1.3.6
sc.tl.umap(adata)
sc.pl.umap(adata, color=["sample", "louvain"])

# SCVI
# tell scvi where to get the data
scvi.data.setup_anndata(adata, layer="counts", batch_key = 'sample')
scvi.data.view_anndata_setup(adata)
"""
"""
import csv

with open('../markergenes/azimuth_pred.csv', newline='') as f:
    reader = csv.reader(f)
    pred = list(reader)
"""
pred = pd.read_csv('../markergenes/azimuth_pred.csv')

adata.obs["barcode"] = adata.obs['sample'].astype(str) + '_' + adata.obs_names

ad_ob = pd.DataFrame(adata.obs)

adata.obs = pd.merge(adata.obs, pred, on = 'barcode', how = 'left')

sc.pl.umap(
    adata,
    color="predicted.celltype_level3",
    # Setting a smaller point size to get prevent overlap
    size=2
)