# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 15:16:47 2024

@author: r816g589
"""
import scanpy as sc
import anndata as ad
import scvelo as scv

# set matplotlib settings to defaults
scv.set_figure_params()

# import data
adata = scv.datasets.pancreas()

scv.pp.filter_and_normalize(adata, min_shared_counts = 20, n_top_genes = 2000)

scv.pp.moments(adata, n_pcs = 30, n_neighbors = 30)

scv.tl.velocity(adata, mode = 'deterministic')

scv.tl.recover_dynamics(adata, n_jobs = 4)

# scv.tl.velocity(adata, mode = 'dynamical')

scv.tl.velocity_graph(adata, n_jobs = -1)
scv.pl.velocity_embedding_stream(adata, basis = 'umap')

top_genes = adata.var['fit_likelihood'].sort_values(ascending = False).index
scv.pl.scatter(adata, basis = top_genes[:10], ncols = 5, frameon = False)

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color = 'latent_time', color_map = 'gnuplot', colorbar = True)
