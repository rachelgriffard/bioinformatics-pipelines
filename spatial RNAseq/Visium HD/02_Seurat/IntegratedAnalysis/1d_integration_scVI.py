import scvi
import scanpy as sc
import time
import pandas as pd
start_time = time.time()
# Load the date
adata = sc.read_h5ad("/path/to/integration/RObjects/04b_ScanPyObject.h5ad")


# Setup and train scVI
scvi.model.SCVI.setup_anndata(adata, batch_key="orig.ident")
model = scvi.model.SCVI(adata, n_latent=10, n_layers=1, gene_likelihood="zinb") # default params
model.train()

# Store the latent representation
adata.obsm["X_scVI"] = model.get_latent_representation()

# Save to .h5ad
adata.write("/path/to/integration/RObjects/04c_IntegratedData_scVI_default.h5ad")

print(time.time() - start_time)
print("Starting clustering...")
# compute neighbor
sc.pp.neighbors(adata, use_rep="X_scVI")
# cluster
sc.tl.leiden(adata)
print(time.time() - start_time)
print("Starting UMAP...")
# umap
sc.tl.umap(adata, min_dist=0.3)
print(time.time() - start_time)

# Save to .h5ad
adata.write("/path/to/integration/RObjects/04c_IntegratedData_scVI_LeidenCluster_UMAP_default.h5ad")

print("Staring UMAP with cell labels...")
sc.pp.neighbors(adata, use_rep="cell.pred")

df = adata.obs
df.to_csv('metadata_scvi.csv')
