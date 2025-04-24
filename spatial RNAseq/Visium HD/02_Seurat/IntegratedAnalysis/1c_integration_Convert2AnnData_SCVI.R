# convert Seurat to h5ad


merged.obj = readRDS('RObjects/02_mergedSeurat.RDS')

merged.obj[["Spatial.008um"]] = JoinLayers(merged.obj[["Spatial.008um"]])

# Subset to only the counts layer
merged.obj.sub = DietSeurat(merged.obj, assays = "Spatial.008um", layers = "counts")

# Set the data layer to be the counts
merged.obj.sub[["Spatial.008um"]]$data = merged.obj.sub[["Spatial.008um"]]$counts

# remove the FOVs just because
merged.obj.sub[["slice1.008um"]] = NULL
merged.obj.sub[["slice1.008um.2"]] = NULL
merged.obj.sub[["slice1.008um.3"]] = NULL
merged.obj.sub[["slice1.008um.4"]] = NULL

# removing all the meta features information which is not needed
merged.obj.sub@assays$Spatial.008um@meta.data = data.frame(c())

# convert back to Class = Assay
merged.obj.sub[["Spatial.008um"]] = as(object = merged.obj.sub[["Spatial.008um"]], Class = "Assay")

# set the number of variable features
merged.obj.sub = FindVariableFeatures(merged.obj.sub, nfeatures = 2000)

# save as h5Seurat
SaveH5Seurat(merged.obj.sub, filename = "RObjects/04b_ScanPyObject.h5Seurat", overwrite=TRUE, verbose=TRUE)

# convert to h5ad
Convert("RObjects/04b_ScanPyObject.h5Seurat", dest = "h5ad", overwrite=TRUE)

# Now you can run the Integrate_scVI.py script which uses the h5ad object you just created