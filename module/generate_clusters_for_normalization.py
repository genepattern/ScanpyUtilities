import sys
import numpy as np
import scanpy as sc

print("normalization input:", sys.argv[1])
print("clustering cells to use as input for scran")

adata = sc.read(sys.argv[1])

sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e6)
sc.pp.log1p(adata)
print("running pca")
sc.pp.pca(adata, n_comps=15)
print("computing neighbors")
sc.pp.neighbors(adata)
print("running louvain clustering")
sc.tl.louvain(adata, key_added='groups', resolution=0.5)

adata.write("temp_clustered_for_scran.h5ad", compression='gzip')#, compression_opts=1)
