import sys
import numpy as np
import scanpy as sc

print("annotation input:", sys.argv[1])
adata = sc.read(sys.argv[1])

if 'n_counts' not in adata.obs.keys():
    adata.obs['n_counts'] = adata.X.sum(1)
if 'log_counts' not in adata.obs.keys():
    adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
if 'n_genes' not in adata.obs.keys():
    adata.obs['n_genes'] = (adata.X > 0).sum(1)

if 'n_counts' not in adata.var.keys():
    rsum = adata.X.sum(0)
    adata.var['n_counts'] = np.array(rsum)[0].tolist()
if 'n_cells' not in adata.var.keys():
    rsum = (adata.X > 0).sum(0)
    adata.var['n_cells'] = np.array(rsum)[0].tolist()

print("annotation output:", sys.argv[2])
adata.write(sys.argv[2], compression='gzip')#, compression_opts=1)
