import sys
import scanpy as sc
from scipy.sparse import csr_matrix

print("normalizing using size factors and log(D + 1)")

adata = sc.read(sys.argv[1])

size_factors = sc.read_csv("temp_size_factors.csv", first_column_names=True)
adata.obs['size_factors'] = size_factors.X

adata.layers['counts'] = adata.X.copy()
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)

adata.X = csr_matrix(adata.X)
adata.write(sys.argv[2], compression='gzip')#, compression_opts=1)
