import sys
import scanpy as sc

input_file = sys.argv[1]
output_file = sys.argv[2]
run_umap = int(sys.argv[3])
run_tsne = int(sys.argv[4])

print("dimension reduction input: ", input_file)
adata = sc.read(input_file)

# should it recompute PCA when asked to do umap or tsne?
#if "X_pca" not in adata.obsm.keys() or "PCs" not in adata.varm.keys():
print("-- computing pca --")
sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata)

if run_umap == 1:
    print("-- computing umap --")
    sc.tl.umap(adata)

if run_tsne == 1:
    print("-- computing tsne --")
    sc.tl.tsne(adata)

print("dimension reduction output: ", output_file)
adata.write(output_file, compression='gzip')#, compression_opts=1)
