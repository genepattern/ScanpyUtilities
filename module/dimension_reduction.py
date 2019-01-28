import sys
import scanpy.api as sc

input_file = sys.argv[1]
output_file = sys.argv[2]
run_umap = int(sys.argv[3])
run_tsne = int(sys.argv[4])

adata = sc.read(input_file)

sc.pp.pca(adata, n_comps=50, svd_solver='arpack')
sc.pp.neighbors(adata)

if run_umap == 1:
    print("computing umap")
    sc.tl.umap(adata)

if run_tsne == 1:
    print("computing tsne")
    sc.tl.tsne(adata)

adata.write(output_file, compression='gzip', compression_opts=1)