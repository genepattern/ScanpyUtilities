import sys
import scanpy as sc

print("cell filtering input:", sys.argv[1])
print("cell filtering output:", sys.argv[2])

adata = sc.read(sys.argv[1])

min_counts = int(sys.argv[3])
max_counts = int(sys.argv[4])
min_genes = int(sys.argv[5])
max_genes = int(sys.argv[6])

print("number of cells: ", adata.shape[0])
if min_counts > 0:
    print("filtering out cells with less than", min_counts, "counts")
    sc.pp.filter_cells(adata, min_counts=min_counts)
    print("cells remaining: ", adata.shape[0])
if max_counts > 0:
    print("filtering out cells with more than", max_counts, "counts")
    sc.pp.filter_cells(adata, max_counts=max_counts)
    print("cells remaining: ", adata.shape[0])
if min_genes > 0:
    print("filtering out cells with less than", min_genes, "genes")
    sc.pp.filter_cells(adata, min_genes=min_genes)
    print("cells remaining: ", adata.shape[0])
if max_genes > 0:
    print("filtering out cells with more than", max_genes, "genes")
    sc.pp.filter_cells(adata, max_genes=max_genes)
    print("cells remaining: ", adata.shape[0])

adata.write(sys.argv[2], compression='gzip', compression_opts=1)
