import sys
import scanpy as sc

print("gene filtering input:", sys.argv[1])
print("gene filtering output:", sys.argv[2])

adata = sc.read(sys.argv[1])

min_counts = int(sys.argv[3])
max_counts = int(sys.argv[4])
min_cells = int(sys.argv[5])
max_cells = int(sys.argv[6])

print("number of genes: ", adata.shape[1])
if min_counts > 0:
    print("filtering out genes with less than", min_counts, "counts")
    sc.pp.filter_genes(adata, min_counts=min_counts)
    print("genes remaining: ", adata.shape[1])
if max_counts > 0:
    print("filtering out genes with more than", max_counts, "counts")
    sc.pp.filter_genes(adata, max_counts=max_counts)
    print("genes remaining: ", adata.shape[1])
if min_cells > 0:
    print("filtering out genes with less than", min_cells, "cells")
    sc.pp.filter_genes(adata, min_cells=min_cells)
    print("genes remaining: ", adata.shape[1])
if max_cells > 0:
    print("filtering out genes with more than", max_cells, "cells")
    sc.pp.filter_genes(adata, max_cells=max_cells)
    print("genes remaining: ", adata.shape[1])

adata.write(sys.argv[2], compression='gzip')#, compression_opts=1)
