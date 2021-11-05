import sys
import numpy as np
import scanpy as sc

sc.settings.figdir = ""
input_file = sys.argv[1]
output_file = sys.argv[2]
clustering_method = sys.argv[5]

adata = sc.read(input_file)

if sys.argv[3] == 1:
    if clustering_method != "none":
        sc.pl.umap(adata, color=clustering_method, save=output_file + "_" + clustering_method + '_clusters_UMAP.png')
    else:
        sc.pl.umap(adata, save=output_file + '_UMAP.png')

if sys.argv[4] == 1:
    if clustering_method != "none":
        sc.pl.tsne(adata, color=clustering_method, save=output_file + "_" + clustering_method + '_clusters_tSNE.png')
    else:
        sc.pl.tsne(adata, save=output_file + '_tSNE.png')
