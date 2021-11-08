import sys
import numpy as np
import scanpy as sc
from matplotlib import pyplot as plt

sc.settings.figdir = "./"
input_file = sys.argv[1]
output_file = sys.argv[2]
clustering_method = sys.argv[5]

adata = sc.read(input_file)

if int(sys.argv[3]) == 1:
    if clustering_method != "none":
        print("Plotting UMAP with" + clustering_method + "clusters.")
        with plt.rc_context():
            sc.pl.umap(adata, color=clustering_method, show=False)
            plt.savefig(output_file + "_" + clustering_method + '_clusters_UMAP.png', bbox_inches="tight")
    else:
        print("Plotting UMAP, no clustering specified.")
        with plt.rc_context():
            sc.pl.umap(adata, show=False)
            plt.savefig(output_file + '_UMAP.png', bbox_inches="tight")

if int(sys.argv[4]) == 1:
    if clustering_method != "none":
        print("Plotting tSNE with" + clustering_method + "clusters.")
        sc.pl.tsne(adata, color=clustering_method, show=False)
        plt.savefig(output_file + "_" + clustering_method + '_clusters_tSNE.png', bbox_inches="tight")
    else:
        print("Plotting tSNE, no clustering specified.")
        sc.pl.tsne(adata, show=False)
        plt.savefig(output_file + '_tSNE.png', bbox_inches="tight")
