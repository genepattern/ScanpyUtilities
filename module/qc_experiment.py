import os, sys, re, argparse, shutil
from optparse import OptionParser

import anndata
import anndata as ad
import scanpy as sc
import seaborn as sns

__author__ = "Anthony S. Castanza"
__email__ = "acastanza@ucsd.edu"
__version__ = "1.0.0"


def main():

    sc.settings.figdir = ""

    adata = sc.read(sys.argv[1])
    plots = ['n_genes_by_counts', 'total_counts']

    if sys.argv[2] != "SKIP":
        with open(sys.argv[2]) as f:
            mito_genes = f.read().splitlines()
        mito_genes = list(set([re.sub('-I$', '', sub) for sub in mito_genes]))
        adata.var['mt'] = [x in mito_genes for x in adata.var_names]
    else:
        if 'gene_name' in adata.var:
            gene_name_uppercase = [x.upper() for x in list(adata.var['gene_name'])]
        else:
            gene_name_uppercase = [x.upper() for x in list(adata.var.index)]
        is_mt = list(map(lambda x:x.startswith('MT-'),gene_name_uppercase))
        mito_genes = list(set(adata.var.index[is_mt]))
        if len(mito_genes) > 0:
            adata.var['mt'] = [x in mito_genes for x in adata.var_names]

    if sys.argv[2] != "SKIP" or 'mt' in adata.var:
        sc.pp.calculate_qc_metrics(
            adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        plots = plots + ['pct_counts_mt']
    else:
        sc.pp.calculate_qc_metrics(
            adata, percent_top=None, log1p=False, inplace=True)

    sc.pl.violin(adata, plots,
                 jitter=0.4, multi_panel=True, save ="_" + sys.argv[3] + "_qc_violin_plots.png")

    if sys.argv[2] != "SKIP" or 'mt' in adata.var:
        sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',
                      save="_" + sys.argv[3] + '_qc_pct_counts_mitochondrial_vs_total_counts.png')

    sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
                  save="_" + sys.argv[3] + '_qc_n_genes_by_counts_vs_total_counts.png')

    sc.pl.scatter(adata, x='total_counts', y='n_cells_by_counts',
                  save="_" + sys.argv[3] + '_qc_n_cells_by_counts_vs_total_counts.png')

    sns.jointplot(
        data=adata.obs,
        x="total_counts",
        y="n_genes_by_counts",
        kind="hex",
    ).savefig(sys.argv[3] + '_qc_total_genes_by_counts_vs_total_counts_with_histograms.png')

if __name__ == '__main__':
    main()
