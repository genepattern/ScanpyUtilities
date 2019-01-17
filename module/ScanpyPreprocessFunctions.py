import h5py
import scanpy.api as sc
import numpy as np
import sys
import os
from scipy.sparse import csr_matrix

def addDataToGroup(group, name, data):
    group.create_dataset(name, data=data, compression='gzip',
        compression_opts=1, shuffle=True)

def read10xH5(in_file, genome):
    adata = sc.read_10x_h5(in_file, genome=genome)
    hf = h5py.File(in_file, 'r')
    all_keys = list(hf[genome].keys())
    keys_10x = ['barcodes', 'data', 'gene_names', 'genes', 'indices', 'indptr', 'shape']
    extra_keys = [x for x in all_keys if x not in keys_10x]
    for key in extra_keys:
        adata.obs[key] = hf[genome][key].value
    return adata

def write10xH5(adata, out_file, genome):
    if os.path.isfile(out_file):
        print("Output file already exists!")
        sys.exit()
    hf = h5py.File(out_file, 'w')
    g = hf.create_group(genome)
    sparse_mat = csr_matrix(adata.X)
    addDataToGroup(g, 'barcodes', adata.obs.axes[0].values.astype('|S'))
    addDataToGroup(g, 'data', sparse_mat.data.astype('i4'))
    addDataToGroup(g, 'gene_names', adata.var.axes[0].values.astype('|S'))
    addDataToGroup(g, 'genes', adata.var['gene_ids'].values.astype('|S'))
    addDataToGroup(g, 'indices', sparse_mat.indices.astype('i8'))
    addDataToGroup(g, 'indptr', sparse_mat.indptr.astype('i8'))
    addDataToGroup(g, 'shape', np.array([sparse_mat.shape[1], sparse_mat.shape[0]], dtype='i4'))
    for key in adata.obs.keys():
        addDataToGroup(g, key, adata.obs[key])
    hf.close()

def annotateData(adata):
    if 'n_counts' not in adata.obs.keys():
        adata.obs['n_counts'] = adata.X.sum(1)
    if 'log_counts' not in adata.obs.keys():
        adata.obs['log_counts'] = np.log(adata.obs['n_counts'])
    if 'n_genes' not in adata.obs.keys():
        adata.obs['n_genes'] = (adata.X > 0).sum(1)

def filterCells(adata, min_counts=None, min_genes=None, max_counts=None, max_genes=None):
    if min_counts is not None:
        sc.pp.filter_cells(adata, min_counts=min_counts)
    if min_genes is not None:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    if max_genes is not None:
        sc.pp.filter_cells(adata, max_genes=max_genes)

def filterGenes(adata, min_counts=None, min_cells=None, max_counts=None, max_cells=None):
    if min_counts is not None:
        sc.pp.filter_genes(adata, min_counts=min_counts)
    if min_cells is not None:
        sc.pp.filter_genes(adata, min_cells=min_cells)
    if max_counts is not None:
        sc.pp.filter_genes(adata, max_counts=max_counts)
    if max_cells is not None:
        sc.pp.filter_genes(adata, max_cells=max_cells)

def findHighlyVariableGenes(adata, n_top_genes, flavor='cell_ranger', log=False):
    sc.pp.filter_genes_dispersion(adata, flavor=flavor, n_top_genes=n_top_genes,
        log=log, subset=False)
   
def writeHighlyVariableSubset(adata, out_file, genome):
    if 'highly_variable' not in adata.var.keys():
        print("no field for highly variable genes")
        sys.exit()
    subset =  adata[:, adata.var['highly_variable']]
    write10xH5(subset, out_file, genome)
