import sys
import scanpy as sc

print("cell filtering input:", sys.argv[1])
print("cell filtering output:", sys.argv[2])

adata = sc.read(sys.argv[1])

min_counts = int(sys.argv[3])
max_counts = int(sys.argv[4])
min_genes = int(sys.argv[5])
max_genes = int(sys.argv[6])
mito_file = sys.argv[7]
mito_pct = int(sys.argv[8])

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
if mito_file != "SKIP":
    with open(mito_file) as f:
        mito_genes = f.read().splitlines()
    mito_genes = list(set([sub.replace('-I', '') for sub in mito_genes]))
    adata.var['mt'] = [x in mito_genes for x in adata.var_names]
    print("filtering out cells with less than", mito_pct, "% mitochondrial fraction")
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    adata = adata[adata.obs.pct_counts_mt < mito_pct, :]
    print("cells remaining: ", adata.shape[0])

adata.write(sys.argv[2], compression='gzip')#, compression_opts=1)
