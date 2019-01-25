import sys
import scanpy.api as sc

input_file = sys.argv[1]
output_basename = sys.argv[2]
n_genes = int(sys.argv[3])

anno_outfile = output_basename + "_high_variance_genes_annotated.h5ad"
subset_outfile = output_basename + "_high_variance_genes_subset.h5ad"
print("high variance genes output (annotated): ", anno_outfile)
print("high variance genes output (subsetted): ", subset_outfile)

adata = sc.read(input_file)
sc.pp.filter_genes_dispersion(adata, flavor='cell_ranger',
    n_top_genes=n_genes, log=False, subset=False)
adata.write(anno_outfile, compression='gzip', compression_opts=1)
subset = adata[:, adata.var['highly_variable']]
subset.write(subset_outfile, compression='gzip', compression_opts=1)