import sys
import numpy as np
import scanpy as sc

input_file = sys.argv[1]
output_file = sys.argv[2]

print("clustering input: ", input_file)
adata = sc.read(input_file)

print("-- computing clusters with method " + sys.argv[3] + " and resolution " + str(sys.argv[4]) + " --")

if sys.argv[3] == "louvain":
    sc.tl.louvain(adata, resolution=float(sys.argv[4]))

if sys.argv[3] == "leiden":
    sc.tl.leiden(adata, resolution=float(sys.argv[4]))

print("clustering output: ", output_file)
adata.write(output_file, compression='gzip')#, compression_opts=1)
