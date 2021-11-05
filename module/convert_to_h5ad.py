# datafile, outfilename, source file type, and genome are the args
# sys.argv[1]

import os, sys
import scanpy as sc
if sys.argv[1] == "h5":
    adata = sc.read_10x_h5(sys.argv[2], genome=sys.argv[4])

if sys.argv[1] == "loom":
    adata = sc.read_loom(sys.argv[2])

adata.write(sys.argv[3], compression='gzip')#, compression_opts=1)
