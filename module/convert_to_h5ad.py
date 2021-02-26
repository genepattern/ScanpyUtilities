# datafile, genome, outfilename are the args
# sys.argv[1]

import sys
import scanpy as sc
adata = sc.read_10x_h5(sys.argv[1], genome=sys.argv[2])
adata.write(sys.argv[3], compression='gzip')#, compression_opts=1)
