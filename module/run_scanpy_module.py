import argparse
import distutils.util
import os
import sys
from ScanpyPreprocessFunctions import *

def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue



parser = argparse.ArgumentParser()
# ~~~~Module Required Arguments~~~~~ #
parser.add_argument("-f", "--data.file", dest="data_file",
                    type=str,
                    help="An h5ad file containing the single cell data")

parser.add_argument("-o", "--output.filename", dest="output_filename",
                    type=str,
                    help="Base filename for the output file",
                    default='OUT')

parser.add_argument("-a", "--annotate",
                    help="Whether to annotate the dataset True/False",
                    type=distutils.util.strtobool)

parser.add_argument("-cmnc", "--cells.min.counts", dest="cells_min_counts",
                    type=check_positive,
                    help="Filter out cells with fewer total counts than min counts")

parser.add_argument("-cmxc", "--cells.max.counts", dest="cells_max_counts",
                    type=check_positive,
                    help="Filter out cells with more total counts than max counts")

parser.add_argument("-cmng", "--cells.min.genes", dest="cells_min_genes",
                    type=check_positive,
                    help="Filter out cells with fewer than min genes expressed")

parser.add_argument("-cmxg", "--cells.max.genes", dest="cells_max_genes",
                    type=check_positive,
                    help="Filter out cells with more than max genes expressed")

parser.add_argument("-gmnc", "--genes.min.counts", dest="genes_min_counts",
                    type=check_positive,
                    help="Filter out genes with fewer total counts than min counts")

parser.add_argument("-gmxc", "--genes.max.counts", dest="genes_max_counts",
                    type=check_positive,
                    help="Filter out genes with more total counts than max counts")

parser.add_argument("-gmnc", "--genes.min.cells", dest="genes_min_cells",
                    type=check_positive,
                    help="Filter out genes expressed in fewer than min cells")

parser.add_argument("-gmxc", "--genes.max.cells", dest="genes_max_cells",
                    type=check_positive,
                    help="Filter out genes expressed in more than max cells")

parser.add_argument("-bcn", "--batch.correct.and.normalize", dest="batch_correct_and_normalize",
                    help="Whether to perform batch correction and normalization True/False.  This step is performed after filtering if filtering is true",
                    type=distutils.util.strtobool)


# ~~~~Development Optional Arguments~~~~~ #
# Reminder: "store_true" args, are False by default and when the flag is passed
# then they become True
parser.add_argument("-v", "--verbose",
                    action="store_true",
                    help="increase output verbosity")
parser.add_argument("-d", "--debug",
                    action="store_true",
                    help="increase output verbosity")
args = parser.parse_args()
if args.verbose:
    print("Ah! The old verbosaroo")

print("~~~~~~~~~~~~~~~~~~~~~~")
print("Using arguments:")
print(args)
#print("Now getting work done.")
print("~~~~~~~~~~~~~~~~~~~~~~")

cells_min_counts = None if not args.cells_min_counts else args.cells_min_counts
cells_max_counts = None if not args.cells_max_counts else args.cells_max_counts
cells_min_genes  = None if not args.cells_min_genes  else args.cells_min_genes
cells_max_genes  = None if not args.cells_max_genes  else args.cells_max_genes
genes_min_counts = None if not args.genes_min_counts else args.genes_min_counts
genes_max_counts = None if not args.genes_max_counts else args.genes_max_counts
genes_min_cells  = None if not args.genes_min_cells  else args.genes_min_cells
genes_max_cells  = None if not args.genes_max_cells  else args.genes_max_cells

fname, fext = os.path.splitext(args.data_file)
if fext == '.h5ad':
    adata = sc.read(args.data_file)
else:
    print("invalid file, must be h5ad - see scanpy documentation", file=sys.stderr)
    sys.exit(1)

if (args.annotate):
    print("Performing annotation")
    annotateData(adata)

if (any([cells_min_counts, cells_max_counts, cells_min_genes, cells_max_genes]):
    print("Filtering cells")
    filterCells(adata, min_counts=cells_min_counts, max_counts=cells_max_counts, min_genes=cells_min_genes, max_genes=cells_max_genes)

if (any([genes_min_counts, genes_max_counts, genes_min_cells, genes_max_cells]):
    print("Filtering genes")
    filterGenes(adata, min_counts=genes_min_counts, max_counts=genes_max_counts, min_cells=genes_min_cells, max_cells=genes_max_cells)

fullyFilteredFile = args.output_filename + "_preprocessed.h5ad"
print("Done - last output is " + fullyFilteredFile)
adata.write(fullyFilteredFile, compression='gzip', compression_opts=1)
