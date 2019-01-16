import argparse
import distutils.util
from ScanpyPreprocessFunctions import *

def check_positive(value):
    ivalue = int(value)
    if ivalue <= 0:
         raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
    return ivalue



parser = argparse.ArgumentParser()
# ~~~~Module Required Arguments~~~~~ #
parser.add_argument("-f", "--data.file",
                    dest="data_file",
                    type=str,
                    help="H5 Data file of scRNASeq data")
parser.add_argument("-mnc", "--min.counts", dest="min_counts",
                    type=check_positive,
                    help="Filter out samples with fewer than Min Counts",
                    default='750')
parser.add_argument("-mng", "--min.genes", dest="min_genes",
                    type=check_positive,
                    help="Filter out samples with less than min genes",
                    default='180')
parser.add_argument("-mnl", "--min.cells", dest="min_cells",
                    type=check_positive,
                    help="Filter out samples with less than min cells",
                    default='20')

parser.add_argument("-mxc", "--max.counts", dest="max_counts",
                    type=check_positive,
                    help="Filter out samples with more than Max Counts")
parser.add_argument("-mxg", "--max.genes", dest="max_genes",
                    type=check_positive,
                    help="Filter out samples with more than max genes")
parser.add_argument("-mxl", "--max.cells", dest="max_cells",
                    type=check_positive,
                    help="Filter out samples with more than max cells")


parser.add_argument("-o", "--output.filename", dest="output_filename",
                    type=str,
                    help="The basename to use for output file",
                    default='OUT')
parser.add_argument("-g", "--group.name", dest="group_name",
                    type=str,
                    help="The group name in the hdf5 file to use",
                    default='genome')

parser.add_argument("-a", "--annotate",
                    help="Whether to annotate the dataset True/False",
                    type=distutils.util.strtobool)

parser.add_argument("-fc", "--filterCells",
                    help="Whether to filter cells.",
                    type=distutils.util.strtobool)

parser.add_argument("-fg", "--filterGenes",
                    help="Whether to filter genes.",
                    type=distutils.util.strtobool)


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

max_counts = None
if args.max_counts:
   max_counts = args.max_counts

if args.max_cells:
   max_cells = args.max_cells

if args.max_genes:
   max_genes = args.max_genes


   adata = read10xH5(args.data_file, 'GRCh38')


if (args.annotate):
    print("Performing annotation")
    annotateData(adata )

if (args.filterCells):
    print("Filtering cells")
    filterCells(adata,  min_counts=args.min_counts, min_genes=args.min_genes, max_counts=max_counts, max_genes=max_genes)

if (args.filterGenes):
    print("Filtering genes")
    filterGenes(adata, min_cells=args.min_cells, max_cells = max_cells)


fullyFilteredFile = args.output_filename + "_preprocessed.h5"
print("Done - last output is " + fullyFilteredFile)

write10xH5(adata, fullyFilteredFile, 'GRCh38')
