suppressMessages(library(Matrix))
suppressMessages(library(rhdf5))
suppressMessages(library(scran))
suppressMessages(library(sva))

# read command line arguments
args <- commandArgs(TRUE)
input_file = args[1]
batch_variable <- args[2]

# read batch ids
obs <- h5read(input_file, "/obs")
batch_ids <- obs[,batch_variable]

# read raw data from file
raw_values <- as.numeric(h5read(input_file, "/X/data", native=TRUE))
raw_indices <- as.integer(h5read(input_file, "/X/indices", native=TRUE))
raw_indptr <- as.integer(h5read(input_file, "/X/indptr", native=TRUE))

# create matrix
data <- sparseMatrix(i=raw_indices+1, p=raw_indptr, x=raw_values)

# batch correct data
data <- as.matrix(data)
cor_data <- ComBat(dat=data, batch=batch_ids)
cor_data <- as(cor_data, "sparseMatrix")

# delete raw data from h5ad file (we need to rename it to raw.*)
h5delete(input_file, "/X")

# write raw data to h5ad file
h5createGroup(input_file, "raw.X")
h5createDataset(input_file, "/raw.X/data", dim=length(raw_values),
    H5type="H5T_NATIVE_FLOAT", level=1, native=TRUE)
h5createDataset(input_file, "/raw.X/indices", dim=length(raw_indices),
    H5type="H5T_NATIVE_INT", level=1, native=TRUE)
h5createDataset(input_file, "/raw.X/indptr", dim=length(raw_indptr),
    H5type="H5T_NATIVE_INT", level=1, native=TRUE)
h5write(raw_values, file=input_file, name="/raw.X/data")
h5write(raw_indices, file=input_file, name="/raw.X/indices")
h5write(raw_indptr, file=input_file, name="/raw.X/indptr")

# write batch corrected data to h5ad file
h5createGroup(input_file, "X")
h5createDataset(input_file, "/X/data", dim=length(cor_data@x),
    H5type="H5T_NATIVE_FLOAT", level=1, native=TRUE)
h5createDataset(input_file, "/X/indices", dim=length(cor_data@i),
    H5type="H5T_NATIVE_INT", level=1, native=TRUE)
h5createDataset(input_file, "/X/indptr", dim=length(cor_data@p),
    H5type="H5T_NATIVE_INT", level=1, native=TRUE)
h5write(cor_data@x, file=input_file, name="/X/data")
h5write(cor_data@i, file=input_file, name="/X/indices")
h5write(cor_data@p, file=input_file, name="/X/indptr")














test_fun <- function(...)
{
    batches <- list(...)
    return(batches)
}