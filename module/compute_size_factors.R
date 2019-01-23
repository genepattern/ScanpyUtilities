suppressMessages(library(Matrix))
suppressMessages(library(rhdf5))
suppressMessages(library(scran))

print("computing size factors with scran")

x <- h5read("temp_clustered_for_scran.h5ad", "/obs")
groups <- x$groups
names(groups) <- x$index

raw_data <- as.numeric(h5read("temp_clustered_for_scran.h5ad", "/X/data"))
indices <- as.integer(h5read("temp_clustered_for_scran.h5ad", "/X/indices"))
indptr <- as.integer(h5read("temp_clustered_for_scran.h5ad", "/X/indptr"))

data <- sparseMatrix(i=indices+1, p=indptr, x=raw_data)

size_factors <- computeSumFactors(data, clusters=groups, min.mean=0.1,
    BPPARAM=BiocParallel::MulticoreParam())
names(size_factors) <- rownames(groups)
write.csv(size_factors, file="temp_size_factors.csv")