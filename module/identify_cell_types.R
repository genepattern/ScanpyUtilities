library(Matrix)
library(rhdf5)
library(monocle)

# read commnad line arguments
args <- commandArgs(trailingOnly=TRUE)
adata_file <- args[1]
marker_file <- args[2]
gene_anno_db <- args[3]
if (!is.null(geno_anno_db) & !(gene_anno_db %in% c("org.Hs.eg.db", "org.Mm.eg.db")))
    stop("unrecognized gene annotation database")
library(gene_anno_db)

# read count data from file
raw_data <- as.numeric(h5read(adata_file, "/X/data"))
indices <- as.integer(h5read(adata_file, "/X/indices"))
indptr <- as.integer(h5read(adata_file, "/X/indptr"))
counts <- sparseMatrix(i=indices+1, p=indptr, x=raw_data)

# read sample info from file
obs <- h5read(adata_file, "/obs")
rownames(obs) <- obs[,1]
obs <- obs[,-1]
pd <- AnnotatedDataFrame(data=obs)

# read gene info from file
var <- h5read(adata_file, "/var")
rownames(var) <- var[,1]
var <- var[,-1]
fd <- new("AnnotatedDataFrame", data=var)

# create cds
cds <- newCellDataSet(counts, phenoData=pd, featureData=fd,
    expressionFamily=negbinomial.size())
cds <- estimateSizeFactors(cds)

# evaluate quality of marker file and return results in an image
if (is.null(gene_anno_db))
    marker_check <- check_markers(cds, marker_file)
else
    marker_check <- check_markers(cds, marker_file, db=get(gene_anno_db),
        cds_gene_id_type="ENSEMBL", marker_file_gene_id_type="SYMBOL")
png("ica_marker_check.png")
plot_markers(marker_check)
dev.off()

# train classifier and apply to data set
if (is.null(gene_anno_db))
{
    classifier <- train_cell_classifier(cds=cds, marker_file=marker_file)
    cds <- classify_cells(cds, classifier, cluster_extend=TRUE)
}
else
{
    classifier <- train_cell_classifier(cds=cds, marker_file=marker_file,
        db=get(gene_anno_db), cds_gene_id_type="ENSEMBL",
        marker_file_gene_id_type="SYMBOL")
    cds <- classify_cells(cds, classifier, db=get(gene_anno_db),
        cluster_extend=TRUE, cds_gene_id_type="ENSEMBL")
}

# write cell types back to adata file
obs$cell_type <- as.character(pd$cell_type)
h5delete(adata_file, "/obs")
h5write(obs, adata_file, "/obs")
