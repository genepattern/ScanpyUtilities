suppressMessages(library(Matrix))
suppressMessages(library(rhdf5))
suppressMessages(library(monocle3))
suppressMessages(library(garnett))

# read commnad line arguments
args <- commandArgs(trailingOnly = TRUE)
adata_file <- args[1]
marker_file <- args[2]
gene_anno_db <- args[3]
marker_gene_id_type <- args[4]

if (!is.na(gene_anno_db)) {
 print(paste("gene annotation database:", gene_anno_db))
 if (!(gene_anno_db %in% c("org.Hs.eg.db", "org.Mm.eg.db")))
  stop("unrecognized gene annotation database")
 suppressMessages(library(gene_anno_db, character.only = TRUE))
}

# read count data from file
raw_data <- as.numeric(h5read(adata_file, "/X/data"))
indices <- as.integer(h5read(adata_file, "/X/indices"))
indptr <- as.integer(h5read(adata_file, "/X/indptr"))
counts <- sparseMatrix(i = indices + 1, p = indptr, x = raw_data)

# read sample info from file
obs <- h5read(adata_file, "/obs")
obs_names <- as.data.frame(obs$"_index")
rownames(obs_names) <- obs_names[, 1]
obs_names <- obs_names[, -1]

# read gene info from file
var <- h5read(adata_file, "/var")
var_names <- as.data.frame(var$"_index")
if ("gene_name" %in% names(var$"__categories")) {
var_gene_short_name <- as.data.frame(var$"__categories")
var_gene_short_name <- as.data.frame(var_gene_short_name$gene_name[(var$gene_name) + 1])
colnames(var_gene_short_name) = "gene_short_name"
var_names <- cbind(var_names, var_gene_short_name)
}
if (all(var_names[, 1][grep("^ENS", var_names[, 1])] == var_names[, 1])) {
 var_names[, 1] <- gsub("\\..*", "", var_names[, 1])
}
rownames(var_names) <- var_names[, 1]
var_names <- var_names[, -1 , drop = FALSE]

# create cds
cds <- new_cell_data_set(counts, cell_metadata = obs_names, gene_metadata = var_names)
cds <- estimate_size_factors(cds)


# If Not Using a Pretrained Classifier
if (tools::file_ext(marker_file) == "txt") {
 # evaluate quality of marker file and return results in an image ##not working
 if (is.na(gene_anno_db)) {
  marker_check <- check_markers(cds, marker_file, db = "none")
 } else {
  marker_check <- check_markers(cds, marker_file, db = get(gene_anno_db), cds_gene_id_type = "ENSEMBL",
   marker_file_gene_id_type = marker_gene_id_type)
 }
 png("ica_marker_check.png")
 plot_markers(marker_check)
 dev.off()
}

# train classifier and apply to data set
if (is.na(gene_anno_db)) {
 if (tools::file_ext(marker_file) == "txt") {

  classifier <- train_cell_classifier(cds = cds, marker_file = marker_file,
   db = "none")
 } else {
  cat("Reading pretrained classifier from RDS...\n")
  classifier <- readRDS(marker_file)
  cat("Done!\n")
 }

 cds <- classify_cells(cds, classifier, cluster_extend = TRUE, db = "none")
} else {
 if (tools::file_ext(marker_file) == "txt") {
  classifier <- train_cell_classifier(cds = cds, marker_file = marker_file,
   db = get(gene_anno_db), cds_gene_id_type = "ENSEMBL", marker_file_gene_id_type = marker_gene_id_type)
 } else {
  cat("Reading pretrained classifier from RDS...\n")
  classifier <- readRDS(marker_file)
  cat("Done!\n")
 }
 cds <- classify_cells(cds, classifier, db = get(gene_anno_db), cluster_extend = TRUE,
  cds_gene_id_type = "ENSEMBL")
}
print(table(pData(cds)$cell_type))

# write cell types back to adata file
obs$cell_type <- as.character(pData(cds)$cell_type)
obs$extended_cell_type <- as.character(pData(cds)$cluster_ext_type)
obs$index <- rownames(obs)
h5delete(adata_file, "/obs")
h5write(obs, adata_file, "/obs")
