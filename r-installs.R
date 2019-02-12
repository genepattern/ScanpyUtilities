#install.packages(c("Matrix", "scran","sva"), repos='http://cran.us.r-project.org')
#suppressMessages(library(Matrix))
#suppressMessages(library(rhdf5))
#suppressMessages(library(scran))
#suppressMessages(library(sva))

install.packages("BiocManager")
source("http://bioconductor.org/biocLite.R")
bioc_packages <- c("rhdf5", "igraph", "sva", "scran")
#biocLite()
#biocLite(bioc_packages)
BiocManager::install(bioc_packages)

sessionInfo()



