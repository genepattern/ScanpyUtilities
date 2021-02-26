# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
FROM r-base:4.0.3

MAINTAINER Ted Liefeld <jliefeld@cloud.ucsd.edu>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN mkdir /build

# install system dependencies
RUN apt-get update --yes
RUN apt-get install build-essential --yes
RUN apt-get install libcurl4-gnutls-dev --yes
RUN apt-get install libhdf5-serial-dev --yes
# RUN apt-get install libigraph0-dev --yes #This should install automatically with python-igraph as the repo version fails
RUN apt-get install libxml2-dev --yes

# install python with conda
RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"

# install R dependencies
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install()'
RUN R -e 'BiocManager::install("XML")'
RUN R -e 'BiocManager::install("matrixStats")'
RUN R -e 'BiocManager::install("RSQLite")'
RUN R -e 'BiocManager::install("rhdf5")'
RUN R -e 'BiocManager::install("igraph")'
RUN R -e 'BiocManager::install("sva")'
RUN R -e 'BiocManager::install("scran")'
RUN R -e 'BiocManager::install("monocle")'
RUN R -e 'BiocManager::install("DelayedArray")'
RUN R -e 'BiocManager::install("DelayedMatrixStats")'
RUN R -e 'BiocManager::install("org.Hs.eg.db")'
RUN R -e 'BiocManager::install("org.Mm.eg.db")'
RUN R -e 'remotes::install_github("cole-trapnell-lab/garnett")'

# install python dependencies
RUN pip install numpy==1.19.2
RUN pip install pandas==1.2.2
RUN pip install scipy==1.5.3
RUN pip install anndata==0.7.6
RUN pip install python-igraph==0.8.2
RUN pip install louvain==0.7.0
RUN pip install scanpy==1.7
RUN pip install cmake==3.18.4
RUN pip install MulticoreTSNE==0.1
RUN pip install loompy==3.0.6

# copy module files
COPY module/* /build/
RUN chmod a+x /build/run_module.sh

# display software versions
RUN python --version
RUN pip --version
RUN R --version

# default command
CMD ["python --version"]
