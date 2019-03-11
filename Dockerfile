# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
FROM r-base:3.5.2

MAINTAINER Ted Liefeld <jliefeld@cloud.ucsd.edu>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN mkdir /build

# install system dependencies
RUN apt-get update --yes && \
    apt-get install build-essential --yes && \
    apt-get install libcurl4-gnutls-dev --yes && \
    apt-get install libhdf5-serial-dev --yes && \
    apt-get install libigraph0-dev --yes

# install python with conda
RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda

ENV PATH="/opt/conda/bin:${PATH}"

# install python dependencies
RUN pip install numpy==1.15.4
RUN pip install pandas==0.23.4
RUN pip install scipy==1.1.0
RUN pip install anndata==0.6.18
RUN pip install python-igraph==0.7.1.post6
RUN pip install louvain==0.6.1
RUN pip install scanpy==1.3.3

# install R dependencies
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install()'
RUN R -e 'BiocManager::install("rhdf5")'
RUN R -e 'BiocManager::install("igraph")'
RUN R -e 'BiocManager::install("sva")'
RUN R -e 'BiocManager::install("scran")'

# copy module files
COPY module/* /build/
RUN chmod a+x /build/run_module.sh

# display software versions
RUN python --version
RUN pip --version
RUN R --version

# default command
CMD ["python --version"]
