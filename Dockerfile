# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
# builds genepattern/scanpyutilities:0.100
FROM r-base:4.0.3

MAINTAINER Ted Liefeld <jliefeld@cloud.ucsd.edu>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN mkdir /build

# for install log files - check here for log files when debugging
RUN mkdir /logs

# install system dependencies
RUN apt-get update --yes
RUN apt-get install build-essential=12.9 --yes | tee /logs/build-essential_install.log
# 7.74.0-1.3+b1 was listed as version installed, but can't apt-get install this version - it's just curl, so moving on - see log for details
RUN apt-get install libcurl4-gnutls-dev --yes | tee /logs/libcurl4-gnutls-dev_install.log
RUN apt-get install libhdf5-serial-dev --yes | tee /logs/llibhdf5-serial-dev_install.log
# RUN apt-get install libhdf5-dev=1.10.6+repack-5 --yes << wouldn't install by version - see log above for details
#RUN apt-get install libigraph0-dev --yes
RUN apt-get install libxml2-dev --yes | tee /logs/libxml2-dev_install.log
# RUN apt-get install libxml2-dev=2.9.10+dfsg-6.7 --yes << wouldn't install by version see log above for details
RUN apt-get install libtool=2.4.6-15 --yes | tee /logs/libtool_install.log
RUN apt-get install flex=2.6.4-8 --yes | tee /logs/flex_install.log
RUN apt-get install bison --yes | tee /logs/bison_install.log
# RUN apt-get install bison=2:3.7.6+dfsg-1 --yes << wouldn't install by version see log above for details

# install python with conda
RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py38_4.10.3-Linux-x86_64.sh && \
    bash Miniconda3-py38_4.10.3-Linux-x86_64.sh -b -p /opt/conda
ENV PATH="/opt/conda/bin:${PATH}"

# install R dependencies
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/remotes/remotes_2.4.0.tar.gz', repo=NULL, type='source')" | tee /logs/remotes_install.log
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.16.tar.gz', repo=NULL, type='source')" | tee /logs/BiocManager_install.log
# RUN R -e "'Bioc"Manager::install()'
RUN R -e "BiocManager::install('XML', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('matrixStats', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('RSQLite', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('rhdf5', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('igraph', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('sva', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('scran', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('monocle', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('DelayedArray', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('DelayedMatrixStats', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('org.Hs.eg.db', version = '3.12', ask = FALSE)"
RUN R -e "BiocManager::install('org.Mm.eg.db', version = '3.12', ask = FALSE)"
#this should probably be locked down to a tag, but I'm not sure which is correct...

RUN R -e "remotes::install_github('cole-trapnell-lab/garnett', ref='monocle3')"

RUN echo "Here goes"

# install python dependencies
# cython is new addition 9/21/21
RUN pip install Cython==0.29.24
RUN pip install numba==0.52.0
RUN pip install numpy==1.20.3
RUN pip install pandas==1.2.2
RUN pip install scipy==1.7.1
RUN pip install anndata==0.7.6
RUN pip install python-igraph==0.9.6
RUN pip install louvain==0.7.0
RUN pip install scanpy==1.8.1
RUN pip install cmake==3.18.2
RUN pip install MulticoreTSNE==0.1
RUN pip install loompy==3.0.6
RUN pip install leidenalg==0.8.7

# copy module files
COPY module/* /build/
RUN chmod a+x /build/run_module.sh

# display software versions
RUN python --version
RUN pip --version
RUN R --version

# default command
CMD ["python", "--version"]
