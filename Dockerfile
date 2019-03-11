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
    apt-get install libigraph0-dev --yes && \
    apt-get install python --yes && \
    wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py

# set up conda environment
RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda

ENV PATH="/opt/conda/bin:/build:${PATH}"

# install python dependencies
ADD requirements.txt /build/requirements.txt
RUN pip install -r /build/requirements.txt

# install R dependencies
RUN R -e 'install.packages("remotes")'
RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install()'
RUN R -e 'BiocManager::install("rhdf5")'
RUN R -e 'BiocManager::install("igraph")'
RUN R -e 'BiocManager::install("sva")'
RUN R -e 'BiocManager::install("scran")'

COPY module/* /build/
RUN chmod a+x /build/run_module.sh

ENV PYTHONPATH /build:$PYTHONPATH

CMD ["python --version"]
