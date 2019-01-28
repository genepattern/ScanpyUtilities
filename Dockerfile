# copyright 2017-2018 Regents of the University of California and the Broad Institute. All rights reserved.
FROM r-base:3.5.2

MAINTAINER Ted Liefeld <jliefeld@cloud.ucsd.edu>

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

#RUN apt-get update && \
#   apt-get install zip --yes && \
#   apt-get install software-properties-common --yes

RUN apt-get update  --yes && \
    apt-get install build-essential --yes && \
    apt-get install python --yes && \
    wget --no-check-certificate https://bootstrap.pypa.io/get-pip.py && \
    python get-pip.py


RUN mkdir /build

RUN mkdir /conda && \
    cd /conda && \
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -b -p /opt/conda

ENV PATH="/opt/conda/bin:/build:${PATH}"

RUN apt-get install libhdf5-serial-dev --yes && \
    apt-get install -y libigraph0-dev

ADD r-installs.R /build/r-installs.R
ADD requirements.txt /build/requirements.txt
RUN pip install -r /build/requirements.txt
RUN Rscript /build/r-installs.R

COPY module/* /build/
RUN chmod a+x /build/run_module.sh

ENV PYTHONPATH /build:$PYTHONPATH

CMD [ "python --version"]
