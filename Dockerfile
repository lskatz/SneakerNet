FROM ubuntu:xenial

LABEL base.image="ubuntu:xenial"
LABEL container.version="1.0.0"
LABEL software="SneakerNet"
LABEL software.version="0.8.8"
LABEL description="QA/QC pipeline for a MiSeq/HiSeq/Ion Torrent run"
LABEL website="https://github.com/lskatz/SneakerNet"
LABEL license="https://github.com/lskatz/SneakerNet/blob/master/LICENSE"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="pjx8@cdc.gov"

### install dependencies ###
RUN apt-get update && apt-get -y --no-install-recommends install \
 perl \
 rsync \
 ssh \
 wget \
 cpanminus \
 make \
 prodigal \
 libmoo-perl \
 liblist-moreutils-perl \
 libjson-perl \
 gzip \
 file \

 && apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/* 

# install perl modules
# MLST:
# 
RUN cpanm .....

# CG-Pipeline
RUN

# Kraken1
RUN

# Krona


# Skesa 2.3.0 - DL binary and rename as 'skesa'
RUN mkdir skesa && \
  cd skesa && \
  wget https://github.com/ncbi/SKESA/releases/download/v2.3.0/skesa.centos6.10 && \
  mv skesa.centos6.10 skesa && \
  chmod +x skesa

# Prodigal - 2.6.2 via apt

# Shovill


# mlst
# dependencies in apt: libmoo-perl liblist-moreutils-perl libjson-perl gzip file
# other dependencies: any2fasta ncbi-blast+
RUN wget https://github.com/tseemann/mlst/archive/v2.19.0.tar.gz && \
 tar -xzf v2.19.0.tar.gz && \
 rm v2.19.0.tar.gz

# any2fasta
RUN cd /usr/local/bin && \
  wget https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta && \
  chmod +x any2fasta

# ncbi-blast+
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
 tar -xzf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
 rm ncbi-blast-2.9.0+-x64-linux.tar.gz

# ColorID


# Get SneakerNet 0.8.8 and make /data
RUN wget https://github.com/lskatz/SneakerNet/archive/v0.8.8.tar.gz && \
  tar -zxf v0.8.8.tar.gz && \
  rm v0.8.8.tar.gz && \
  mkdir /data

# download minikraken db
RUN 

# set PATH and perl local settings
ENV PATH="${PATH}:\
/mlst-2.19.0/bin:\
/ncbi-blast-2.9.0+/bin:\
/skesa:\
" \
    LC_ALL=C

WORKDIR /data
