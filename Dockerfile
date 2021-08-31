FROM ubuntu:bionic

LABEL base.image="ubuntu:bionic"
LABEL dockerfile.version="1.0.0"
LABEL software="SneakerNet"
LABEL software.version="0.20.1"
LABEL description="SneakerNet QA/QC system for primary genomic data"
LABEL website="https://github.com/lskatz/SneakerNet"
LABEL license="https://github.com/lskatz/SneakerNet/LICENSE.md"
LABEL maintainer1="Lee Katz"
LABEL maintainer1.email="gzu2@cdc.gov"
LABEL maintainer2="Curtis Kapsak"
LABEL maintainer2.email="pjx8@cdc.gov"

# https://github.com/StaPH-B/docker-builds/blob/master/mash/2.2/Dockerfile
FROM staphb/mash:2.2 AS mash
# https://github.com/StaPH-B/docker-builds/blob/master/skesa/2.4.0/Dockerfile
FROM staphb/skesa:2.4.0 AS skesa
# etc on the URL sources
FROM staphb/kraken:1.1.1-no-db AS kraken
FROM staphb/mlst:2.19.0 AS mlst
FROM staphb/prokka:1.14.5 AS prokka
FROM staphb/shovill:1.1.0 AS shovill
FROM staphb/seqtk:1.3 AS seqtk
#FROM staphb/staramr:0.7.1 AS staramr
FROM staphb/salmid:0.1.23 AS salmid
#FROM rust:1.46.0-slim-buster AS rust

# Other sources
FROM mgibio/samtools:1.9 AS samtools
FROM flowcraft/krona:2.7-1 AS krona
#FROM mickaelsilva/chewbbaca_py3 AS chewbbaca

# EDIT: this bioperl container uses perl/5.18 which doesn't match our perl v5.26.1
# Bring in libraries
#FROM bioperl/bioperl:release-1-7-2 AS bioperl

# "Import" ubuntu:bionic one more time as a hack.
# I hope this step kind of clears the slate so that the COPY
# commands don't fubar.
FROM ubuntu:bionic

# Let's get all the executables into our path
COPY --from=mash      /mash-Linux64-v2.2/* /usr/local/bin/
COPY --from=skesa     /skesa/*             /usr/local/bin/
COPY --from=kraken    /opt/kraken/*        /usr/local/bin/
COPY --from=kraken    /opt/bin/*           /usr/local/bin/
COPY --from=mlst      /mlst-2.19.0         /mlst-2.19.0    
COPY --from=mlst      /usr/local/bin/any2fasta              /usr/local/bin/
COPY --from=prokka    /bedtools2/bin       /usr/local/bin/
COPY --from=prokka    /prokka-1.14.5/bin   /usr/local/bin/
COPY --from=prokka    /barrnap-0.9/bin     /usr/local/bin/
COPY --from=seqtk     /seqtk-1.3           /usr/local/bin/
#COPY --from=staramr   /usr/local/bin       /usr/local/bin/
#COPY --from=salmid    /usr/local/bin       /usr/local/bin/
COPY --from=samtools  /opt/samtools/bin/samtools            /usr/local/bin/
COPY --from=krona     /NGStools/KronaTools-2.7              /NGStools/KronaTools-2.7
#COPY --from=chewbbaca /NGStools/clustalw-2.1-linux-x86_64-libcppstatic  /NGStools
#COPY --from=chewbbaca /NGStools/Prodigal                                /NGStools
#COPY --from=chewbbaca /NGStools/prodigal_training_files                 /NGStools
#COPY --from=chewbbaca /usr/local/bin/*     /usr/local/bin/
COPY --from=mlst      /ncbi-blast-2.9.0+   /ncbi-blast-2.9.0+/

#COPY --from=rust      /usr/local/rustup    /usr/local/rustup
#COPY --from=rust      /usr/local/cargo     /usr/local/cargo

# Libraries
#COPY --from=staramr   /usr/local/lib/python3.6             /usr/local/lib/python3.6/
#COPY --from=blast     /lib/x86_64-linux-gnu /lib/x86_64-linux-gnu
#COPY --from=blast     /lib64               /lib64
#COPY --from=blast     /usr/lib/x86_64-linux-gnu  /usr/lib/x86_64-linux-gnu

# Taking a risk using python3.5 libraries in a python3.6 folder
#COPY --from=salmid    /usr/local/lib/python3.5             /usr/local/lib/python3.6/
#COPY --from=chewbbaca /usr/local/lib/python3.5             /usr/local/lib/python3.6/
#COPY --from=bioperl   /usr/lib/            /usr/lib
#COPY --from=bioperl   /usr/local/lib/      /usr/local/lib
#COPY --from=bioperl   /usr/share           /usr/share

# Tons of Shovill executables
COPY --from=shovill   /usr/local/bin/       /usr/local/bin/
COPY --from=shovill   /SPAdes-3.14.1-Linux  /SPAdes-3.14.1-Linux
COPY --from=shovill   /kmc                 /usr/local/bin/
COPY --from=shovill   /Lighter-1.1.1       /usr/local/bin/
COPY --from=shovill   /trimmomatic         /trimmomatic
COPY --from=shovill   /bwa/bwa-0.7.17      /usr/local/bin/
COPY --from=shovill   /megahit/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin       /usr/local/bin/
COPY --from=shovill   /velvet/velvet-1.2.10                 /usr/local/bin/
COPY --from=shovill   /flash/FLASH-1.2.11  /usr/local/bin/
COPY --from=shovill   /shovill/shovill-1.1.0                /shovill/shovill-1.1.0
COPY --from=shovill   /pilon               /pilon
COPY --from=shovill   /samclip             /usr/local/bin/


# System installations that I think I need
# https://gist.github.com/ryanwilsonperkin/0daf26385813196291bf32492802a4ca
# https://github.com/ilikenwf/apt-fast
#RUN echo deb http://ppa.launchpad.net/apt-fast/stable/ubuntu bionic main >> /etc/apt/sources.list.d/apt-fast.list && \
# echo deb-src http://ppa.launchpad.net/apt-fast/stable/ubuntu bionic main >> /etc/apt/sources.list.d/apt-fast.list && \
# apt-key adv --keyserver keyserver.ubuntu.com --recv-keys A2166B8DE8BDC3367D1901C11EE2FF37CA8DA16B && \
# apt-get update && \
# apt-get install -y apt-fast

RUN apt-get update && \
 apt-get  install -y --no-install-recommends \
 perl \
 git \
 rsync \
 vim \
 less \
 ssh \
 wget \
 curl \
 bsdmainutils \
 cpanminus \
 make \
 prodigal \
 libmoo-perl \
 liblist-moreutils-perl \
 libjson-perl \
 gzip \
 file \
 zlib1g-dev \
 g++ \ 
 gawk \
 libexpat1-dev \
 sendmail \ 
 zip \
 python \
 python-setuptools \
 build-essential \
 python3 \
 python3-pip \
 python3-setuptools \
 python3-venv \
 pigz \
 gcc \
 libpthread-stubs0-dev \
  default-jre \
 unzip \
 bzip2 \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 libssl-dev \
 libfindbin-libs-perl \
 psmisc \
 libatlas-base-dev \
 mafft \
 libpython3-dev \
 locales \
 && \
 apt-get autoclean && rm -rf /var/lib/apt/lists/* 
#python-matplotlib ipython python-pandas python-sympy python-nose
# python-numpy python-scipy  \

# separate installation so that I can take advantage of cache while testing
RUN apt-get update && \
 apt-get install -y libgd-perl libgd-dev

# Set LC_ALL env
# https://github.com/hpcng/singularity/issues/11#issuecomment-325235446
RUN echo "LC_ALL=en_US.UTF-8" >> /etc/environment && \
 echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen && \
 echo "LANG=en_US.UTF-8" > /etc/locale.conf && \
 locale-gen en_US.UTF-8


# Perl libraries
RUN cpanm --force --notest \
  LWP::Protocol::https \
  IO::Socket::SSL \
  DBD::SQLite \
  DBI \
  File::Which 

# Unsure if I need these
#  Bio::Sketch::Mash \
#  Bio::Kmer \
 
# All the manual installations I have to do go after this line

RUN git clone https://github.com/lskatz/CG-Pipeline.git

# SneakerNet version (corresponds to github version and $VERSION)
# TODO I'd rather set this in .gitlab-ci.yml and so here is a parameter expansion format
#ENV SNVER=${SNVER:-0.20.1}

# Get SneakerNet and make /data
# apt deps: sendmail-base zip bsdmainutils (for column command)
# perl modules listed in cpanm comments above (some installed there, remaining installed w cpanm command below)
# 
# The command at the end t/00_env.t sets up an environment including generating contaminated reads
# for further unit testing.
RUN export SNVER=${SNVER:-0.20.1} && \
 echo "SNVER is $SNVER" && \
 wget https://github.com/lskatz/SneakerNet/archive/v${SNVER}.tar.gz && \
 tar -zxf v${SNVER}.tar.gz && \
 rm v${SNVER}.tar.gz && \
 cd /SneakerNet-${SNVER} && \
 cpanm --installdeps --notest --force . && \
 perl Makefile.PL && \
 make && \
 sed -i 's+/opt/kraken/full-20140723+/kraken-database+g' config/settings.conf && \
 mkdir /data /kraken-database && \
 perl t/00_env.t

ENV PATH="${PATH}:\
/CG-Pipeline/scripts:\
/KronaTools-2.7.1/bin:\
/SneakerNet-${SNVER:-0.20.1}/scripts:/SneakerNet-${SNVER:-0.20.1}/SneakerNet.plugins:\
/mlst-2.19.0/bin/:\
/NGStools/KronaTools-2.7/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:\
/NGStools/chewBBACA:/NGStools/chewBBACA/utils:/NGStools/prodigal_training_files:/NGStools/clustalw-2.1-linux-x86_64-libcppstatic:\
/usr/local/bin/Trimmomatic-0.38:\
/colorid:\
/pilon:\
/trimmomatic:\
/usr/local/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:\
/shovill/shovill-1.1.0/bin:\
/SPAdes-3.14.1-Linux/bin:\
/ncbi-blast-2.9.0+/bin:\
"\
 LC_ALL=en_US.UTF-8 \
 RUSTUP_HOME=/usr/local/rustup CARGO_HOME=/usr/local/cargo RUST_VERSION=1.46.0 \
 BLASTDB=/blast/blastdb 

## pip installations after this line
RUN python3 -m pip install --upgrade pip

# Staramr: lifting the code from staphb
# https://github.com/StaPH-B/docker-builds/blob/master/staramr/0.7.1/Dockerfile
RUN pip3 install staramr==0.7.1 pandas==0.25.3 && \
  staramr db update -d && \
  staramr db info

# Pip installations after I set the path
# SalmID 0.1.23
# apt deps: python-setuptools python3 python3-pip curl build-essential file git python3-venv
RUN pip3 install poetry && \
 git clone https://github.com/hcdenbakker/SalmID.git --branch 0.1.23 --single-branch && \
 cd SalmID && \
 poetry build -vvv && \
 pip3 install dist/salmid*.whl

# Chewbbaca
# Taking code from here: https://hub.docker.com/r/mickaelsilva/chewbbaca_py3/dockerfile
WORKDIR /NGStools/
#GET training files and Prodigal 
RUN git clone https://github.com/hyattpd/Prodigal.git && \
  pip3 install biopython plotly SPARQLWrapper chewbbaca && \
  cd /NGStools/Prodigal && \
  make install
WORKDIR /NGStools/
RUN git clone https://github.com/mickaelsilva/prodigal_training_files && \
  wget www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz && \
  tar -zxf clustalw-2.1-linux-x86_64-libcppstatic.tar.gz && \
  rm clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
# Reset the working directory after chewbbaca installation
WORKDIR /

# Rust for at least colorid
#RUN curl https://sh.rustup.rs -sSf | bash -s -- -y

# colorid: the correct version is not versioned but is at this hash
RUN mkdir colorid && \
  cd colorid && \
  wget --no-check-certificate https://github.com/hcdenbakker/colorid/releases/download/v0.1.4.3/colorid_Linux64v0.1.4.3 && \
  chmod +x colorid_Linux64v0.1.4.3 && \
  mv colorid_Linux64v0.1.4.3 /usr/local/bin/colorid

# Trying to avoid an error where LC_ALL gets somehow undefined before this step
#   bash: warning: setlocale: LC_ALL: cannot change locale (en_US.UTF-8)
ENV LC_ALL=en_US.UTF-8

WORKDIR /data

