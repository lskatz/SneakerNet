FROM ubuntu:xenial

LABEL base.image="ubuntu:xenial"
LABEL container.version="1"
LABEL software="SneakerNet"
LABEL software.version="0.8.8"
LABEL description="QA/QC pipeline for a MiSeq/HiSeq/Ion Torrent run"
LABEL website="https://github.com/lskatz/SneakerNet"
LABEL license="https://github.com/lskatz/SneakerNet/blob/master/LICENSE"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="pjx8@cdc.gov"

### install dependencies ###
# vim installed temporarily for debugging
# libexpat1-dev needed for cpanm to install XML::Parser and other modules
RUN apt-get update && apt-get -y --no-install-recommends install \
 perl \
 git \
 rsync \
 ssh \
 wget \
 curl \
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
 bioperl \
 libexpat1-dev \
 sendmail \ 
 zip \
 python \
 python3 \
 python3-pip \
 python3-setuptools \
 vim && \
 apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/* 

### install perl modules as root w cpanm ###
# MLST:
# kraken1: Getopt::Std
# SN: Config::Simple local::lib version::vpp Bio::FeatureIO XML::DOM XML::Parser XML::DOM::XPath
RUN cpanm --notest --force Getopt::Std \
  Config::Simple \
  local::lib \
  version::vpp \
  XML::DOM \
  XML::Parser \
  XML::DOM::XPath \
  Bio::FeatureIO 

# CG-Pipeline 
# again can't download because CDC spoofs security certificates which causes errors...GRRRRR
RUN GIT_SSL_NO_VERIFY=true git clone https://github.com/lskatz/CG-Pipeline.git

# Jellyfish - kraken dep
# apt deps: gawk
RUN wget --no-check-certificate https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-1.1.12.tar.gz && \
  tar -zxf jellyfish-1.1.12.tar.gz && \
  rm -rf jellyfish-1.1.12.tar.gz && \
  cd jellyfish-1.1.12 && \
  ./configure --prefix=/opt/ && \
  make -j 4 && \
  make install

# Kraken1
# apt deps: wget zlib1g-dev make g++ rsync cpanminus
# cpan deps: Getopt::Std
RUN wget --no-check-certificate https://github.com/DerrickWood/kraken/archive/v1.1.1.tar.gz && \
  tar -xzf v1.1.1.tar.gz && \
  rm -rf v1.1.1.tar.gz && \
  cd kraken-1.1.1 && \
  mkdir /opt/kraken && \
  ./install_kraken.sh /opt/kraken/

# Krona
# apt deps: curl
RUN wget --no-check-certificate https://github.com/marbl/Krona/releases/download/v2.7.1/KronaTools-2.7.1.tar && \
    tar -xf KronaTools-2.7.1.tar && \
    rm KronaTools-2.7.1.tar && \
    cd KronaTools-2.7.1 && \
    ./install.pl --prefix . && \
    ./updateTaxonomy.sh

# Skesa 2.3.0 - DL binary and rename as 'skesa'
RUN mkdir skesa && \
  cd skesa && \
  wget --no-check-certificate https://github.com/ncbi/SKESA/releases/download/v2.3.0/skesa.centos6.10 && \
  mv skesa.centos6.10 skesa && \
  chmod +x skesa

# Prodigal - 2.6.2 via apt

# Shovill
# only used in crypto


# SPAdes 3.14.0 (needed for shovill and metagenomics workflow)
# apt deps: python
RUN wget http://cab.spbu.ru/files/release3.14.0/SPAdes-3.14.0-Linux.tar.gz && \
  tar -xzf SPAdes-3.14.0-Linux.tar.gz && \
  rm -r SPAdes-3.14.0-Linux.tar.gz

# mlst 2.16.2 (had to downgrade since later versions of mlst require perl 5.26.0 which is not available on apt for ubuntu:xenial)
# dependencies in apt: libmoo-perl liblist-moreutils-perl libjson-perl gzip file
# other dependencies: any2fasta ncbi-blast+
RUN wget --no-check-certificate https://github.com/tseemann/mlst/archive/v2.16.2.tar.gz &&\
 tar -xzf v2.16.2.tar.gz &&\
 rm v2.16.2.tar.gz

# any2fasta
RUN cd /usr/local/bin && \
  wget --no-check-certificate https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta && \
  chmod +x any2fasta

# ncbi-blast+ 2.9.0
RUN wget --no-check-certificate ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
 tar -xzf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
 rm ncbi-blast-2.9.0+-x64-linux.tar.gz

# staramr
# apt deps: python3 python3-pip python3-setuptools git ncbi-blast+ (blast installed manually)
# update pip3 and install staramr 0.5.1
RUN python3 -m pip install -U pip && \
 pip3 install staramr==0.5.1

# ColorID
# precompiled binary here https://github.com/hcdenbakker/colorid/releases

# Get SneakerNet 0.8.14 and make /data
# apt deps: sendmail-base zip 
# perl modules listed in cpanm comments above (some installed there, remaining installed w cpanm command below)
RUN wget --no-check-certificate https://github.com/lskatz/SneakerNet/archive/v0.8.14.tar.gz && \
  tar -zxf v0.8.14.tar.gz && \
  rm v0.8.14.tar.gz && \
  cd /SneakerNet-0.8.14 && \
  cpanm --installdeps --notest --force . && \
  mkdir /data

#### TEMP COMMENTED OUT TO SAVE BUILD TIME ####
# minikraken db
#RUN mkdir /kraken-database && \
#  cd /kraken-database && \
#  wget https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz && \
#  tar -zxf minikraken_20171019_4GB.tgz && \
#  rm -rf minikraken_20171019_4GB.tgz 

# set PATH and perl local settings
ENV PATH="${PATH}:\
/mlst-2.16.2/bin:\
/ncbi-blast-2.9.0+/bin:\
/skesa:\
/opt/kraken:\
/opt/bin:\
/CG-Pipeline/scripts:\
/KronaTools-2.7.1/bin:\
/SPAdes-3.14.0-Linux/bin\
" \
    LC_ALL=C \
    ls='ls --color=auto'

# check SN dependencies
RUN ./SneakerNet-0.8.14/scripts/SneakerNet.checkdeps.pl default

WORKDIR /data
