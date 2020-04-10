FROM ubuntu:xenial

LABEL base.image="ubuntu:xenial"
LABEL container.version="1"
LABEL software="SneakerNet"
LABEL software.version="0.9.8"
LABEL description="QA/QC pipeline for a MiSeq/HiSeq/Ion Torrent run"
LABEL website="https://github.com/lskatz/SneakerNet"
LABEL license="https://github.com/lskatz/SneakerNet/blob/master/LICENSE"
LABEL maintainer="Curtis Kapsak"
LABEL maintainer.email="pjx8@cdc.gov"

### install dependencies ###
# libexpat1-dev needed for cpanm to install XML::Parser and other perl modules
RUN apt-get update && apt-get -y --no-install-recommends install \
 perl \
 git \
 rsync \
 vim \
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
 bioperl \
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
 openjdk-9-jre \
 unzip \
 bzip2 \
 libncurses5-dev \
 libbz2-dev \
 liblzma-dev \
 libcurl4-gnutls-dev \
 libssl-dev \
 libfindbin-libs-perl && \
 apt-get clean && apt-get autoclean && rm -rf /var/lib/apt/lists/* 

### install perl modules as root w cpanm ###
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
# let's try without the prefix: GIT_SSL_NO_VERIFY=true
RUN git clone https://github.com/lskatz/CG-Pipeline.git

# Jellyfish 1.1.12 (kraken dep)
# apt deps: gawk
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v1.1.12/jellyfish-1.1.12.tar.gz && \
 tar -zxf jellyfish-1.1.12.tar.gz && \
 rm -rf jellyfish-1.1.12.tar.gz && \
 cd jellyfish-1.1.12 && \
 ./configure --prefix=/opt/ && \
 make -j 4 && \
 make install

# Kraken 1.1.1
# apt deps: wget zlib1g-dev make g++ rsync cpanminus
# cpan deps: Getopt::Std
RUN wget https://github.com/DerrickWood/kraken/archive/v1.1.1.tar.gz && \
 tar -xzf v1.1.1.tar.gz && \
 rm -rf v1.1.1.tar.gz && \
 cd kraken-1.1.1 && \
 mkdir /opt/kraken && \
 ./install_kraken.sh /opt/kraken/

# Krona 2.7.1
# apt deps: curl
RUN wget https://github.com/marbl/Krona/releases/download/v2.7.1/KronaTools-2.7.1.tar && \
 tar -xf KronaTools-2.7.1.tar && \
 rm KronaTools-2.7.1.tar && \
 cd KronaTools-2.7.1 && \
 ./install.pl --prefix . && \
 ./updateTaxonomy.sh

# Skesa 2.3.0 - DL binary and rename as 'skesa'
RUN mkdir skesa && \
 cd skesa && \
 wget https://github.com/ncbi/SKESA/releases/download/v2.3.0/skesa.centos6.10 && \
 mv skesa.centos6.10 skesa && \
 chmod +x skesa

# Prodigal - 2.6.2 via apt

# Seqtk 1.3 (shovill dep)
RUN wget https://github.com/lh3/seqtk/archive/v1.3.tar.gz && \
 tar -zxf v1.3.tar.gz && \
 rm v1.3.tar.gz && \
 cd seqtk-1.3/ && \
 make && \
 make install

# SPAdes 3.14.0 (needed for shovill and metagenomics workflow)
# apt deps: python
RUN wget http://cab.spbu.ru/files/release3.14.0/SPAdes-3.14.0-Linux.tar.gz && \
 tar -xzf SPAdes-3.14.0-Linux.tar.gz && \
 rm -r SPAdes-3.14.0-Linux.tar.gz

# Mash 2.2 (shovill dep)
RUN wget https://github.com/marbl/Mash/releases/download/v2.2/mash-Linux64-v2.2.tar && \
 tar -xvf mash-Linux64-v2.2.tar && \
 rm -rf mash-Linux64-v2.2.tar

# lighter 1.1.1 (shovill dep)
RUN wget https://github.com/mourisl/Lighter/archive/v1.1.1.tar.gz && \
 tar -zxf v1.1.1.tar.gz && \
 rm -rf v1.1.1.tar.gz && \
 cd Lighter-1.1.1 && \
 make

# trimmomatic 0.38 (shovill dep)
RUN mkdir trimmomatic && \
 cd trimmomatic && \
 wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.38.zip && \
 unzip Trimmomatic-0.38.zip && \
 rm -rf Trimmomatic-0.38.zip && \
 chmod +x Trimmomatic-0.38/trimmomatic-0.38.jar && \
 echo "#!/bin/bash" >> trimmomatic && \
 echo "exec java -jar /trimmomatic/Trimmomatic-0.38/trimmomatic-0.38.jar """"$""@"""" " >> trimmomatic && \
 chmod +x trimmomatic

# BWA 0.7.17 (shovill dep)
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 && \
 tar -xjf bwa-0.7.17.tar.bz2 && \
 rm bwa-0.7.17.tar.bz2 && \
 cd bwa-0.7.17 && \
 make

# Samtools 1.9 (shovill dep)
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
 tar -xjf samtools-1.9.tar.bz2 && \
 rm samtools-1.9.tar.bz2 && \
 cd samtools-1.9 && \
 ./configure && \
 make && \
 make install

# MEGAHIT 1.1.4 (shovill dep)
RUN mkdir megahit && \
 cd megahit && \
 wget https://github.com/voutcn/megahit/releases/download/v1.1.4/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz && \
 tar -xzf megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz && \
 rm megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin.tar.gz

# Velvet 1.2.10 (shovill dep)
RUN wget https://github.com/dzerbino/velvet/archive/v1.2.10.tar.gz && \
 tar -xzf v1.2.10.tar.gz && \
 rm -rf v1.2.10.tar.gz && \
 cd velvet-1.2.10 && \
 make

# Flash 1.2.11 (shovill dep)
RUN wget https://sourceforge.net/projects/flashpage/files/FLASH-1.2.11.tar.gz && \
 tar -zxf FLASH-1.2.11.tar.gz && \
 rm -rf FLASH-1.2.11.tar.gz && \
 cd FLASH-1.2.11 && \
 make

# Pilon 1.22 (shovill dep)
RUN mkdir pilon && \
 cd pilon && \
 wget https://github.com/broadinstitute/pilon/releases/download/v1.22/pilon-1.22.jar && \
 chmod +x pilon-1.22.jar && \
 echo "#!/bin/bash" >> pilon && \
 echo "exec java -jar /pilon/pilon-1.22.jar """"$""@"""" " >> pilon && \
 chmod +x pilon

# Samclip
RUN mkdir samclip && \
 cd samclip && \
 wget https://raw.githubusercontent.com/tseemann/samclip/master/samclip && \
 chmod +x samclip

# Shovill 1.0.4
# apt deps: pigz zlib1g-dev make gcc g++ libpthread-stubs0-dev openjdk-9-jre unzip bzip2 libncurses5-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libfindbin-libs-perl
RUN wget https://github.com/tseemann/shovill/archive/v1.0.4.tar.gz && \
 tar -xzf v1.0.4.tar.gz && \
 rm v1.0.4.tar.gz

# mlst 2.16.2 
# (had to downgrade since later versions of mlst require perl 5.26.0 which is not available on apt for ubuntu:xenial)
# dependencies in apt: libmoo-perl liblist-moreutils-perl libjson-perl gzip file
# other dependencies: any2fasta ncbi-blast+
RUN wget https://github.com/tseemann/mlst/archive/v2.16.2.tar.gz &&\
 tar -xzf v2.16.2.tar.gz &&\
 rm v2.16.2.tar.gz

# any2fasta
RUN cd /usr/local/bin && \
 wget https://raw.githubusercontent.com/tseemann/any2fasta/master/any2fasta && \
 chmod +x any2fasta

# ncbi-blast+ 2.9.0
# blast version in apt for ubuntu:xenial is 2.2.31 (from 2014?)
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
 tar -xzf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
 rm ncbi-blast-2.9.0+-x64-linux.tar.gz

# staramr 0.5.1
# apt deps: python3 python3-pip python3-setuptools git ncbi-blast+ (blast installed manually)
# update pip3 and install staramr 0.5.1
RUN python3 -m pip install -U pip && \
 pip3 install staramr==0.5.1

# ColorID 1.4.3
RUN mkdir colorid && \
 cd colorid && \
 wget https://github.com/hcdenbakker/colorid/releases/download/v0.1.4.3/colorid_Linux64v0.1.4.3 && \
 mv colorid_Linux64v0.1.4.3 colorid && \
 chmod +x colorid

# SalmID 0.1.23
# apt deps: python-setuptools python3 python3-pip curl build-essential file git python3-venv
RUN pip3 install poetry && \
 git clone https://github.com/hcdenbakker/SalmID.git --branch 0.1.23 --single-branch && \
 cd SalmID && \
 poetry build -vvv && \
 pip3 install dist/salmid*.whl

# chewBBACA 2.1.0
# apt deps: prodigal
# BLAST 2.5.0+ or above required
# python deps (installed via pip3 cmd below) numpy scipy biopython plotly SPARQLWrapper
RUN pip3 install chewbbaca==2.1.0

# Get SneakerNet 0.9.8 and make /data
# apt deps: sendmail-base zip bsdmainutils (for column command)
# perl modules listed in cpanm comments above (some installed there, remaining installed w cpanm command below)
ENV SNVER=0.9.8
RUN wget https://github.com/lskatz/SneakerNet/archive/v${SNVER}.tar.gz && \
 tar -zxf v${SNVER}.tar.gz && \
 rm v${SNVER}.tar.gz && \
 cd /SneakerNet-${SNVER} && \
 cpanm --installdeps --notest --force . && \
 perl Makefile.PL && \
 make && \
 sed -i 's+/opt/kraken/full-20140723+/kraken-database/minikraken_20171013_4GB+g' config/settings.conf && \
 mkdir /data

# minikraken db
RUN mkdir /kraken-database && \
 cd /kraken-database && \
 wget  https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz && \
 tar -zxf minikraken_20171019_4GB.tgz && \
 rm -rf minikraken_20171019_4GB.tgz 

# set PATH and perl local settings
ENV PATH="${PATH}:\
/mlst-2.16.2/bin:\
/ncbi-blast-2.9.0+/bin:\
/skesa:\
/opt/kraken:\
/opt/bin:\
/CG-Pipeline/scripts:\
/KronaTools-2.7.1/bin:\
/SPAdes-3.14.0-Linux/bin:\
/mash-Linux64-v2.2:\
/Lighter-1.1.1:\
/trimmomatic:\
/bwa-0.7.17:\
/megahit/megahit_v1.1.4_LINUX_CPUONLY_x86_64-bin:\
/velvet-1.2.10:\
/FLASH-1.2.11:\
/pilon:\
/samclip:\
/shovill-1.0.4/bin:\
/colorid:\
/SneakerNet-${SNVER}/scripts:\
/SneakerNet-${SNVER}/SneakerNet.plugins\
" \
 LC_ALL=C

# check SN dependencies for each workflow
RUN ./SneakerNet-${SNVER}/scripts/SneakerNet.checkdeps.pl default && \
 ./SneakerNet-${SNVER}/scripts/SneakerNet.checkdeps.pl metagenomics && \
 ./SneakerNet-${SNVER}/scripts/SneakerNet.checkdeps.pl cryptosporidium && \
 ./SneakerNet-${SNVER}/scripts/SneakerNet.checkdeps.pl iontorrent

WORKDIR /data
