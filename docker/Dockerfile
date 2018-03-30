#################################################################
# Dockerfile
#
# Software:         MEA
# Software Version: 1.0
# Description:      Methylomic and Epigenomic Allele-specific analysis pipeline
# Website:          https://github.com/julienrichardalbert/MEA
# Provides:         MEA
# Base Image:       biowardrobe2/mea:v0.0.1
# Build Cmd:        docker build --rm -t biowardrobe2/mea:v0.0.1 -f mea-Dockerfile .
# Pull Cmd:         docker pull biowardrobe2/mea:v0.0.1
# Run Cmd:          docker run --rm -ti biowardrobe2/mea:v0.0.1 /bin/bash
#################################################################

### Base Image
FROM biowardrobe2/scidap:v0.0.3
MAINTAINER Michael Kotliar "misha.kotliar@gmail.com"
ENV DEBIAN_FRONTEND noninteractive

################## BEGIN INSTALLATION ######################

WORKDIR /tmp


ENV VERSION_HTSLIB 1.4
ENV URL_HTSLIB "https://github.com/samtools/htslib/releases/download/${VERSION_HTSLIB}/htslib-${VERSION_HTSLIB}.tar.bz2"

ENV VERSION_SAMTOOLS 1.4
ENV URL_SAMTOOLS "https://github.com/samtools/samtools/releases/download/${VERSION_SAMTOOLS}/samtools-${VERSION_SAMTOOLS}.tar.bz2"

ENV VERSION_BWA 0.7.17
ENV URL_BWA "https://github.com/lh3/bwa/archive/v${VERSION_BWA}.tar.gz"

ENV VERSION_BOWTIE 1.2.0
ENV URL_BOWTIE "https://github.com/BenLangmead/bowtie/archive/v${VERSION_BOWTIE}.tar.gz"

ENV VERSION_BOWTIE2 2.3.0
ENV URL_BOWTIE2 "https://github.com/BenLangmead/bowtie2/archive/v${VERSION_BOWTIE2}.tar.gz"

ENV VERSION_UCSC_APPS v358
ENV URL_UCSC_APPS "http://hgdownload.soe.ucsc.edu/admin/exe/userApps.${VERSION_UCSC_APPS}.src.tgz"

# Hardcoded version v2.r900
ENV URL_SHAPEIT "https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r900.glibcv2.17.linux.tar.gz"

ENV VERSION_VCFTOOL 0.1.15
ENV URL_VCFTOOL "https://github.com/vcftools/vcftools/archive/v${VERSION_VCFTOOL}.tar.gz"

ENV VERSION_BEDTOOLS 2.26.0
ENV URL_BEDTOOLS "https://github.com/arq5x/bedtools2/releases/download/v${VERSION_BEDTOOLS}/bedtools-${VERSION_BEDTOOLS}.tar.gz"

ENV VERSION_BEDOPS 2.4.32
ENV URL_BEDOPS "https://github.com/bedops/bedops/archive/v${VERSION_BEDOPS}.tar.gz"

ENV VERSION_BISMARK 0.19.0
ENV URL_BISMARK "https://github.com/FelixKrueger/Bismark/archive/${VERSION_BISMARK}.tar.gz"

ENV VERSION_TOPHAT 2.1.1
ENV URL_TOPHAT "https://github.com/infphilo/tophat/archive/v${VERSION_TOPHAT}.tar.gz"

ENV VERSION_STAR 2.5.3a
ENV URL_STAR "https://github.com/alexdobin/STAR/archive/${VERSION_STAR}.tar.gz"

ENV VERSION_MEA 1.0
ENV URL_MEA "https://github.com/michael-kotliar/mea/archive/v${VERSION_MEA}.tar.gz"


### Updating the system

RUN apt-get clean all &&\
    apt-get update &&\
    apt-get install -y libncurses5-dev libssl-dev libmysqlclient-dev libpng12-dev uuid-dev liblzma-dev libbz2-dev libboost-all-dev libtbb-dev &&\

### Installing samtools
    cd /tmp &&\
    wget -q $URL_HTSLIB &&\
    bzip2 -d htslib-${VERSION_HTSLIB}.tar.bz2 &&\
    tar -xf htslib-${VERSION_HTSLIB}.tar &&\
    cd htslib-${VERSION_HTSLIB} &&\
    ./configure &&\
    make -j 4 &&\
    make install &&\
    cd .. &&\
    rm -rf ./* &&\

    wget -q $URL_SAMTOOLS &&\
    bzip2 -d samtools-${VERSION_SAMTOOLS}.tar.bz2 &&\
    tar -xf samtools-${VERSION_SAMTOOLS}.tar &&\
    cd samtools-${VERSION_SAMTOOLS} &&\
    ./configure &&\
    make -j 4 &&\
    make install &&\
    cd .. &&\
    rm -rf ./samtools-${VERSION_SAMTOOLS} &&\
    rm -rf ./samtools-${VERSION_SAMTOOLS}.tar &&\

### Installing BWA
    cd /tmp &&\
    wget -q -O - $URL_BWA | tar -zxv &&\
    cd bwa-${VERSION_BWA} &&\
    make -j 4 &&\
    cd .. &&\
    cp ./bwa-${VERSION_BWA}/bwa /usr/local/bin/ &&\
    rm -rf ./* &&\

### Installing bowtie
    cd /tmp &&\
    wget -q -O - $URL_BOWTIE | tar -zxv && \
    cd bowtie-${VERSION_BOWTIE} && \
    make -j 4 && \
    cd .. && \
    cp ./bowtie-${VERSION_BOWTIE}/bowtie /usr/local/bin/ && \
    cp ./bowtie-${VERSION_BOWTIE}/bowtie-* /usr/local/bin/ && \
    rm -rf ./* && \

### Installing bowtie2
    cd /tmp &&\
    wget -q -O - $URL_BOWTIE2 | tar -zxv && \
    cd bowtie2-${VERSION_BOWTIE2} && \
    make -j 4 && \
    cd .. && \
    cp ./bowtie2-${VERSION_BOWTIE2}/bowtie2 /usr/local/bin/ && \
    cp ./bowtie2-${VERSION_BOWTIE2}/bowtie2-* /usr/local/bin/ && \
    rm -rf ./* && \

### Installing ucsc userApps
    cd /tmp &&\
    wget -q -O - $URL_UCSC_APPS | tar -zxv --strip-components=2 &&\
    make -j 4 &&\
    cp ./bin/* /usr/local/bin/ &&\
    rm -rf ./* &&\

### Installing SHAPEIT
    cd /tmp &&\
    wget -q -O - $URL_SHAPEIT | tar -zxv --strip-components=1 &&\
    cp ./bin/shapeit /usr/local/bin/ &&\
    chmod a+x /usr/local/bin/shapeit &&\
    rm -rf ./* &&\

### Installing VCFTools
    cd /tmp &&\
    wget -q -O - $URL_VCFTOOL | tar -zxv &&\
    cd vcftools-${VERSION_VCFTOOL} &&\
    ./autogen.sh &&\
    ./configure --prefix=/usr/local &&\
    make -j 4 install &&\
    cd .. &&\
    rm -rf ./* &&\

### Installing Bedtools
    cd /tmp &&\
    wget -q -O - $URL_BEDTOOLS | tar -zxv && \
    cd bedtools2 && \
    make -j 4 && \
    make -j 4 install &&\
    cd .. &&\
    rm -rf ./* &&\

### Installing Bedops
    cd /tmp &&\
    wget -q -O - $URL_BEDOPS | tar -zxv && \
    cd bedops-${VERSION_BEDOPS} && \
    make -j 4 && \
    make -j 4 install &&\
    cp ./bin/* /usr/local/bin/ &&\
    cd .. &&\
    rm -rf ./* &&\

### Installing bismark
    cd /tmp &&\
    wget -q -O /opt/Bismark-${VERSION_BISMARK}.tar.gz ${URL_BISMARK} && \
    cd /opt && \
    tar xzf ./Bismark-${VERSION_BISMARK}.tar.gz &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/NOMe_filtering /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/bam2nuc /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark2bedGraph /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark2report /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark2summary /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark_genome_preparation /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/bismark_methylation_extractor /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/coverage2cytosine /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/deduplicate_bismark /usr/local/bin/ &&\
    ln -s /opt/Bismark-${VERSION_BISMARK}/filter_non_conversion /usr/local/bin/ &&\
    rm /opt/Bismark-${VERSION_BISMARK}.tar.gz &&\
    rm -rf /tmp/* &&\

### Installing tophat2
    cd /tmp &&\
    wget -q -O - $URL_TOPHAT | tar -zxv && \
    cd tophat-${VERSION_TOPHAT} && \
    ./autogen.sh && \
    ./configure && \
    make && \
    make install &&\
    cd .. &&\
    rm -rf ./* &&\

### Installing STAR
    cd /tmp &&\
    wget -q -O - $URL_STAR | tar -zxv && \
    cd STAR-${VERSION_STAR}/source && \
    make -j 4 && \
    cd .. && \
    cp ./bin/Linux_x86_64_static/STAR /usr/local/bin/ && \
    cd .. && \
    rm -rf ./* && \

### Installing MEA
    cd /tmp &&\
    wget -q -O - $URL_MEA | tar -zxv &&\
    cd mea-${VERSION_MEA}/src &&\
    ls | grep -v "mea.config\|mea.setting" | xargs chmod +x && \
    cp ./* /usr/local/bin/ &&\
    cd ../.. &&\
    rm -rf ./* &&\
    strip /usr/local/bin/*; true

