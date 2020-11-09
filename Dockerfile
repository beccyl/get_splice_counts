FROM ubuntu:18.04

MAINTAINER Rebecca Louise Evans <rebecca.louise.evans@gmail.com>

LABEL description="RetroHunter" version="0.4"

USER root

RUN apt-get clean \
    && apt-get update -qq \
    && apt-get install -qq \
        bzip2 \
        gcc \
        g++ \
        git \
        gzip \
        libbz2-dev \
        libcurl4-openssl-dev \
        libcrypto++-dev \
        libgsl0-dev \
        liblzma-dev \
        libncurses5-dev \
        libperl-dev \
        libssl-dev \
        libz-dev \
        make \
        perl \
        wget \
        xz-utils \
        zlib1g-dev \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt

ENV HTS_VERSION 1.11

# Set Standard settings
ENV PATH /usr/bin:/bin:/usr/local/bin

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/${HTS_VERSION}/htslib-${HTS_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTS_VERSION}.tar.bz2 \
    && rm htslib-${HTS_VERSION}.tar.bz2  \
    && cd htslib-${HTS_VERSION}/ \
    && ./configure && make install \
    && /sbin/ldconfig

#RUN git clone https://github.com/nadiadavidson/get_splice_counts.git
RUN git clone https://github.com/beccyl/get_splice_counts.git

# use the Makefile for docker
RUN mv get_splice_counts/Makefile.docker get_splice_counts/Makefile

RUN make -C get_splice_counts

RUN chmod 755 get_splice_counts/retro_hunter

ENV PATH "${PATH}:/opt/get_splice_counts"

# set the working directory
WORKDIR /data

# copy reference file
COPY GRCh38_gencodev22_with_id.exon_flank.fasta GRCh38_gencodev22_with_id.exon_flank.fasta

#CMD ["/bin/bash"]
CMD ["retro_hunter", "/data/GRCh38_gencodev22.exon_flank.fasta", "/data/batch_run/ERR11223344/ERR11223344", "/data/batch_run/ERR11223344/ERR11223344.bam"] 
