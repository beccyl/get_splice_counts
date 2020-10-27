FROM ubuntu:18.04

MAINTAINER Rebecca Louise Evans <rebecca.louise.evans@gmail.com>

LABEL description="RetroHunter" version="0.4"

#RUN find /var/lib/apt/lists -type f -exec rm -rf '{}' \;

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
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        libz-dev \
        make \
        wget \
        xz-utils \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local

ENV SAMTOOLS_VERSION 1.9

# Set Standard settings
ENV PATH /usr/bin:/bin:/usr/local/bin

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
    && cd samtools-${SAMTOOLS_VERSION}/ \
    && ./configure && make && make install

ENV PATH "/usr/local/samtools-${SAMTOOLS_VERSION}/:${PATH}"

# install htslib
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar -xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
    && rm -rf htslib-${SAMTOOLS_VERSION}.tar.bz2  \
    && cd htslib-${SAMTOOLS_VERSION}/ \
    && ./configure && make && make install

ENV PATH "/usr/local/htslib-${SAMTOOLS_VERSION}/:${PATH}"

# link libraries
RUN /sbin/ldconfig

# environment variables needed by Makefile
ENV CPPFLAGS "-I/usr/local/samtools-${SAMTOOLS_VERSION} -I/usr/local/htslib-${SAMTOOLS_VERSION}"
ENV LDFLAGS "-L/usr/local/samtools-${SAMTOOLS_VERSION} -L/usr/local/htslib-${SAMTOOLS_VERSION}"

#RUN git clone https://github.com/nadiadavidson/get_splice_counts.git

# fix the BamReader header file
COPY BamReader.h get_splice_counts/BamReader.h
COPY ReferenceReader.h get_splice_counts/ReferenceReader.h
COPY retro_hunter.c++ get_splice_counts/retro_hunter.c++

# fix the Makefile
COPY Makefile.docker get_splice_counts/Makefile

RUN cd get_splice_counts && make

RUN chmod 755 get_splice_counts/retro_hunter

ENV PATH "/usr/local/get_splice_counts/:${PATH}"

# set the working directory
WORKDIR /data

# copy reference file
COPY GRCh38_gencodev22_with_id.exon_flank.fasta GRCh38_gencodev22_with_id.exon_flank.fasta

#CMD ["/bin/bash"]
CMD ["retro_hunter", "/data/GRCh38_gencodev22.exon_flank.fasta", "/data/batch_run/ERR11223344/ERR11223344", "/data/batch_run/ERR11223344/ERR11223344.bam"] 
