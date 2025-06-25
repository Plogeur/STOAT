FROM ubuntu:24.04

MAINTAINER xian-hui.chang@inrae.fr

ARG THREADS=8

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND=noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN=true

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    sudo \
    python3 \
    python3-dev \
    build-essential \
    cmake \
    wget \
    git \
    g++ \
    gcc \
    make \
    doxygen \
    pkg-config \
    libjansson-dev \
    graphviz \
    zlib1g-dev \
    libbz2-dev \
    libncurses-dev \
    libbz2-dev  \
    protobuf-compiler \
    libprotoc-dev \
    libprotobuf-dev \
    automake \
    gettext \
    autopoint \
    libtool \
    jq \
    bsdmainutils \
    bc \
    rs \
    parallel \
    npm \
    curl \
    unzip \
    redland-utils \
    librdf-dev \
    bison \
    flex \
    gawk \
    lzma-dev \
    liblzma-dev \
    liblz4-dev \
    libffi-dev \
    libcairo-dev \
    libboost-all-dev \
    libzstd-dev \
    pybind11-dev \
    python3-pybind11 \
    samtools \
    gtk-doc-tools \
    libgoogle-perftools-dev \
    tabix \
    bcftools \
    pandoc \
    raptor2-utils \
    rasqal-utils \
    python3-pip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /work
# vg
RUN git clone --recursive https://github.com/vgteam/vg.git
WORKDIR /work/vg
RUN make get-deps

RUN make -j $THREADS

ENV PATH=$PATH:/work/vg/bin/

# htslib for libvgio

# Get htslib
WORKDIR /deps
RUN wget https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xvjf htslib-1.22.tar.bz2
WORKDIR /deps/htslib-1.22
RUN ./configure --prefix=/usr/ &&\
    make && \
    make install

# pangwas
WORKDIR /work
RUN git clone --recursive https://forgemia.inra.fr/xian-hui.chang/pangwas.git
RUN mkdir pangwas/build
WORKDIR /work/pangwas/build
RUN cmake .. && \
    make

ENV PATH=$PATH:/work/pangwas/bin/

## python and libbdsg
#WORKDIR /work
#
#RUN git clone --recursive https://github.com/vgteam/libbdsg.git && \
#    mkdir libbdsg/build && \
#    cd libbdsg/build && \
#    cmake .. && \
#    make -j 4
#
#ENV PYTHONPATH $PYTHONPATH:/work/libbdsg/lib
#
