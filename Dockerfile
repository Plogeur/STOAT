FROM python:3.12

# Prevent dpkg from trying to ask any questions, ever
ENV DEBIAN_FRONTEND noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN true

## python, snakemake and awscli
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    git screen wget curl gcc less nano \
    sudo \
    python3-dev \
    python3-pip \
    python3-setuptools \
    python3-wheel \
    make \
    g++ \
    pkg-config \
    build-essential \ 
    cmake \
    doxygen \
    && rm -rf /var/lib/apt/lists/*

RUN pip3 install --upgrade pip
RUN pip3 install --no-cache-dir pybind11

ENV TZ=America/Los_Angeles
WORKDIR /build

# jansson
RUN wget --no-check-certificate https://github.com/akheron/jansson/releases/download/v2.14/jansson-2.14.tar.bz2 && \
tar -xjf jansson-2.14.tar.bz2

RUN cd jansson-2.14 && ./configure && make  -j 6 && make install && \
cd .. && rm -rf jansson-2.14 jansson-2.14.tar.bz2

## libbdsg
RUN git clone --recursive https://github.com/vgteam/libbdsg.git && \
    cd libbdsg && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 6

RUN pip3 install setuptools_scm
RUN pip3 install libbdsg/.

# Clone the STOAT repository and set it as the working directory
RUN git clone https://github.com/Plogeur/STOAT && \
    cd STOAT && \
    pip3 install . && \  
    cd .. && rm -rf STOAT

# # Clone the STOAT-NODE repository and set it as the working directory
# RUN git clone https://github.com/Plogeur/STOAT --branch stoat2bin && \
#     cd STOAT && \ 
#     pip3 install . && \
#     cd .. && rm -rf STOAT

# init home
WORKDIR /home

# docker build . -t stoat