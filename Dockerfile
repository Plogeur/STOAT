FROM ubuntu:22.04

# Éviter les interactions pendant l'installation
ENV DEBIAN_FRONTEND=noninteractive

# Mettre à jour le système et installer les dépendances de base
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    pkg-config \
    wget \
    libhts-dev \
    libboost-all-dev \
    libjansson-dev \
    protobuf-compiler \
    libprotoc-dev \
    libprotobuf-dev \
    valgrind \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /home
# Clone the STOAT C++ repository and set it as the working directory
RUN git clone --recursive https://github.com/Plogeur/STOAT \
    && cd STOAT \ 
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j$(nproc)
