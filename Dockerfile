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
    libeigen3-dev \
    libboost-all-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libsdsl-dev \
    libjansson-dev \
    && rm -rf /var/lib/apt/lists/*

# Installation de BDSG depuis les sources
WORKDIR /tmp
RUN git clone --recursive https://github.com/vgteam/libbdsg.git \
    && cd libbdsg \
    && mkdir build \
    && cd build \
    && cmake -DRUN_DOXYGEN=OFF .. \
    && make -j$(nproc) \
    && make install \
    && cd / \
    && rm -rf /tmp/libbdsg

# Installation de Catch2
WORKDIR /tmp
RUN git clone https://github.com/catchorg/Catch2.git \
    && cd Catch2 \
    && cmake -Bbuild -H. -DBUILD_TESTING=OFF \
    && cmake --build build/ \
    && cmake --install build/ \
    && cd / \
    && rm -rf /tmp/Catch2

# Clone the STOAT C++ repository and set it as the working directory
RUN git clone https://github.com/Plogeur/STOAT --branch stoat_cxx \
    && cd STOAT \ 
    && mkdir build \
    && cd build \
    && cmake .. \
    && make -j 4

# Définir le répertoire de travail
WORKDIR /home
