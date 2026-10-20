# Recipe of the published Singularity image (biorsurfaces.sif).
# The May 2024 image contains Ubuntu 22.04.3, Intel oneAPI 2024.0 (MKL 2024.0,
# TBB 2021.11), Clang 14 and CMake 3.22. The original recipe used the
# "latest" tag; the tag below is the oneAPI 2024.0 release for Ubuntu 22.04.
#
# Build from the repository root:
#   docker build -f singularity/Dockerfile -t biorsurfaces .
#   docker save biorsurfaces > biorsurfaces.tar
#   singularity build biorsurfaces.sif docker-archive://biorsurfaces.tar
FROM intel/oneapi-hpckit:2024.0.1-devel-ubuntu22.04

# Set environment variables
ENV CC=clang
ENV CXX=clang++

# Install dependencies
RUN apt-get update && apt-get install -y \
    git \
    xorg-dev\
    libglu1-mesa-dev\
    freeglut3-dev\
    mesa-common-dev\
    libomp-dev\
    cmake \
    build-essential \
    clang\
    libx11-dev \
    libnetcdf-dev \
    libxrandr-dev \
    libxinerama-dev \
    libnetcdf-c++4-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the source tree (this repository)
COPY . /app/biorsurfaces/

RUN cd /app/biorsurfaces && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j$(nproc)
# Copy the executable into the container
RUN mkdir /app/simulation &&\
    cp /app/biorsurfaces/build/bin/biorsurfaces /app/simulation

WORKDIR /app/simulation
