FROM ubuntu:20.04

LABEL description="methylink"
LABEL maintainer="Mei Wu"

ENV WORK_DIR=/work

# Define versions
ARG HTSLIB_VERSION=1.18

# Install essential packages
RUN apt-get update && \
    apt-get install -y \
        cmake \
        build-essential \
        wget && \
    apt-get clean

# Install htslib
WORKDIR /opt/htslib
RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar -xvf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    rm -r htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure && \
    make


WORKDIR /opt/methylink
COPY ./* .

RUN python -m pip install .