FROM continuumio/miniconda3:4.9.2
LABEL authors="Gisela Gabernet <gisela.gabernet@gmail.com>,Edmund Miller <edmund.a.miller@gmail.com>" \
    description="Docker image containing Space Ranger"
# Disclaimer: this container is not provided nor supported by 10x Genomics.

# Install procps and curl and clean apt cache
RUN apt-get update --allow-releaseinfo-change \
    && apt-get install -y procps curl \
    && apt-get clean -y && rm -rf /var/lib/apt/lists/*

# Copy pre-downloaded spaceranger file
ENV SPACERANGER_VER 3.0.0
COPY spaceranger-$SPACERANGER_VER.tar.gz /opt/spaceranger-$SPACERANGER_VER.tar.gz

# Install spaceranger
RUN \
    cd /opt && \
    tar -xzvf spaceranger-$SPACERANGER_VER.tar.gz && \
    export PATH=/opt/spaceranger-$SPACERANGER_VER:$PATH && \
    ln -s /opt/spaceranger-$SPACERANGER_VER/spaceranger /usr/bin/spaceranger && \
    rm -rf /opt/spaceranger-$SPACERANGER_VER.tar
