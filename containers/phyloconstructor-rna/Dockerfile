FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="PhyloConstructor RNA tools"
LABEL org.opencontainers.image.description="SRA Toolkit, pigz, BBMap and SPAdes/rnaSPAdes for PhyloConstructor"
LABEL org.opencontainers.image.version="1.0.0"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        --strict-channel-priority \
        sra-tools=3.4.1 \
        bbmap=39.00 \
        spades=4.2.0 \
    && micromamba clean --all --yes

# Writable SRA configuration directory.
ENV NCBI_SETTINGS=/tmp/.ncbi/user-settings.mkfg

# Tests performed during construction.
RUN fasterq-dump --version \
    && prefetch --version \
    && gzip -h \
    && bbduk.sh --version \
    && rnaspades.py --version

WORKDIR /data

CMD ["bash"]