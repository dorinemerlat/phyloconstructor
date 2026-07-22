FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="PhyloConstructor sequence tools"
LABEL org.opencontainers.image.description="PhyKIT, SeqKit, CD-HIT and AGAT for PhyloConstructor"
LABEL org.opencontainers.image.version="1.0.0"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install --yes --name base \
        --channel conda-forge \
        --channel bioconda \
        --strict-channel-priority \
        phykit=2.0.1 \
        seqkit=2.9.0 \
        cd-hit=4.8.1 \
        agat=1.6.1 \
    && micromamba clean --all --yes

# Verification of the commands used by the workflow.
RUN phykit create_concat --help >/dev/null \
    && seqkit version \
    && command -v cd-hit \
    && command -v agat_sp_extract_sequences.pl

WORKDIR /data

CMD ["bash"]