FROM mambaorg/micromamba:1.5.5

USER root

RUN apt-get update && \
    apt-get install -y procps zip curl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER ./container/environment.yml /tmp/

COPY --chown=$MAMBA_USER:$MAMBA_USER ./scripts /scripts

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

COPY --chown=$MAMBA_USER:$MAMBA_USER ./scripts/vep.py /opt/conda/lib/python3.13/site-packages/multiqc/modules/vep/vep.py

RUN curl -sSL -o /scripts/clinvar_tab_delimited_index.html https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/ && \
    curl -sSL -o /scripts/variant_summary.txt.gz https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

ENV MULTIQC_CONFIG_PATH="/scripts/multiqc_config.yaml"

ENV PATH="/opt/conda/bin:$PATH"