FROM mambaorg/micromamba:1.5.5

USER root

RUN apt-get update && \
    apt-get install -y procps zip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER ./container/environment.yml /tmp/

COPY --chown=$MAMBA_USER:$MAMBA_USER ./scripts /scripts

RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1

ENV MULTIQC_CONFIG_PATH="/scripts/multiqc_config.yaml"

ENV PATH="/opt/conda/bin:$PATH"