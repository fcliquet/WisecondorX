FROM mambaorg/micromamba:2.0-noble

LABEL maintainer="fcliquet"
LABEL org.opencontainers.image.source="https://github.com/fcliquet/WisecondorX"

COPY --chown=$MAMBA_USER:$MAMBA_USER conda.yml /tmp/conda.yml
RUN micromamba install -y -n base -f /tmp/conda.yml \
    && micromamba clean -a -y

USER root
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

COPY . /app
WORKDIR /app
RUN pip install --no-cache-dir .

USER $MAMBA_USER
WORKDIR /home/$MAMBA_USER

ENTRYPOINT ["wisecondorx"]
