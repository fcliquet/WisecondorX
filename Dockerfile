FROM continuumio/miniconda3:latest

LABEL maintainer="fcliquet"
LABEL org.opencontainers.image.source="https://github.com/fcliquet/WisecondorX"

# Install R + Bioconductor packages + faiss via conda
RUN conda install -y -c conda-forge -c bioconda \
        python=3.12 \
        r-base \
        r-jsonlite \
        bioconductor-dnacopy \
        pysam \
        faiss-cpu \
    && conda clean -a -y

# Install WisecondorX via pip
COPY . /app
WORKDIR /app
RUN pip install --no-cache-dir .

# Smoke test: verify key tools are available
RUN wisecondorx --help \
    && python -c "import wisecondorx; print('WisecondorX installed')" \
    && python -c "import faiss; print(f'faiss {faiss.__version__}')" \
    && python -c "import pysam; print(f'pysam {pysam.__version__}')" \
    && Rscript -e "library(DNAcopy); cat('DNAcopy OK\n')" \
    && Rscript -e "library(jsonlite); cat('jsonlite OK\n')"

WORKDIR /data
ENTRYPOINT ["wisecondorx"]
