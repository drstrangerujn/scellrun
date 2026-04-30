# Minimal scellrun container — single-stage, conda-managed Python 3.11.
#
# Build:
#   docker build -t scellrun .
#
# Run end-to-end pipeline on an h5ad in the current directory:
#   docker run --rm -v "$PWD":/work scellrun analyze data.h5ad --profile joint-disease --no-ai
#
# This image is intentionally bare-bones: no GPU, no jupyter, no Anthropic SDK
# (install with `pip install scellrun[ai]` if you want LLM-augmented stages).
# v1.0 ships the build recipe; tagged-image push to a registry is a v1.1 task.
FROM continuumio/miniconda3:latest

RUN conda create -n scellrun python=3.11 -y && \
    echo "source activate scellrun" > ~/.bashrc
ENV PATH=/opt/conda/envs/scellrun/bin:$PATH

RUN pip install --no-cache-dir scellrun

WORKDIR /work
ENTRYPOINT ["scellrun"]
CMD ["--help"]
