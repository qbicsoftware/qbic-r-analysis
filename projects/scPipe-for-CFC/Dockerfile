FROM continuumio/miniconda

LABEL name="scPipe-for-CFC-r-analysis"
LABEL authors="gurpreet-bioinfo"
LABEL organization="Quantitative Biology Center"
LABEL description="scPipe-for-CFC"
LABEL maintainer="gurpreet-bioinfo <sven.fillinger@qbic.uni-tuebingen.de>"  
LABEL version="0.1"
LABEL github="https://github.com/qbicsoftware"

# Install procps so that Nextflow can poll CPU usage
RUN apt-get update && apt-get install -y procps && apt-get clean -y 
# Update the base version of conda
RUN conda update -n base conda
# Copy the list of R packages into the container's /tmp dir
COPY environment.yml /tmp
RUN conda update -n base conda && \
    conda env create -f /tmp/environment.yml && \
    conda clean -a

ENV PATH /opt/conda/envs/scPipe-for-CFC-ranalysis-0.1/bin:$PATH
