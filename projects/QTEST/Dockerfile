FROM continuumio/miniconda

LABEL authors="Sven Fillinger"
LABEL organization="Quantitative Biology Center"
LABEL description="Put a short analyses description here"
LABEL maintainer="Sven Fillinger <sven.fillinger@qbic.uni-tuebingen.de>"  
LABEL version="0.1dev"
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

ENV PATH /opt/conda/envs/qbicsoftware-QABCD-ranalyses-0.1dev/bin:$PATH
