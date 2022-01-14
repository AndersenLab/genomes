# docker build --file gatk4.Dockerfile -t swantonlab/manta .
FROM continuumio/miniconda3
RUN apt-get update && apt-get install -y procps && \
    apt-get clean
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge
RUN conda create -n genome \
                        bioconda::bwa=0.7.17 \
                        bioconda::gatk4=4.1.7.0 \
                        bioconda::samtools=1.10 \
                        bioconda::picard=2.22.4 \
                        xsv=0.13.0 \
                        snpeff=4.3.1t \
                        parallel=20200322 \
                        r-tidyverse \
                        r-data.table \
                        r-vroom \
                        bedtools \ 
                        bedops \
                        ucsc-bigbedtobed
    && conda clean -a
ENV PATH /opt/conda/envs/genome/bin:$PATH
RUN conda env export --name genome > genome.yml

LABEL Name="genome" Author="Daniel Cook"