FROM continuumio/miniconda3

COPY environment.yml .
RUN conda env create -f environment.yml