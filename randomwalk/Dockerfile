# Dockerfile VERSION = v0.1
# docker build -t random_walk_workflow$VERSION . 

FROM continuumio/miniconda3

RUN conda create --name multixrank python=3.10 -y
RUN python3 -m pip install multixrank
RUN python3 -m pip install scipy==1.8

WORKDIR /random_walk/

#to COPY the remote file at working directory in container
COPY multiXrank_covid_full_4-layers_PCC_0.1_threshold.py .

ENTRYPOINT [ "python3", "multiXrank_covid_full_4-layers_PCC_0.1_threshold.py"]