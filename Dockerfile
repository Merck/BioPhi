FROM continuumio/miniconda

RUN apt-get update -y && apt-get install -y --no-install-recommends \
        build-essential

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge

WORKDIR /opt/biophi

COPY environment.yml .
COPY Makefile .

RUN make env-update ENV_NAME=base

COPY . .

RUN make env-setup ENV_NAME=base

RUN useradd docker \
  && mkdir /home/docker \
  && chown docker:docker /home/docker \
  && addgroup docker staff
USER docker

CMD [ "biophi", "web", "--host", "0.0.0.0" ]
