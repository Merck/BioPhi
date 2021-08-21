FROM continuumio/miniconda3

RUN apt-get update -y && apt-get install -y --no-install-recommends \
        build-essential

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
