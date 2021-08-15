ENV_NAME := biophi-dev
SHELL := /bin/bash
CONDA_ACTIVATE = eval "$$(conda shell.bash hook)" && conda activate $(ENV_NAME)
DOCKER_CELERY_CPUS := 4

env: env-create env-setup

env-create: environment.yml
	conda env create -n $(ENV_NAME) -f $<

env-update: environment.yml
	conda env update -n $(ENV_NAME) -f $<

env-setup:
	$(CONDA_ACTIVATE); pip install -e .

pytest:
	$(CONDA_ACTIVATE); \
        biophi -h; \
        pytest tests/pytest

docker-build:
	docker-compose build

docker-run:
	docker-compose up

docker-python:
	docker-compose run worker python

web:
	$(CONDA_ACTIVATE); \
		FLASK_APP=biophi.common.web.views \
		FLASK_ENV=development \
		flask run

celery:
	$(CONDA_ACTIVATE); \
		FLASK_ENV=development \
		celery -A biophi.common.web.tasks.celery worker --loglevel=INFO
