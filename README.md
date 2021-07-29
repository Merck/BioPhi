<p align="center">
    <br>
    <img height="90" src="biophi/common/web/static/img/logo-light/2x/biophi_logo@2x.png?raw=true?raw=true">
    <br>
</p>

BioPhi is an open-source antibody design platform. 
It features methods for automated antibody humanization (Sapiens), humanness evaluation (OASis) and an interface for computer-assisted antibody sequence design.

BioPhi is available at: [http://biophi.dichlab.org](http://biophi.dichlab.org)

Learn more in the BioPhi, Sapiens and OASis in our pre-print:

> TBD

## Running BioPhi on your machine

If you don't want to use the [public BioPhi server](http://biophi.dichlab.org), you can run BioPhi on your own machine.

<details>
<summary>See more</summary>

### 1. Download OASis database

To run BioPhi with OASis humanness evaluation locally, 
you will need to download the OASis database (22GB uncompressed).

### 2a. Run simplified server using Docker

If you have [Docker](https://www.docker.com/products/docker-desktop), 
you can run a simplified BioPhi server using:
```bash
docker run TBD
```

### 2b. Run simplified server using Conda

You can also install BioPhi using [Bioconda](https://bioconda.github.io/user/install.html):

```bash
# Recommended: Create a separate BioPhi environment
conda create -n biophi python=3.8
conda activate biophi

# Install BioPhi
conda install biophi

# Run simplified BioPhi server (not for live deployment!)
biophi web
```

**Note:** This is simplified usage for local use only. See Deploying your own BioPhi server *(TODO LINK)* to learn about 
deploying BioPhi properly on a server.

</details>

## BioPhi command-line interface

BioPhi also provides a command-line interface that enables bulk processing.

<details>
    <summary>See more</summary>

```bash
TBD
```
  
</details>

## Contributing

BioPhi is an open and extensible platform, contributions are welcome. 
Submit any feature requests on the [Issues](https://github.com/Merck/biophi/issues) page.

## Development

BioPhi is composed of three services that need to be running at the same time:

- `web`: Flask web server that handles both the frontend and the backend of the web application
- `celery`: Asynchronous worker service(s) that process long-running tasks
- `redis`: In-memory database for storing celery queue tasks and results

### Run BioPhi dev server through Docker Compose

Running through Docker Compose is easiest in terms of setup, but web server autoreload is not supported,
so you will have to restart the services after each code update.

<details>
    <summary>See more</summary>

#### 1. Install Docker

See https://docs.docker.com/docker-for-mac/install/

#### 2. Build all images using Docker Compose

```bash
# Run using Makefile
make docker-build
# or directly using
docker-compose build
```

#### 3. Run all services using Docker Compose

```bash
# Run using Makefile
make docker-run
# or directly using
docker-compose up
```

To build and run, you can use:
```bash
# Run using Makefile
make docker-build docker-run
# or directly using
docker-compose up --build
```

#### 4. Handle code updates

After your code is updated, you will need to stop the services, run build and start again. 
See the next section for info on running locally with flask auto-reload.

</details>


### Run BioPhi dev server using Conda

Running each service locally using Conda will enable flask auto-reload, 
which is useful if you are going back and forth between your IDE and the browser.

<details>
    <summary>See more</summary>

#### 1. Install Conda

Install [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/download.html) 
or one of the alternatives ([Miniconda](https://docs.conda.io/en/latest/miniconda.html), 
[Miniforge](https://github.com/conda-forge/miniforge))

#### 2. Install Redis server

Install and run [Redis server](https://redis.io/download). 
On Mac, you can [install Redis using Brew](https://medium.com/@petehouston/install-and-config-redis-on-mac-os-x-via-homebrew-eb8df9a4f298).

#### 3. Setup environment

```bash
# Install dependencies using the provided Makefile
make env
# Or directly using
conda env create -n biophi-dev -f environment.yml
conda activate biophi-dev
pip install -e .
```

#### 4. Run all services

You will have to run each service in a separate terminal (Use Cmd+T to open a new tab):

```bash
# Run Redis server (this depends on your installation, the server might already be running)
redis-server

# In a separate terminal, run celery worker queue
make celery

# In a separate terminal, run flask web server
make web
```

See the provided 

#### 5. Handle code updates

After your code is updated, the flask web service should refresh automatically. 
However, the celery service needs to be stopped and started manually, 
so you will need to do that if you update code that is executed from the workers.
</details>

## Deploying your own BioPhi server

You can deploy your own internal BioPhi server. 
You will need to run the three separate services - the flask web server, 
the celery worker and the redis database.

This will depend on your platform and your cloud provider, the easiest deployment is using [Docker Compose](https://docs.docker.com/compose/gettingstarted/)
through the provided [docker-compose.yml](docker-compose.yml) file.

## Acknowledgements

BioPhi is based on antibody repertoires from the Observed Antibody Space:

> Kovaltsuk, A., Leem, J., Kelm, S., Snowden, J., Deane, C. M., & Krawczyk, K. (2018). Observed Antibody Space: A Resource for Data Mining Next-Generation Sequencing of Antibody Repertoires. The Journal of Immunology, 201(8), 2502–2509. https://doi.org/10.4049/jimmunol.1800708

Antibody numbering is performed using ANARCI:

> Dunbar, J., & Deane, C. M. (2016). ANARCI: Antigen receptor numbering and receptor classification. Bioinformatics, 32(2), 298–300. https://doi.org/10.1093/bioinformatics/btv552