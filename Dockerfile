FROM continuumio/miniconda3:4.8.2

SHELL ["/bin/bash", "-c"]
LABEL maintainer="Lukas Lüftinger <lukas.lueftinger@ares-genetics.com>"
LABEL description="A container for running the phenDB web server."
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN mkdir -p /apps
ADD . /apps/phenDB
RUN mkdir -p /apps/phenDB/data/  # to use external data folder, mount here
WORKDIR /apps/phenDB

# install environment
RUN apt-get update --fix-missing \
    && apt-get install -y build-essential mariadb-server libmariadbclient-dev git hmmer sudo wget zip unzip gzip

RUN conda env update -n base -f conda.yaml && conda clean -a -y
RUN git clone https://github.com/univieCUBE/phenotrex.git ../phenotrex && pip install ../phenotrex

# set up database and web server
RUN source source/maintenance_scripts/variables.sh \
    && echo "port = 33060" >> /etc/mysql/my.cnf \
    && bash devel_scripts/initial_setup.sh \
    && chmod +x source/pipeline/bin/* \
    && chmod +x source/maintenance_scripts/* \
    && chmod +x source/web_server/manage.py

CMD ["source/maintenance_scripts/run.sh"]
