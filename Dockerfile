FROM continuumio/miniconda3:4.8.2

LABEL maintainer="Lukas Lüftinger <lukas.lueftinger@ares-genetics.com>"
LABEL description="A container for running the phenDB web server."
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TINI_VERSION v0.19.0

ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT [ "/usr/bin/tini", "--" ]

RUN mkdir -p /apps
ADD . /apps/phenDB
WORKDIR /apps/phenDB

RUN apt-get update --fix-missing \
    && apt-get install -y apache2 apache2-dev mariadb-server libmariadbclient-dev git hmmer
RUN conda env update -n base -f conda.yaml && conda clean -a -y
RUN git clone https://github.com/univieCUBE/phenotrex.git ../phenotrex && pip install ../phenotrex
RUN service start mysql
RUN mysql < devel_scripts/set_up_dev.sql
RUN source source/maintenance_scripts/variables.sh
