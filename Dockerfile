FROM continuumio/miniconda3:4.8.2

LABEL maintainer="Lukas Lüftinger <lukas.lueftinger@ares-genetics.com>"
LABEL description="A container for running the phenDB web server."
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TINI_VERSION v0.19.0
SHELL ["/bin/bash", "-c"]

ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT [ "/usr/bin/tini", "--" ]

RUN mkdir -p /apps
ADD . /apps/phenDB
WORKDIR /apps/phenDB

# install environment
RUN apt-get update --fix-missing \
    && apt-get install -y apache2 apache2-dev mariadb-server libmariadbclient-dev git hmmer sudo
RUN conda env update -n base -f conda.yaml && conda clean -a -y
RUN git clone https://github.com/univieCUBE/phenotrex.git ../phenotrex && pip install ../phenotrex
# SHAP version currently broken: phenotrex requires 0.35 but the models require 0.37
RUN pip install shap==0.37.0
# load environment variables - prefer ones already present (pass to docker build)
RUN source source/maintenance_scripts/variables.sh

# set up database
RUN bash devel_scripts/inital_setup.sh
RUN ln -s source/maintenance_scripts/phenDB.sh ./phenDB.sh

