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

# ensure /data is present
ARG ENOG_DESCRIPTIONS_FILENAME
RUN mount | grep '/data' || exit 1
RUN [ -d "/data/models" ] || exit 1
RUN [ -f "/data/${ENOG_DESCRIPTIONS_FILENAME}" ] || exit 1

# install environment
RUN apt-get update --fix-missing \
    && apt-get install -y apache2 apache2-dev mariadb-server libmariadbclient-dev git hmmer
RUN conda env update -n base -f conda.yaml && conda clean -a -y
RUN git clone https://github.com/univieCUBE/phenotrex.git ../phenotrex && pip install ../phenotrex
# SHAP version currently broken: phenotrex requires 0.35 but the models require 0.37
RUN pip install shap==0.37.0
# load environment variables - prefer ones already present (pass to docker build)
RUN source source/maintenance_scripts/variables.sh

# set up database
RUN service start mysql
RUN mysql < devel_scripts/set_up_dev.sql
RUN python3 source/web_server/manage.py makemigrations phenotypePredictionApp
RUN python3 source/web_server/manage.py migrate

# get variables
RUN source source/maintenance_scripts/variables.sh

# add taxonomy data
RUN mkdir -p /apps/miniconda3/opt/krona/taxonomy
RUN python3 source/maintenance_scripts/add_taxonomy_to_db.py /apps/miniconda3/opt/krona/taxonomy

# add Enog data
RUN python3 source/maintenance_scripts/add_enogs_to_db.py ${ENOG_ANNOTATION_FILE}

# add Models
RUN python3 source/maintenance_scripts/add_models_to_db.py ${PHENDB_MODEL_DIR}


RUN ln -s source/maintenance_scripts/phenDB.sh ./phenDB.sh
RUN ln -s /data ./data

