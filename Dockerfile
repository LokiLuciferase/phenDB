FROM continuumio/miniconda3:4.8.2

SHELL ["/bin/bash", "-c"]
LABEL maintainer="Lukas Lüftinger <lukas.lueftinger@ares-genetics.com>"
LABEL description="A container for running the phenDB web server."
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV TINI_VERSION v0.19.0
EXPOSE 80
EXPOSE 3306

ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini
ENTRYPOINT [ "/usr/bin/tini", "--" ]

RUN mkdir -p /apps
ADD . /apps/phenDB
WORKDIR /apps/phenDB

# install environment
RUN apt-get update --fix-missing \
    && apt-get install -y build-essential mariadb-server libmariadbclient-dev git hmmer sudo wget

# if data directory is not given, download from CUBE fileshare
RUN if [ ! -d /apps/phenDB/data ]; then \
    wget https://fileshare.csb.univie.ac.at/phenDB/phenDB_data/latest.tar.gz \
    && tar -xvf latest.tar.gz && rm latest.tar.gz && mv latest/ data/; fi

RUN conda env update -n base -f conda.yaml && conda clean -a -y
RUN git clone https://github.com/univieCUBE/phenotrex.git ../phenotrex && pip install ../phenotrex

# set up database and web server
RUN source source/maintenance_scripts/variables.sh \
    && bash devel_scripts/initial_setup.sh \
    && ln -s source/maintenance_scripts/phenDB.sh ./phenDB.sh

CMD ["source", "source/maintenance_scripts/variables.sh", "&&", "/bin/bash", "phenDB.sh"]
