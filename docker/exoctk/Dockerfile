# syntax=docker/dockerfile:1

# Determine whether to clone or copy
ARG dev_source

# Get python version
ARG python_version

# Create Container
FROM python:${python_version}-alpine AS base

# Get remaining arguments
ARG exoctk_version

# set up working directory
WORKDIR /exoctk

# Set up container environment
ENV EXOCTK_DATA="/exoctk/exoctk_data/"
ENV MODELGRID_DIR="/exoctk/exoctk_data/modelgrid/"
ENV EXOCTK_VERSION=${exoctk_version}
ENV FLASK_APP=app_exoctk.py
ENV FLASK_RUN_HOST="0.0.0.0"
ENV FLASK_RUN_PORT="5000"
ENV SHARED_DATA_DIR="/exoctk/data"

# Install necessary base software
RUN apk update
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories
RUN apk --no-cache --update-cache add gcc gfortran build-base wget freetype-dev libpng-dev openblas-dev git hdf5-dev
RUN apk --no-cache --update-cache add nano

# Install numpy with custom batman
WORKDIR /exoctk
RUN pip install --force-reinstall --no-cache-dir numpy==2.2.6
RUN pip uninstall batman-package
RUN git clone https://github.com/lkreidberg/batman.git
WORKDIR /exoctk/batman
RUN pip install .
WORKDIR /exoctk

# Install very specific versions of packages
RUN pip install svo-filters==0.4.4

ARG exoctk_repo

# Clone the repository (if the environment is set to do that)
FROM base AS clone

# Clone exoctk repository
WORKDIR /exoctk
ADD https://api.github.com/repos/${exoctk_repo}/git/refs/heads/$EXOCTK_VERSION version.json
RUN git clone https://github.com/${exoctk_repo}
WORKDIR /exoctk/exoctk
RUN git checkout $EXOCTK_VERSION
RUN pip install -e .

# Copy the local repository (if the environment is set to do that)
FROM base AS copy

# Add local exoctk repository
ADD ./ /exoctk/exoctk/
WORKDIR /exoctk/exoctk
RUN pip install -e .

# Bring it all back together
FROM ${dev_source} AS final

# Set up the port that exoctk will listen on
EXPOSE 5000

# Add the entrypoint scripts
WORKDIR /exoctk
COPY ./docker/exoctk/exoctk_entrypoint.sh exoctk_entrypoint.sh
RUN chmod +x exoctk_entrypoint.sh
COPY ./docker/exoctk/worker_entrypoint.sh worker_entrypoint.sh
RUN chmod +x worker_entrypoint.sh

# Run exoctk
WORKDIR /exoctk/exoctk/exoctk
CMD ["tail", "-f", "/dev/null"]
