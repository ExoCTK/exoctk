# syntax=docker/dockerfile:1

<<<<<<< HEAD
=======
# Determine whether to clone or copy
ARG dev_source

>>>>>>> upstream/main
# Get python version
ARG python_version

# Create Container
<<<<<<< HEAD
FROM python:${python_version}-alpine
=======
FROM python:${python_version}-alpine AS base
>>>>>>> upstream/main

# Get remaining arguments
ARG exoctk_version

# set up working directory
WORKDIR /exoctk

# Set up container environment
ENV EXOCTK_DATA="/exoctk/exoctk_data/"
ENV MODELGRID_DIR="/exoctk/exoctk_data/modelgrid/"
ENV EXOCTK_VERSION=${exoctk_version}
<<<<<<< HEAD
=======
ENV FLASK_APP=app_exoctk.py
ENV FLASK_RUN_HOST="0.0.0.0"
ENV FLASK_RUN_PORT="5000"
ENV SHARED_DATA_DIR="/exoctk/data"
>>>>>>> upstream/main

# Install necessary base software
RUN apk update
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories
RUN apk --no-cache --update-cache add gcc gfortran build-base wget freetype-dev libpng-dev openblas-dev git hdf5-dev
<<<<<<< HEAD

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

# Clone exoctk repository
RUN git clone https://github.com/ExoCTK/exoctk.git
=======
RUN apk --no-cache --update-cache add nano curl firefox geckodriver firefox-esr

WORKDIR /exoctk

# Install redis and celery
RUN pip install gunicorn
RUN pip install 'batman-package@git+https://github.com/lkreidberg/batman.git@3c69c43'
RUN pip install svo-filters==0.5.3
RUN pip install selenium
RUN pip install pytest-selenium
RUN pip install tabulate

ARG exoctk_repo

# Clone the repository (if the environment is set to do that)
FROM base AS clone

# Clone exoctk repository
WORKDIR /exoctk
ADD https://api.github.com/repos/${exoctk_repo}/git/refs/heads/$EXOCTK_VERSION version.json
RUN git clone https://github.com/${exoctk_repo}
>>>>>>> upstream/main
WORKDIR /exoctk/exoctk
RUN git checkout $EXOCTK_VERSION
RUN pip install -e .

<<<<<<< HEAD
# Set up the port that exoctk will listen on
EXPOSE 5000

# Run exoctk
WORKDIR /exoctk
# CMD ["tail", "-f", "/dev/null"]
CMD ["python", "/exoctk/exoctk/exoctk/exoctk_app/app_exoctk.py"]
=======
# Clone the repository from a detached-head release (if set to do that)
FROM base AS release

# Clone exoctk repository
WORKDIR /exoctk
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
COPY ./docker/exoctk/test_exoctk_page.py test_exoctk_page.py
COPY ./docker/exoctk/exoctk_test_script.sh exoctk_test_script.sh
RUN chmod +x exoctk_test_script.sh

# Run exoctk
WORKDIR /exoctk/exoctk/exoctk
CMD ["tail", "-f", "/dev/null"]
>>>>>>> upstream/main
