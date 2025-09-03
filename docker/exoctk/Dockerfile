# syntax=docker/dockerfile:1

# Get python version
ARG python_version

# Create Container
FROM python:${python_version}-alpine

# Get remaining arguments
ARG exoctk_version

# set up working directory
WORKDIR /exoctk

# Set up container environment
ENV EXOCTK_DATA="/exoctk/exoctk_data/"
ENV MODELGRID_DIR="/exoctk/exoctk_data/modelgrid/"
ENV EXOCTK_VERSION=${exoctk_version}

# Install necessary base software
RUN apk update
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories
RUN apk --no-cache --update-cache add gcc gfortran build-base wget freetype-dev libpng-dev openblas-dev git hdf5-dev

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
WORKDIR /exoctk/exoctk
RUN git checkout $EXOCTK_VERSION
RUN pip install -e .

# Set up the port that exoctk will listen on
EXPOSE 5000

# Run exoctk
WORKDIR /exoctk
# CMD ["tail", "-f", "/dev/null"]
CMD ["python", "/exoctk/exoctk/exoctk/exoctk_app/app_exoctk.py"]
