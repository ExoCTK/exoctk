# syntax=docker/dockerfile:1

# Get python version and create final container
ARG python_version
FROM python:${python_version}-alpine

ENV FORTGRID_DIR=/data/exoctk/fortney
ENV pandeia_refdata=/data/pandeia
ENV PYSYN_CDBS=/data/trds

# Get remaining arguments
ARG pandexo_version
ARG pandeia_version

# set up working directory
WORKDIR /pandexo

# Set up container environment
ENV PANDEXO_VERSION=${pandexo_version}
ENV PANDEIA_VERSION=${pandeia_version}

# Install necessary base software
RUN apk update
RUN echo "http://dl-8.alpinelinux.org/alpine/edge/community" >> /etc/apk/repositories
RUN apk --no-cache --update-cache add gcc gfortran build-base wget freetype-dev libpng-dev openblas-dev hdf5-dev git
RUN apk --no-cache --update-cache add linux-headers

# Install numpy with custom batman
RUN pip install --force-reinstall --no-cache-dir numpy==2.2.6
RUN pip uninstall batman-package
RUN git clone https://github.com/lkreidberg/batman
WORKDIR /pandexo/batman
RUN pip install .
WORKDIR /pandexo

# Install pandexo
RUN pip install numpy==1.26.4
RUN pip install scipy==1.14.1
# ***** NOTE *****
# In the future, it would be useful to be able to either directly install this
# with pip (i.e. pip install pandexo==VERSION), or do a git checkout of a given
# version (i.e. git checkout VERSION)
RUN git clone --recursive https://github.com/natashabatalha/PandExo.git
WORKDIR /pandexo/PandExo
RUN pip install .
RUN pip install pandeia.engine==${pandeia_version}

WORKDIR /pandexo
CMD ["start_pandexo", "--workers=10"]
