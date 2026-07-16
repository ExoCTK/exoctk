# syntax=docker/dockerfile:1

# Get python version and create final container
ARG python_version
FROM python:${python_version}-alpine

ENV FORTGRID_DIR=/data/exoctk/fortney/fortney_models.db
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

# Install batman
RUN pip install 'batman-package@git+https://github.com/lkreidberg/batman.git@3c69c43'

# Install pandexo
ARG pandexo_branch
# This endpoint accepts either a branch name or an immutable commit SHA and
# changes whenever the resolved commit changes, invalidating Docker's cache.
ADD https://api.github.com/repos/natashabatalha/PandExo/commits/${pandexo_branch} /tmp/pandexo_ref.json
RUN pip install --no-cache-dir "pandexo.engine[all]@git+https://github.com/natashabatalha/PandExo.git@${pandexo_branch}"

# Optionally replace PandExo's declared Pandeia dependency with a local
# release-candidate checkout.  Compose supplies an inert fallback context
# when PANDEIA_SOURCE is unset.
RUN --mount=type=bind,from=pandeia_rc,target=/tmp/pandeia \
    if [ -f /tmp/pandeia/pyproject.toml ] || [ -f /tmp/pandeia/setup.py ]; then \
        cp -a /tmp/pandeia /tmp/pandeia-build && \
        pip install --no-cache-dir /tmp/pandeia-build && \
        rm -rf /tmp/pandeia-build; \
    fi

# Add the entrypoint scripts
COPY ./pandexo_entrypoint.sh pandexo_entrypoint.sh
RUN chmod +x pandexo_entrypoint.sh

# Run Pandexo
CMD ["tail", "-f", "/dev/null"]
