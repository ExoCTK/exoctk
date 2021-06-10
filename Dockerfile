FROM ubuntu:latest
SHELL ["/bin/bash", "--login", "-c"]

# Update and install dependencies
RUN apt-get update
RUN apt-get install -y bash
RUN apt-get install -y curl

# Install miniconda3
RUN curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
RUN chmod 700 ./Miniconda3-latest-Linux-x86_64.sh
RUN bash ./Miniconda3-latest-Linux-x86_64.sh -b -p /root/miniconda3
ENV PATH /root/miniconda3/bin:$PATH

# Install initial conda environment
RUN conda create -n exoctk-3.7 python=3.7
RUN source activate exoctk-3.7
RUN conda install git

# Install the exoctk-3.7 conda environment
RUN cd /root/
RUN git clone https://github.com/ExoCTK/exoctk.git
RUN conda env update -f exoctk/env/environment-3.7.yml
RUN python exoctk/setup.py develop

# Test the installation
RUN python -c "import exoctk"

# Get the ExoCTK data package
