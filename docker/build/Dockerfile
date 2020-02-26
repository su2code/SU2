FROM ubuntu:18.04
ENV LANG C.UTF-8
RUN apt-get update && apt-get install -y \
    python3 \
    pkg-config \
    python3-pip \
    git \
    build-essential \
    python3-numpy \
    python3-scipy \
    python3-mpi4py \
    swig \
    libopenmpi-dev \
    openmpi-bin \
    ccache \
 && rm -rf /var/lib/apt/lists/* \
 && update-alternatives --install /usr/bin/python python /usr/bin/python3 10 \
 && /usr/sbin/update-ccache-symlinks \
 && echo 'export PATH="/usr/lib/ccache:$PATH"' | tee -a ~/.bashrc 

# Copies your code file from your action repository to the filesystem path `/` of the container
COPY compileSU2.sh /compileSU2.sh

# Code file to execute when the docker container starts up (`entrypoint.sh`)
ENTRYPOINT ["/compileSU2.sh"]
