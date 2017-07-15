FROM ubuntu:14.04.3
MAINTAINER Tomaz Curk <tomazc@gmail.com>

# thanks to https://github.com/bschiffthaler/ngs/blob/master/base/Dockerfile
# and https://github.com/AveraSD/ngs-docker-star/blob/master/Dockerfile

RUN useradd -m -d /home/icuser icuser

# update system
RUN apt-get update && apt-get upgrade -y && \
    apt-get install -y \
    build-essential \
    gfortran \
    libatlas-base-dev \
    wget \
    g++ \
    make \
    python3 \
    python3-pip \
    python3-setuptools \
    python-virtualenv \
    python-pip \
    git && \
    apt-get build-dep -y python3-matplotlib && \
    apt-get install -y \
    python3-matplotlib

RUN apt-get autoclean -y && \
    apt-get autoremove -y


#################
### samtools
RUN apt-get install -y \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    samtools


#################
#### iCount
RUN chown -R icuser.icuser /home/icuser

USER icuser
WORKDIR /home/icuser
RUN virtualenv -p python3 /home/icuser/.icountenv

USER root
RUN echo "source /home/icuser/.icountenv/bin/activate" >> /etc/bash.bashrc
USER icuser

WORKDIR /home/icuser

#################
### prerequisites for iCount
RUN /home/icuser/.icountenv/bin/pip install numpy
RUN /home/icuser/.icountenv/bin/pip install pandas
RUN /home/icuser/.icountenv/bin/pip install cutadapt
RUN /home/icuser/.icountenv/bin/pip install pysam
RUN /home/icuser/.icountenv/bin/pip install pybedtools
RUN /home/icuser/.icountenv/bin/pip install numpydoc
RUN /home/icuser/.icountenv/bin/pip install sphinx
RUN /home/icuser/.icountenv/bin/pip install matplotlib

CMD ["/bin/bash"]
