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
    zlib1g-dev \
    libzmq-dev \
    make \
    python3 \
    python3-pip \
    python3-setuptools \
    python3-matplotlib \
    python-virtualenv \
    python-pip \
    git

RUN apt-get autoclean -y && \
    apt-get autoremove -y

# matplotlib
RUN apt-get install -y \
    pkg-config \
    libfreetype6-dev \
    libpng12-dev \
    zlib1g-dev \
    python3-matplotlib


#################
### samtools
RUN apt-get install -y samtools


#################
### bedtools, need at least version 2.26, where merge command reports strand
# RUN apt-get install -y bedtools
WORKDIR /tmp/bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz
RUN tar -zxvf bedtools-2.26.0.tar.gz
WORKDIR /tmp/bedtools/bedtools2
RUN make
RUN make install


#################
### RNA-star
WORKDIR /tmp/STAR
RUN wget https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz
RUN tar -xvzf STAR_2.4.2a.tar.gz
WORKDIR /tmp/STAR/STAR-STAR_2.4.2a/source
RUN make STAR
RUN mkdir -p /home/icuser/bin && cp STAR /home/icuser/bin
WORKDIR /tmp
RUN rm -rfv STAR


#################
#### iCount
RUN chown -R icuser.icuser /home/icuser

USER icuser
WORKDIR /home/icuser
ADD requirements.txt /home/icuser
RUN virtualenv -p python3 /home/icuser/.icountenv
RUN .icountenv/bin/pip install --upgrade -r requirements.txt

USER root
ADD . /home/icuser/iCount_src
RUN chown -R icuser.icuser /home/icuser

USER icuser
WORKDIR /home/icuser/iCount_src
RUN ../.icountenv/bin/python setup.py develop

USER root
RUN echo "source /home/icuser/.icountenv/bin/activate" >> /etc/bash.bashrc
USER icuser

RUN mkdir /home/icuser/storage

EXPOSE 6543

CMD ["/bin/bash"]
