FROM continuumio/anaconda3
MAINTAINER Tomaz Curk <tomazc@gmail.com>

# suppress prompt for tzdata
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Ljubljana

# thanks to https://github.com/bschiffthaler/ngs/blob/master/base/Dockerfile
# and https://github.com/AveraSD/ngs-docker-star/blob/master/Dockerfile

RUN conda update -n base -c defaults conda -y

RUN conda update -n base -c defaults conda
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

### samtools
RUN conda install -c bioconda -y "samtools>=1.10"

### bedtools, need at least version 2.26, where merge command reports strand
RUN conda install -c bioconda -y bedtools

### RNA-star
RUN conda install -c bioconda -y star

#################
#### iCount
#################

#RUN useradd -m -d /home/icuser icuser

#RUN chown -R icuser.icuser /home/icuser

#USER icuser
#WORKDIR /home/icuser
#RUN virtualenv -p python3 /home/icuser/.icountenv

#USER root
# to speed-up building of Docker images
#RUN /home/icuser/.icountenv/bin/pip install numpy pandas pysam pybedtools numpydoc matplotlib

#ADD . /home/icuser/iCount_src
#RUN chown -R icuser.icuser /home/icuser

#USER icuser
#WORKDIR /home/icuser/iCount_src

#RUN ../.icountenv/bin/pip install -e .[docs,test]

#USER root
#RUN echo "source /home/icuser/.icountenv/bin/activate" >> /etc/bash.bashrc
#USER icuser

#RUN mkdir /home/icuser/storage

#ENV PATH /home/icuser/bin:$PATH

#WORKDIR /home/icuser

#CMD ["/bin/bash"]
