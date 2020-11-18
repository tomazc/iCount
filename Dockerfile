FROM continuumio/anaconda3
MAINTAINER Tomaz Curk <tomazc@gmail.com>

# suppress prompt for tzdata
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Ljubljana

# thanks to https://github.com/bschiffthaler/ngs/blob/master/base/Dockerfile
# and https://github.com/AveraSD/ngs-docker-star/blob/master/Dockerfile

SHELL ["/bin/bash", "-c"]

RUN conda update -n base -c defaults conda -y
RUN conda install -c conda-forge mamba -y

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

### samtools
#RUN conda install -c bioconda -y "samtools>=1.10"

### bedtools, need at least version 2.26, where merge command reports strand
# RUN conda install -c bioconda -y bedtools

### RNA-star
# RUN conda install -c bioconda -y star

#################
#### iCount
#################

RUN useradd -m -d /home/icuser icuser
RUN chown -R icuser.icuser /home/icuser

USER icuser
WORKDIR /home/icuser
RUN mkdir /home/icuser/storage
RUN git clone https://github.com/tomazc/iCount.git --branch snakemake

RUN conda create -c conda-forge -c bioconda -n iCount_pipeline3 -y
RUN conda init bash
RUN echo "conda activate iCount_pipeline3" >> ~/.bashrc
#RUN conda env update --file iCount/conda_iCount.yaml # needs ~ 4 GB RAM, otherwise killed
#RUN pip install ./iCount

#ENV PATH /home/icuser/bin:$PATH
#WORKDIR /home/icuser
#CMD ["/bin/bash"]
