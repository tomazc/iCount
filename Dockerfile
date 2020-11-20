FROM continuumio/anaconda3:2020.07

MAINTAINER Tomaz Curk <tomazc@gmail.com>

# suppress prompt for tzdata
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Europe/Ljubljana

# thanks to https://github.com/bschiffthaler/ngs/blob/master/base/Dockerfile
# and https://github.com/AveraSD/ngs-docker-star/blob/master/Dockerfile

RUN apt-get install -y vim
RUN apt-get install -y nano

SHELL ["/bin/bash", "--login", "-c"]

RUN conda update -n base -c defaults conda -y
RUN conda install -c conda-forge mamba -y

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

#################
#### iCount
#################

RUN useradd -m -d /home/icuser icuser
ADD . /home/icuser/iCount
RUN chown -R icuser.icuser /home/icuser

USER icuser
WORKDIR /home/icuser
RUN mkdir /home/icuser/storage

RUN conda create -c conda-forge -c bioconda -n iCount_pipeline3 -y
RUN conda init bash
RUN echo "conda activate iCount_pipeline3" >> ~/.bashrc

SHELL ["conda", "run", "-n", "iCount_pipeline3", "/bin/bash", "-c"]

### needs ~ 4 GB RAM, otherwise killed
RUN conda env update --file iCount/conda_iCount.yaml

RUN pip install ./iCount

ENV PATH /home/icuser/bin:$PATH
WORKDIR /home/icuser
CMD ["/bin/bash"]
