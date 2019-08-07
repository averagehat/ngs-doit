FROM centos:centos7
MAINTAINER Michael Panciera

ENV PYTHON_VERSION 3.7

RUN yum -y update && \
    yum -y install curl bzip2 git   # git to allow git+ pip installs, bzip2 for conda

ADD . /src

WORKDIR /src 

RUN  bash install.sh /miniconda

ENV PATH=/miniconda/bin/:$PATH

RUN conda clean --all --yes && \ 
    rm miniconda3.sh && \
    rpm -e --nodeps curl bzip2 && \ 
    yum clean all # this inherited image should `yum clean all` automatically
