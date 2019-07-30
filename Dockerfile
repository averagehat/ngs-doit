FROM centos:centos7
MAINTAINER Michael Panciera

ENV PYTHON_VERSION 3.7

RUN yum -y update && \
    yum -y install curl bzip2

ADD . /ngs_doit

WORKDIR /ngs_doit

RUN bash install.sh /ngs_doit/miniconda

ENV PATH=/ngs_doit/miniconda/bin/:$PATH

RUN conda clean --all --yes && \ 
    rm miniconda3.sh && \
    rmp -e --nodeps curl bzip2 && \ 
    yum clean all # this inherited image should `yum clean all` automatically
