#################################################################
# Dockerfile
#
# Description:      Docker container with BWA-0.7.15 | samtools-1.3 | sambamba-0.6.3 | varscan-2.4.2 to call variants from deep sequencing (CAPP-Seq) data.
# Base Image:       ubuntu
# Pull Cmd:         docker pull anu9109/capp-seq
# Run Cmd:          docker run -it anu9109/capp-seq
# Run tools as:     bwa | samtools | sambamba | 'java -jar /opt/software/varscan/VarScan.v2.4.2.jar'
#################################################################

# Set the base image to Ubuntu
FROM ubuntu

################## BEGIN INSTALLATION ###########################

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
  apt-get install -y \
  git \
  autoconf \
  curl \
  gcc \
  make \
  gawk \
  g++ \
  perl \
  pkg-config \ 
  zlib1g-dev \
  wget \
  libncurses5-dev \
  libcurl4-gnutls-dev \
  libgnutls-dev \
  libssl-dev \
  libexpat1-dev \
  libgd-gd2-perl \
  cpanminus \
  build-essential \
  libgd-dev \
  nettle-dev \
  bioperl \
  default-jre \ 
  default-jdk && \
  apt-get clean && \
  apt-get purge && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*


# create /opt/software
RUN mkdir -p /opt/software

# install BWA-0.7.15
RUN cd /opt/software/ && \
  wget https://github.com/lh3/bwa/releases/download/v0.7.15/bwa-0.7.15.tar.bz2 && \
  tar -xvjf /opt/software/bwa-0.7.15.tar.bz2 && \
  cd /opt/software/bwa-0.7.15 && \
  make && \
  rm /opt/software/bwa-0.7.15.tar.bz2

# install samtools-1.3
RUN cd /opt/software/ && \
  wget https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 && \
  tar -xvjf /opt/software/samtools-1.3.tar.bz2 && \
  cd /opt/software/samtools-1.3 && \
  make && \
  make install && \
  rm /opt/software/samtools-1.3.tar.bz2

# install sambamba-0.6.3
RUN cd /opt/software/ && \
  wget https://github.com/lomereiter/sambamba/releases/download/v0.6.3/sambamba_v0.6.3_linux.tar.bz2 && \
  tar -xvjf sambamba_v0.6.3_linux.tar.bz2 && \
  mv /opt/software/sambamba_v0.6.3 /opt/software/sambamba && \
  rm /opt/software/sambamba_v0.6.3_linux.tar.bz2

# install varscan-2.4.2 
RUN cd /opt/software/ && \
  git clone https://github.com/dkoboldt/varscan.git 

ENV PATH=$PATH:/opt/software:/opt/software/varscan


##################### INSTALLATION END ##########################

# File Author / Maintainer
MAINTAINER Anu Amallraja <anu9109@gmail.com>
