FROM ubuntu:16.04

LABEL name="dlohic"
LABEL base.image="ubuntu:16.04"
LABEL version="1"
LABEL software="DLO-HiC-Tools"
LABEL software.version="0.1.0"
LABEL website="https://github.com/Nanguage/DLO-HiC-Tools"
LABEL documentation="https://github.com/Nanguage/DLO-HiC-Tools/blob/master/doc/source/install.rst"
LABEL license="https://github.com/Nanguage/DLO-HiC-Tools/blob/master/LICENSE"
LABEL tags="Bioinformatics,Genomics,Hi-C,DLO-Hi-C"

MAINTAINER nanguage@yahoo.com

ENV MINICONDA https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
ENV DLOHIC https://github.com/Nanguage/DLO-HiC-Tools/archive/v0.1.1.tar.gz
ENV JUICERTOOLS http://hicfiles.tc4ga.com.s3.amazonaws.com/public/juicer/juicer_tools.1.8.9_jcuda.0.8.jar

ENV LD_ALL C.UTF-8
ENV LANG C.UTF-8

# add apt source and install requirement packages
RUN mv /etc/apt/sources.list /etc/apt/sources.list.bkp && \
    bash -c 'echo -e "deb mirror://mirrors.ubuntu.com/mirrors.txt xenial main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-updates main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-backports main restricted universe multiverse\n\
deb mirror://mirrors.ubuntu.com/mirrors.txt xenial-security main restricted universe multiverse\n\n" > /etc/apt/sources.list' && \
    cat /etc/apt/sources.list.bkp >> /etc/apt/sources.list && \
    cat /etc/apt/sources.list
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y  \
        autotools-dev   \
        automake        \
        cmake           \
        curl            \
        grep            \
        sed             \
        dpkg            \
        fuse            \
        git             \
        wget            \
        zip             \
        openjdk-8-jre   \
        build-essential \
        pkg-config      \
        python          \
	python-dev      \
        python-pip      \
        bzip2           \
        ca-certificates \
        libglib2.0-0    \
        libxext6        \
        libsm6          \
        libxrender1     \
        git             \
        mercurial       \
        subversion      \
        zlib1g-dev &&   \
        apt-get clean && \
        apt-get purge && \
        rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# download juicertools
RUN wget ${JUICERTOOLS}

# install miniconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet ${MINICONDA} -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

ENV PATH=/opt/conda/bin:$PATH

# add channels
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda upgrade conda

# install dependency
RUN conda install --yes cython numpy scipy matplotlib pandas h5py jsonschema
RUN conda install --yes -c coreutils bioconda bwa samtools bedtools pairix tabix pysam

# install DLO-HiC-Tools
RUN wget ${DLOHIC} -O DLO-HiC-Tools.tar.gz && \
    tar zxvf DLO-HiC-Tools.tar.gz && \
    mkdir DLO-HiC-Tools && \
    tar zxvf DLO-HiC-Tools.tar.gz -C DLO-HiC-Tools --strip-components=1 && \
    cd DLO-HiC-Tools && \
    python setup.py install


CMD ["/bin/bash"]

WORKDIR /data
