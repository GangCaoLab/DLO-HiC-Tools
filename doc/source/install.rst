Installation
============

Requirements
------------
This package is only support UNIX-like system, and Python version must be >= 3.4, 
and following softwares are required:

- BWA (>=0.7.15)
- samtools (>=1.6)
- coreutils (>=8.25)
- pairix
- tabix
- cooler
- mafft
- java (optional, for create .hic file)

Recommend install and manage requirements with `Anaconda <https://conda.io/miniconda.html>`_, just use following commands::

    $ conda config --add channels defaults  # add channels
    $ conda config --add channels bioconda
    $ conda config --add channels conda-forge
    $ conda create -n dlohic python=3.6  # create virtual environment
    $ source activate dlohic
    (dlohic) $ conda install --yes numpy scipy matplotlib pandas cython h5py jsonschema graphviz
    (dlohic) $ conda install -c bioconda --yes coreutils bwa samtools mafft pairix tabix cooler pysam java-jdk
    (dlohic) $ conda install -c conda-forge --yes bzip2

Install DLO-HiC-Tools
---------------------

Install from source code
^^^^^^^^^^^^^^^^^^^^^^^^
Download the source code firsyly, then install with setup.py ::

    $ (dlohic) git clone https://github.com/Nanguage/DLO-HiC-Tools.git
    $ (dlohic) cd DLO-HiC-Tools/
    $ (dlohic) python setup.py install

Using Docker
------------

You can also using the docker, pull the dlohic image: ::

    $ docker pull nanguage/dlohic


Run a container, with mount current directory in file system to the '/data' in the container: ::

    $ docker run -v $PWD:/data -ti nanguage/dlohic:latest

