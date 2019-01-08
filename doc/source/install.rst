Installation
============

Requirements
------------
This package is only support UNIX-like system, and Python version must be >= 3.4, 
and following softwares are required:

- BWA
- samtools
- bedtools
- pairix
- tabix
- cooler

Recommend install and manage requirements with `Anaconda <https://conda.io/miniconda.html>`_, just use following commands::

    $ conda create -n dlohic python=3.6  # create virtual environment
    $ source activate dlohic
    (dlohic) $ conda install numpy scipy matplotlib pandas cython h5py jsonschema
    (dlohic) $ conda install -c bioconda --yes bwa samtools bedtools pairix tabix cooler pysam

Install DLO-HiC-Tools
---------------------

Install from source code
^^^^^^^^^^^^^^^^^^^^^^^^
Download the source code firsyly, then install with setup.py ::

    $ git clone https://github.com/Nanguage/DLO-HiC-Tools.git
    $ cd DLO-HiC-Tools/
    $ python setup.py install

Using Docker
------------

You can also using the docker, pull the dlohic image: ::

    $ docker pull nanguage/dlohic


Run a container, with mount current directory in file system to the '/data' in the container: ::

    $ docker run -v $PWD:/data -ti nanguage/dlohic:latest

