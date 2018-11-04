Installation
============

Requirements
------------
This package is only support UNIX-like system, and Python version must be >= 3.4, 
and following softwares are required:

- BWA
- samtools
- pairix
- tabix

Recommend install them with `Anaconda <https://conda.io/miniconda.html>`_, just use following commands::

    $ conda -c bioconda install bwa samtools pairix tabix

Install Python packages
^^^^^^^^^^^^^^^^^^^^^^^
Python packages requirements check and install process is integrated in setup.py,
but should install install Cython moudlue firstly.::

    $ pip install cython

And suggest install and update numpy, scipy and matplotlib with conda: ::

    $ conda install numpy scipy matplotlib

Install DLO-HiC-Tools
---------------------

Install from source code
^^^^^^^^^^^^^^^^^^^^^^^^
Download the source code firsyly, then install with setup.py ::

    $ git clone https://github.com/Nanguage/DLO-HiC-Tools.git
    $ cd DLO-HiC-Tools/
    $ python3 setup.py install

Using Docker
------------

You can also using the docker, pull the dlohic image: ::

    $ docker pull nanguage/dlohic


Run a container, with mount current directory in file system to the '/data' in the container: ::

    $ docker run -v $PWD:/data -ti nanguage/dlohic:latest

