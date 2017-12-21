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

Install DLO-HiC-Tools
---------------------

Install from source code
^^^^^^^^^^^^^^^^^^^^^^^^
Download the source code firsyly, then install with setup.py ::

    $ git clone https://github.com/Nanguage/DLO-HiC-Tools.git
    $ cd DLO-HiC-Tools/
    $ python3 setup.py install
