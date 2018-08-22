Quick Start
===========

DLO-HiC-Tools is an integrated automatic pipeline for
DLO-HiC data analyze. Before we start, you should install package 
and some dependency. see `<install.rst>`__ page.

DLO-HiC-Tools implemented some tools for DLO-HiC data analyze and quality control.
You can list them and get more usage detail by `dlohic` command:

.. code-block::

    $ dlohic --help
    Usage: __main__.py [OPTIONS] COMMAND [ARGS]...

        DLO HiC Tools command line interface.

    Options:
      --log-file TEXT       The log file, default output to stdout.
      --debug / --no-debug  If debug will print all information, default True.
      --help                Show this message and exit.

    Commands:
      PET_span_dist         Count the distribution of PET span.
      bedpe2pairs           Transform bedpe format file to pairs format,...
      build_bedpe           Build bedpe file from fastq or sam/bam file.
      extract_PET           Extract the PETs sequences on both sides of...
      extract_rest_sites    Extract all restriction sites from fasta...
      fragment_length_dist  Draw the distribution figure(kde/box plot),...
      interactions_qc       Count ratio of: inter-chromosome...
      noise_reduce          Remove DLO-HiC noise (self-ligation).
      pipeline              Generate integrated main processes...
      remove_redundancy     Remove the redundancy within pairs.

Use the pipeline
----------------

DLO-HiC-Tools implemented the pipeline with `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_.
If you want generate the HiC matrix at once, you can use this pipeline.
But one thing you should noticed is that you should remove the sequencing adapters firstly,
see data_preprocessing_.

The pipeline contain two files, a Snakemake file and a configuration file.

Firstly, generate them by `dlohic pipeline` command. It will generate the necessary files
under your working directory.

.. code-block::

    $ dlohic pipeline
    dlo_hic.tools.helper.pipeline INFO    @ 08/21/18 21:40:30: Generate config file at ./pipeline_config.ini
    dlo_hic.tools.helper.pipeline INFO    @ 08/21/18 21:40:30: Generate pipeline (Snakemake file) at ./Snakefile
    $ ls
    pipeline_config.ini  Snakefile

Then edit the `pipeline_config.ini` file, just fill the sections with it's comment information.
And you can reference the example data and the configuration file, see `here <https://github.com/Nanguage/DLO-HiC-Tools/tree/master/example/cli>`_.

After all necessary information wrote in to the configuration file, You can run the pipeline.

.. code-block::

    $ snakemake -j 16

The `-j` parameter indicate the number of jobs, it depend on your cpu cores number.
Or, if you use the cluster like pbs system, you can run like this:

.. code-block::

    $ snakemake --cluster qsub -j 16


In addition, you can visualize the pipeline(need `Graphviz <https://www.graphviz.org/>`_ installed):

.. code-block::

    $ snakemake all --dag | dot -Tpng > pipeline.png

.. image:: ../img/pipeline.png

More information about the Snakemake, please read it's document.


.. _data_preprocessing:

Data preprocessing
^^^^^^^^^^^^^^^^^^

Before run the pipeline, you should remove the adapters from the raw fastq file.
You can use the `Trimmomatic <http://www.usadellab.org/cms/?page=trimmomatic>`_ do this, for example:

.. code-block::

    $ java -jar $trimmo_jar SE test.fq test.clean.fq ILLUMINACLIP:$adapter_dir/TruSeq3-SE.fa:3:30:10:1:TRUE


Visualize the result:
---------------------

The pipeline can generate two kind of result matrix format, `.hic` and `.cool`.
You can use some software to visualize them, for example:

- `Juicebox <https://github.com/theaidenlab/Juicebox>`_ provided the GUI for `.hic` file.
- `CoolBox <https://github.com/Nanguage/CoolBox>`_ for `.cool` file visualization in the Jupyter nootebook.
