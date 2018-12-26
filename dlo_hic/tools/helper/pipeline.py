import os
from os.path import join
import logging

import click


log = logging.getLogger(__name__)

HERE = os.path.dirname(os.path.abspath(__file__))

EXAMPLE_CONFIG = join(HERE, "../../templates/pipeline_config.ini")
SNAKEFILE = join(HERE, "../../templates/pipeline.snake")

def gen_file(template, filename, tag):
    if not os.path.exists(filename):
        log.info("Generate %s at %s"%(tag, filename))
        with open(template) as fi, open(filename, 'w') as fo:
            fo.write(fi.read())
    else:
        log.info("%s aleardy exist."%filename)

def gen_snakefile():
    snake = "./Snakefile"
    tag = "pipeline (Snakemake file)"
    gen_file(SNAKEFILE, snake, tag)

def gen_config():
    example_config = "./pipeline_config.ini"
    tag = "config file"
    gen_file(EXAMPLE_CONFIG, example_config, tag)



@click.command(name="pipeline")
def _main():
    """
    Generate integrated main processes pipeline(Snakemake file) and example config file.

    After pipeline and config generated,
    edit config with editor, then run pipeline with Snakemake.

    About Snakemake: http://snakemake.readthedocs.io/en/latest/

    example:

    $ dlo_hic pipeline
    $ ls
    pipeline-config.ini Snakefile
    $ vi pipeline-config.ini # edit pipeline config file.
    $ snakemake all # run pipeline

    """
    gen_config()
    gen_snakefile()


if __name__ == "__main__":
    eval("_main()")