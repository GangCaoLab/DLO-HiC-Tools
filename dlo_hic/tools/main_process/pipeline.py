import os
from os.path import join
import logging

import click


log = logging.getLogger(__name__)

HERE = os.path.dirname(os.path.abspath(__file__))

EXAMPLE = join(HERE, "../../templates/pipeline_config.ini")
SNAKEFILE = join(HERE, "../../templates/pipeline.snake")


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
    example_config = "./pipeline_config.ini"

    log.info("Generate example config file at %s."%example_config)
    with open(EXAMPLE) as fi, open(example_config, "w") as fo:
        fo.write(fi.read())

    snake = "./Snakefile"

    log.info("Generate pipeline (Snakemake file) at %s"%snake)
    with open(SNAKEFILE) as fi, open(snake, 'w') as fo:
        fo.write(fi.read())
