import os
from os.path import join
import logging
import re

import click

from dlo_hic.utils.parse_text import is_comment
from dlo_hic.utils.pipeline import format_ini_config
from dlo_hic import __version__


log = logging.getLogger(__name__)

HERE = os.path.dirname(os.path.abspath(__file__))

EXAMPLE_CONFIG = join(HERE, "../../templates/pipeline_config.ini")
SNAKEFILE = join(HERE, "../../templates/pipeline.snake")


def read_version(path):
    with open(path) as f:
        line = f.readline().strip()
    m = re.match("# pipeline version: (.*)", line)
    if m is None:
        return
    version = m.groups()[0]
    return version


def copy_file(source, target):
    with open(source) as fi, open(target, 'w') as fo:
        fo.write("# pipeline version: {}\n\n".format(__version__))
        fo.write(fi.read())


def backup_path(prefix, suffix='.bak'):
    i = 0
    f = prefix + suffix + ".{}".format(i)
    while os.path.exists(f):
        i += 1
        f = prefix + suffix + ".{}".format(i)
    return f


def gen_file(gen_func, target, tag):
    if not os.path.exists(target):
        log.info("Generate %s at '%s'"%(tag, target))
        gen_func()
    else:
        log.info("'%s' aleardy exist."%target)
        old_ver = read_version(target)
        if old_ver != __version__:
            backup = backup_path(target)
            log.info("'{}' is in old version, backup to '{}'".format(target, backup))
            os.system("mv {} {}".format(target, backup))
            gen_file(gen_func, target, tag)


def gen_snakefile():
    target = "./Snakefile"
    tag = "pipeline (Snakemake file)"
    def gen():
        copy_file(SNAKEFILE, target)
    gen_file(gen, target, tag)


def gen_config():
    target = "./pipeline_config.ini"
    tag = "config file"
    def gen():
        with open(EXAMPLE_CONFIG) as f:
            template = f.read()
        from dlo_hic.config import DEFAULT_PIPELINE_CONFIG
        contents = format_ini_config(template, DEFAULT_PIPELINE_CONFIG)
        with open(target, 'w') as f:
            f.write("# pipeline version: {}\n\n".format(__version__))
            f.write(contents)
    gen_file(gen, target, tag)


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