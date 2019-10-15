"""
DLO HiC command line interface.
"""
import matplotlib
matplotlib.use("Agg")

import sys
import os
from os.path import join, abspath, dirname
import logging
import importlib

import click

from .config import LOGGING_FMT, LOGGING_DATE_FMT
import dlo_hic

log = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stderr)
handler.setFormatter(logging.Formatter(fmt=LOGGING_FMT, datefmt=LOGGING_DATE_FMT))
log.addHandler(handler)
log.setLevel(logging.DEBUG)


CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])

@click.group(name="dlohic", context_settings=CONTEXT_SETTINGS)
@click.option("--log-file", help="The log file, default output to stdout.")
@click.option("--debug", is_flag=True,
    help="If debug will print all information, default False.")
@click.version_option(version=dlo_hic.__version__)
def cli(log_file, debug):
    """ DLO HiC Tools command line interface. """
    log = logging.getLogger() # root logger
    if log_file:
        handler = logging.FileHandler(log_file)
    else:
        handler = logging.StreamHandler(sys.stderr)
    fomatter = logging.Formatter(
        fmt=LOGGING_FMT,
        datefmt=LOGGING_DATE_FMT)
    handler.setFormatter(fomatter)
    log.addHandler(handler)

    if debug:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)


def search_tools():
    """ get all tools's module name
    like:
        dlo_hic.tools.main_process.build_bedpe
        dlo_hic.tools.main_process.bedpe2pairs
        dlo_hic.tools.helper.extract_rest_sites
        ...
    """
    tools_dir = join(dirname(abspath(__file__)), "tools") + os.path.sep
    tools = []
    for subdir, dirs, files in os.walk(tools_dir):
        for file in files:
            filepath = join(subdir, file)
            if filepath.endswith('.py') and "__init__" not in filepath:
                tool = filepath.replace(tools_dir, "")
                tool = tool.replace(".py", "")
                tool = tool.replace(os.path.sep, '.')
                tool = "dlo_hic.tools." + tool
                tools.append(tool)
    return tools


# add command to cli
for tool in search_tools():
    try:
        module = importlib.import_module(tool)
        command = module._main
        cli.add_command(command)
    except ImportError as e:
        log.error(str(e))
    except AttributeError:
        log.debug("%s command line interface not implemented."%tool)
