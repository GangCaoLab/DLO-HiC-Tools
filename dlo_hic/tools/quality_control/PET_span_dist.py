import io
import gzip
import logging

import click

from dlo_hic.utils.parse_text import Bedpe, Pairs


log = logging.getLogger(__name__)


def log_statistic_info(infos, log_file=None):
    """
    log statistic infos to sperate file or stdout.
    """
    if not log_file:
        for info in infos:
            log.info(info.name + ":\n" + str(info))
    else:
        if log_file == '-':
            f = sys.stdout
        else:
            f = open(log_file, 'a')

        for info in infos:
            f.write(info.name + ":\n")
            f.write(info.to_csv(sep="\t"))

        if log_file != '-':
            f.close()


def open_file(file_name):
    if file_name.endswith(".gz"):
        f = io.TextIOWrapper(gzip.open(file_name))
    else:
        f = open(file_name)
    return f


@click.command(name="PET_span_dist")
@click.argument("input")
@click.option("--log-file",
    default="PET_span_dist.txt",
    help="Sperate to store quantiles information. default: PET_span_dist.txt")
@click.option("--file-format",
    default="pairs",
    type=click.Choice(["pairs", "bedpe"]),
    help="The file format of interactions data. " +\
         "support gziped file. default pairs format")
@click.option("--figure", '-f',
    help="The name of output figure.")
@click.option("--fig-type",
    type=click.Choice(['box', 'kde']),
    default="box",
    help="Figure type, box plot or kde density plot.")
@click.option("--dpi", default=600,
    help="dpi of output figure")
def _main(input, log_file, file_format, figure, fig_type, dpi):
    """
    Count the distribution of PET span.
    Output 
    """
    fmt = Pairs if file_format == 'pairs' else Bedpe

    with open_file(input) as f:
        pass


main = _main.callback


if __name__ == "__main__":
    _main()
