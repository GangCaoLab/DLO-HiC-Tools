import io
import re
import sys
import gzip
import random
import logging

import click
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from dlo_hic.utils.parse_text import Bedpe, Pairs, is_comment


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
@click.option("--box-plot",
    default="PET_span.box.png",
    help="The file name of output box-plot figure.")
@click.option("--kde-plot",
    default="PET_span.kde.png",
    help="The file name of output kde-plot figure.")
@click.option("--dpi", default=600,
    help="dpi of output figure")
@click.option("--sample", "-s",
    type=int,
    help="Sample the input file, if specified.")
@click.option("--seed", default=1,
    help="The seed for random value generation. default 1")
def _main(input, log_file,  box_plot, kde_plot, dpi, sample, seed):
    """
    Count the distribution of PET span.

    Input:
        bedpe file or pairs file, support gziped file.

    Output:
        1. quantiles information
        2. box plot
        3. kde plot
        4. PET span text file (optional)
    """

    # conform input file format
    if input.endswith("pairs") or input.endswith("pairs.gz"):
        fmt = Pairs
        log.info("input pairs file.")
    elif input.endswith("bedpe") or input.endswith("bedpe.gz"):
        fmt = Bedpe
        log.info("input bedpe file.")
    else:
        raise NotImplementedError("Only support pairs and bedpe file format.")

    if sample:
        random.seed(seed)
        num_lines = sum(1 for line in open(input))
        rowskeep = set(random.sample(range(num_lines), sample))

    spans = []

    with open_file(input) as f:
        for i, line in enumerate(f):
            if sample:
                if i not in rowskeep:
                    continue
            if is_comment(line):
                continue
            record = fmt(line)
            span = abs(record.pos1 - record.pos2)
            spans.append(span)

    ser = pd.Series(spans)
    ser.name = re.sub(".bedpe.gz$|.bedpe$|.pairs.gz$|.pairs", "", input)

    # draw box plot
    ser.plot.box()
    plt.savefig(box_plot, dpi=dpi)

    # draw kde plot
    ser.plot.kde()
    plt.savefig(kde_plot, dpi=dpi)

    # write description(quantiles information)
    ser.to_csv(log_file, sep="\t")


main = _main.callback


if __name__ == "__main__":
    _main()
