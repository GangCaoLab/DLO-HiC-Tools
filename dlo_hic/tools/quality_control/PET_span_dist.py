import os
from os.path import split, splitext, join
import logging

import numpy as np
import click
import matplotlib.pyplot as plt

from dlo_hic.utils.filetools import open_file, sample_file
from dlo_hic.utils.dataframe import interactions_dataframe


log = logging.getLogger(__name__)


def save_hist(hist, output):
    values = hist[0].tolist()
    ranges = hist[1].tolist()
    ranges = [(ranges[i], ranges[i+1]) for i in range(len(ranges) - 1)]
    with open(output, "a") as f:
        outline = "[Hist]\n"
        f.write(outline)
        outline = "range_start\trange_end\tcount\n"
        for (start, end), cnt in zip(ranges, values):
            outline = "{}\t{}\t{}\n".format(start, end, cnt)
            f.write(outline)
        f.write("\n")


def plot_hist(hist, input, outfig):
    ys = hist[0]
    ranges = hist[1]
    xs = ranges[:-1] + np.diff(ranges)/2

    fig, ax = plt.subplots(figsize=(7, 5))
    xs = xs[ys != 0]
    xs = np.log10(xs)
    ys = ys[ys != 0]
    ys = np.log10(ys)
    plt.scatter(xs, ys, c="#66ccff", alpha=0.7)
    file_name = os.path.split(input)[1]
    plt.title("PET span distribution: {}".format(file_name))
    plt.xlabel("$log_{10}(PET\ span)$")
    plt.ylabel("$log_{10}(number\ of\ PETs)$")
    plt.tight_layout()
    plt.savefig(outfig, dpi=300)


def save_describe(span, output):
    with open(output, 'w') as f:
        outline = "[Describe]\n"
        f.write(outline)
        des = span.describe()
        for key, val in zip(des.index, des):
            outline = "{}\t{}\n".format(key, val)
            f.write(outline)
        f.write("\n")


@click.command(name="PET_span_dist")
@click.argument("input")
@click.argument("output")
@click.argument("outfig")
@click.option("--sample", "-s",
    type=int,
    default=10000,
    show_default=True,
    help="Random sample lines from file. Use -1 to indicate not random sample use whole file.")
@click.option("--hist-bins",
    default=1000,
    show_default=True,
    help="How may bins in hist plot.")
def _main(input, output, outfig, sample, hist_bins):
    """
    Count the distribution of PET span.

    Arguments
    ---------
    input : str
        Path to bedpe file or pairs file, support gziped file.
    output : str
        Path to output file, which store the statistic information.
    outfig : str
        Path to output figure.

    Also plot a log-log distribution fig, store at the output file's dir.

    """

    if sample > 0:
        tmp = input + ".sample"
        sample_file(input, tmp, sample)
        df = interactions_dataframe(input)
    else:
        df = interactions_dataframe(input)

    df['span'] = abs(df.pos2 - df.pos1)
    df[df.chr1 != df.chr2] = np.nan
    #import ipdb; ipdb.set_trace()

    save_describe(df.span, output)

    hist = plt.hist(df.span, bins=hist_bins)
    plt.clf()
    save_hist(hist, output)
    
    plot_hist(hist, input, outfig)

    if sample:
        os.remove(tmp)


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
