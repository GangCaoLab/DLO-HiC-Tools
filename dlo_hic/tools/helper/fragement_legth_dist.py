import os
import io
import sys
import gzip
import itertools
import logging

import click
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import dlo_hic.utils.plot.style as style


log = logging.getLogger(__name__)


def read_rest_sites(fname):
    """
    read restrition sites file
    file must be sorted according to column1 and column2
    (sort -k1,1 -k2,2n)
    return all chromosome's resitriction site positions
    """
    chr2positions = {}

    def read_file(f):
        for line in f:
            items = line.strip().split()
            chr_, pos = items[0], int(items[1])
            chr2positions.setdefault(chr_, [])
            chr2positions[chr_].append(pos)

    if fname.endswith(".gz"):
        # gzip file
        with io.TextIOWrapper(gzip.open(fname, 'r')) as f:
            read_file(f)
    else:
        # normal text file
        with open(fname) as f:
            read_file(f)

    return chr2positions


def possible_frag_len(positions):
    """
    get all possible fragments length of a chromosome
    :positions: all resition site's position
    """
    #comb = itertools.combinations(positions, 2)
    #posi_len = [abs(i - j) for i, j in comb]
    posi_len = np.abs(np.diff(np.asarray(positions)))
    return posi_len


def get_all_frag_len(chr2positions):
    chr2lens = {}
    for chr_, poses in chr2positions.items():
        chr2lens[chr_] = possible_frag_len(poses)
    return chr2lens


def save_frag_len(fname, chr2lens):
    """
    save posible lengthes to file
    """
    with open(fname, 'w') as f:
        for chr_, lens in chr2lens.items():
            lens = [str(l) for l in lens]
            line = chr_ + "\t" + "\t".join(lens) + "\n"
            f.write(line)


def load_frag_len(fname):
    """
    load fragment lengthes of eack chr
    """
    chr2lens = {}
    with open(fname) as f:
        for line in f:
            items = line.strip().split()
            chr_ = items[0]
            lens = [int(i) for i in items[1:]]
            chr2lens[chr_] = lens
    return chr2lens


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


@click.command(name="fragment_length_dist")
@click.argument("rest-sites-files", nargs=-1,
    type=click.Path())
@click.option("--labels",
    help="Label list of each plot, like HindIII,MseI,MbolI ")
@click.option("--keep/--no-keep", default=False,
    help="Keep fragment lengths file. default False")
@click.option("--log-file", default="frag_len.txt",
    help="Sperate file for record statistic information of fragements length. " +\
         "default 'frag_len.txt', use '-' to stdout.")
@click.option("--fig-type",
    type=click.Choice(['box', 'kde']),
    default="box",
    help="Figure type, box plot or kde density plot.")
@click.option("--figure", '-f',
    help="The name of output figure. If not specified will show with `pyplot.show`")
@click.option("--dpi", default=600,
    help="dpi of output figure")
@click.option("--rmext/--no-rmext", default=False,
    help="Remove extream values(values are > 0.95 or < 0.05 quantile).")
@click.option("--right-lim", default=2000,
    help="Right limit for kde plot on x axis. default 2000")
def _main(rest_sites_files, labels,
          keep, log_file,
          fig_type, figure, dpi,
          rmext, right_lim):
    """
    Draw the distribution figure(kde/box plot),
    and output statistic information(quantiles) of fragment legth.
    """
    if fig_type == 'kde' and not rmext:
        log.warning("kde plot should remove extream values. try '--rmext' option.")

    if not rest_sites_files:
        sys.exit(1)

    # process "," joined arg list
    if len(rest_sites_files) == 1 and ',' in rest_sites_files[0]:
        rest_sites_files = rest_sites_files[0].split(",")

    # process label list
    if labels:
        labels = labels.split(",")
        if len(labels) != len(rest_sites_files):
        #
        # if labels can not match to files
        # use file name as it's label
            labels = rest_sites_files
    else:
        labels = rest_sites_files

    fig, ax = plt.subplots()

    all_ser = [] # list for store all pd.Series
    for fname, label in zip(rest_sites_files, labels):

        log.info("dealing %s(label:%s) ..."%(fname, label))

        if os.path.exists(fname + '.frag_len.txt'):
            #
            # if .frag_len.txt file exists, load it directly
            frag_len_f = fname+".frag_len.txt"
            msg = "fragment length file %s founded, load it directly."%frag_len_f
            log.info(msg)
            chr2lens = load_frag_len(fname + '.frag_len.txt')
        else:
            chr2poses = read_rest_sites(fname)
            chr2lens = get_all_frag_len(chr2poses)

        if keep:
            # save fragment lengthes to file
            out_f = fname + '.frag_len.txt'
            log.info("store fragment length infomation to %s"%out_f)
            save_frag_len(out_f, chr2lens)

        all_lens = [] # merge each chr's frag len together
        for chr_, lens in chr2lens.items():
            all_lens.extend(lens)
        all_ser.append(pd.Series(all_lens, name=label))

    # output statistic description of fraglen datasets
    infos = [s.describe() for s in all_ser]
    log_statistic_info(infos)
    log_statistic_info(infos, log_file=log_file)

    df = pd.DataFrame({s.name:s for s in all_ser}) # Series list to Dataframe
    for k in list(df): # iter over Dataframe remove extream value and plot kde
        s = df[k]
        if rmext:
            #
            # remove extreme value(> 0.95 or < 0.05)
            buttom = s.quantile(0.05)
            ceil = s.quantile(0.95)
            s = s[s < ceil]
            s = s[s > buttom]
            df[k] = s
        if fig_type == 'kde':
            sns.kdeplot(s, ax=ax)

    if fig_type == 'kde':
        ax.set_xlim(left=0, right=right_lim)
    elif fig_type == 'box':
        sns.boxplot(data=df, ax=ax)
        ax.set_yscale('log')

    if figure:
        plt.savefig(figure, dpi=dpi)
    else:
        plt.show()


main = _main.callback


if __name__ == "__main__":
    eval("_main()")