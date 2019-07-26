import os
import logging

log = logging.getLogger(__name__)

import click
import matplotlib.pyplot as plt

from dlo_hic.utils.plot.contactmap import *


@click.command(name="plot_contact_map")
@click.argument("path", type=click.Path())
@click.argument("label")
@click.option("--binsize",
    default='auto', 
    help="Binsize(resolution) of output figure's matrix. Use 'auto' to inference it automatically.")
@click.option("--balance",
    is_flag=True,
    help="Balance the output matrix or not.")
@click.option("--global-mat", "-g",
    is_flag=True,
    help="Plot the global contact map(all chromosomes together).")
@click.option("--chrom-list", "-l",
    default="MAIN",
    show_default=True,
    help="The chromosomes to be ploted, use ',' split them, like: 'chr1,chr2,chr3' ."
         "Use 'MAIN' to indecate all main chromosomes(chr1-chrX).")
@click.option("--fig-format", '-f',
    default='.png',
    show_default=True,
    help="Output figure format, like '.png', '.svg' ...")
@click.option("--figure-size", '-s',
    default="10,10",
    help="Output figure size, format in width,height, like: '10,10'")
@click.option("--dpi", default=600,
    help="Dpi of output figure")
def _main(path, label, binsize, balance, global_mat, chrom_list, fig_format, figure_size, dpi):
    """
    \b
    Draw the Hi-C contact map.

    \b
    Args
    ----
    path : str
        Path to input contact map file. Input file should
        in '.mcool' (Multi-resolution Cooler file) or '.hic' file format.

    label : str
        Label of output figure(single figure) or
        directorie(multiple figure of each chromosomes).
    """
    if chrom_list.upper() == 'MAIN':
        chroms = MAIN_CHROM
    else:
        chroms = chrom_list.split(",")

    path = str(path)
    log.info("input:\t" + path)
    log.info("chromosomes:\t" + ",".join(chroms))
    log.info("balance:\t" + str(balance))
    log.info("fig size:\t" + figure_size)
    log.info("fig format:\t" + fig_format)
    log.info("DPI:\t" + str(dpi))

    fig_size = tuple([float(i) for i in figure_size.split(",")])

    if global_mat:
        log.info("Plot global contact map.")
        fig, ax = plot_global_map(path, binsize, balance, chroms)
        outfig = "".join([label, fig_format])
        fig.set_size_inches(fig_size)
        plt.savefig(outfig, dpi=dpi)
        log.info("Figure save to " + outfig)
        plt.close(fig)
    else:
        log.info("Plot each chromosome's contact map.")
        outpath = label
        if not os.path.exists(outpath):
            os.mkdir(outpath)
        for chrom in chroms:
            outfig = os.path.join(outpath, "".join([chrom, fig_format]))
            fig, ax = plot_chrom(path, chrom, binsize, balance)
            fig.set_size_inches(fig_size)
            plt.savefig(outfig, dpi=dpi)
            plt.close(fig)
        log.info("Figures save to dir: " + outpath)


main = _main.callback


if __name__ == "__main__":
    eval("_main()")