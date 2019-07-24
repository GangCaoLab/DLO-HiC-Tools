import sys
import logging
log = logging.getLogger(__name__)

import matplotlib.pyplot as plt
import click

from dlo_hic.utils.infer_adapter import infer_adapter, plot_char_prob_stacked_bar


@click.command(name="infer_adapter")
@click.argument("fq_path", nargs=1,
    type=click.Path())
@click.option("--log-file", default="./adapter.log",
    help="File for store information about adapter inference.")
@click.option("-n", default=50,
    show_default=True,
    help="Number of fastq records used for inference adapter sequence.")
@click.option("--start-pos", default=70,
    show_default=True,
    help="Start position of adapter search.")
@click.option("--prob-thresh", default=0.75,
    show_default=True,
    help="Threshold of occurace probability in adapter inference.")
@click.option("--figure", '-f',
    help="The name of output figure. If not specified will show with `pyplot.show`")
@click.option("--fig-size", default=(12, 5),
    show_default=True,
    help="Figure size.")
@click.option("--dpi", default=600,
    show_default=True,
    help="dpi of output figure")
def _main(fq_path, log_file, n, start_pos, prob_thresh, figure, fig_size, dpi):
    """
    Inference adapter sequence, and plot the decision graph.
    """
    probs, flag_seq, adapter_seq = infer_adapter(fq_path, n, start_pos, prob_thresh)
    log.info("Flag sequence:")
    log.info(flag_seq)
    log.info("Adapter:")
    log.info(adapter_seq)
    with open(log_file, 'w') as f:
        f.write("flag_seq\t{}\n".format(flag_seq))
        f.write("adapter_seq\t{}\n".format(adapter_seq))
    fig, _ = plot_char_prob_stacked_bar(probs, start_pos)
    fig.set_size_inches(*fig_size)
    if figure:
        plt.savefig(figure, dpi=dpi)
    else:
        plt.show()
    return adapter_seq


main = _main.callback


if __name__ == "__main__":
    eval("_main()")