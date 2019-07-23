import sys

import matplotlib.pyplot as plt
import click

from dlo_hic.utils.infer_adapter import infer_adapter, plot_char_prob_stacked_bar


@click.command(name="infer_adapter")
@click.argument("fq_path", nargs=1,
    type=click.Path())
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
def _main(fq_path, n, start_pos, prob_thresh, figure, fig_size, dpi):
    """
    Inference adapter sequence, and plot the decision graph.
    """
    probs, flag_seq, adapter_seq = infer_adapter(fq_path, n, start_pos, prob_thresh)
    print("Flag sequence:", file=sys.stderr)
    print(flag_seq, file=sys.stderr)
    print("Adapter:", file=sys.stderr)
    print(adapter_seq)
    fig, _ = plot_char_prob_stacked_bar(probs, start_pos)
    fig.set_size_inches(*fig_size)
    if figure:
        plt.savefig(figure, dpi=dpi)
    else:
        plt.show()


main = _main.callback


if __name__ == "__main__":
    eval("_main()")