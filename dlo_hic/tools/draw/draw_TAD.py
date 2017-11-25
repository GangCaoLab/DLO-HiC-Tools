# -*- coding: utf-8 -*-

"""
draw_TAD.py
~~~~~~~~~~~~~~

Command line tool for draw TAD.

"""

import sys
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from dlo_hic.IO.matrix import load_hicmat
from dlo_hic.plot.TAD import plot_TAD


def argument_parser():
    parser = argparse.ArgumentParser(
            description="Command line tool for draw TAD.")
    
    parser.add_argument("input",
            help="The input hicmatrix file(.npz).")
    
    parser.add_argument("output",
            help="The name of output image file. like: TAD.jpg")

    parser.add_argument("chr",
            help="The chromosome of TAD")

    parser.add_argument("--start",
            help="The start position of TAD")

    parser.add_argument("--end",
            help="The end position of TAD")

    parser.add_argument("--transform", "-t",
            default='log10',
            help="transform scale, options: log, log2, log10(default) etc")

    parser.add_argument("--colorbar",
            action="store_true",
            help="show colorbar if specified.")

    parser.add_argument("--size", "-s",
            nargs=2,
            type=int,
            default=(20, 14),
            help="The size of output figure in inch, e.g."
                 "output a 10x20 figure, do like: --size 10 20"
                 "default 20x14")

    parser.add_argument("--dpi", "-d",
            type=int,
            default=600,
            help="The dpi(dots per inch) of output figure, default 300")

    return parser


def main(input, output, chr_, start, end, transform, colorbar, size, dpi):
    hicmat = load_hicmat(input)

    if start:
        start = int(start)
    if end:
        end = int(end)

    if (transform == 'None' or transform == 'False'):
        transform = None
    else:
        transform = eval("np.{}".format(args.transform))

    plot_TAD(hicmat, chr_, start, end, figsize=size, transform=transform, cbar=colorbar)
    plt.savefig(output, dpi=dpi)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    input = args.input
    output = args.output
    chr_ = args.chr
    start = args.start
    end = args.end
    transform = args.transform
    colorbar = args.colorbar
    size = args.size
    dpi = args.dpi

    main(input, output, chr_, start, end, transform, colorbar, size, dpi)
