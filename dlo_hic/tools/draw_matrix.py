# -*- coding: utf-8 -*-

"""
draw_matrix.py
~~~~~~~~~~~~~~

Command line tool for draw matrix.

"""

import sys
import argparse

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from dlo_hic.IO.matrix import load_hicmat


def argument_parser():
    parser = argparse.ArgumentParser(
            description="Command line tool for draw matrix.")
    
    parser.add_argument("input",
            help="The input hicmatrix file(.npz).")
    
    parser.add_argument("output",
            help="The name of output image file. like: hic.jpg")

    parser.add_argument("--chrA_chrB",
            nargs=2,
            help="if specified only draw interaction between chrA and chrB. "
                 "like: '--chrA_chrB chr1 chr2' will only draw interaction "
                 "matrix between chr1 and chr2.")

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


def main(input, output, chra_chrb, transform, colorbar, size, dpi):
    hicmat = load_hicmat(input)
    if chra_chrb:
        chr_a, chr_b = chra_chrb
        hicmat = hicmat[chr_a, chr_b]

    if (transform == 'None' or transform == 'False'):
        transform = None
    else:
        transform = eval("np.{}".format(args.transform))

    hicmat.plot(figsize=size, transform=transform, cbar=colorbar)
    plt.savefig(output, dpi=dpi)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    input = args.input
    output = args.output
    chra_chrb = args.chrA_chrB
    transform = args.transform
    colorbar = args.colorbar
    size = args.size
    dpi = args.dpi

    main(input, output, chra_chrb, transform, colorbar, size, dpi)
