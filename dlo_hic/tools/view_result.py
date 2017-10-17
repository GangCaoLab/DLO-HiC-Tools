# -*- coding: utf-8 -*-

"""
view_result.py
~~~~~~~~~~~~~~

Command line tool for generate hic result matrix more easily.

"""

import sys
import argparse

import numpy as np
import matplotlib
matplotlib.use("AGG")
import matplotlib.pyplot as plt

from hicmatrix import HicMatrix, HicChrMatrix


def argument_parser():
    parser = argparse.ArgumentParser(
            description="Command line tool for generate hic result matrix more easily.")
    
    parser.add_argument("--input", "-i",
            required=True,
            help="The prefix of input hicmatrix files(.npy/.npz and .pyobj).")
    
    parser.add_argument("--output", "-o",
            required=True,
            help="The name of output image file.")
    
    parser.add_argument("--size", "-s",
            nargs=2,
            type=int,
            default=(20, 14),
            help="The size of output figure in inch, e.g."
                 "output a 10x20 figure, do like: --size 10 20"
                 "default 20x14")

    parser.add_argument("--dpi", "-d",
            type=int,
            default=300,
            help="The dpi(dots per inch) of output figure, default 300")

    parser.add_argument("--chromosome", "-c",
            help="view one chromosome if specified.")

    parser.add_argument("--cross",
            nargs=2,
            help="view the interact matrix of first chromosome and second chromosome."
                 "e.g. view interaction between chr1 and chr2: --cross chr1 chr2")

    parser.add_argument("--transform", "-t",
            help="transform scale, options: log, log2, log10 etc")

    parser.add_argument("--colorbar",
            action="store_true",
            help="show colorbar if specified.")

    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    hicmat = HicChrMatrix.load(args.input)
    if args.chromosome:
        hicmat = hicmat[args.chromosome, args.chromosome]
    elif args.cross:
        chr_a, chr_b = args.cross
        hicmat = hicmat[chr_a, chr_b]

    if args.transform:
        trans = eval("np.{}".format(args.transform))
    else:
        trans=False

    hicmat.plot(figsize=args.size, transform=trans, cbar=args.colorbar)
    plt.savefig(args.output, dpi=args.dpi)


if __name__ == "__main__":
    main()
