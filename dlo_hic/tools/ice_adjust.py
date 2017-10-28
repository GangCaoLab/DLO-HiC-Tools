# -*- coding: utf-8 -*-

"""
ice_adjust.py
~~~~~~~~~~~~~

Perform ICE(Iterative Correction and Eigenvector decomposition) on hic result matrix.

"""

import sys
import argparse

import matplotlib
import matplotlib.pyplot as plt

from dlo_hic.IO.matrix import load_hicmat
from dlo_hic.normalization import ice


def argument_parser():
    parser = argparse.ArgumentParser(" Perform ICE(Iterative Correction and Eigenvector decomposition) on hic result matrix.")

    parser.add_argument("input",
            help="input matrix file(.npz) name")

    parser.add_argument("output_prefix",
            help="prefix of output file.")

    parser.add_argument("--percentage", "-p",
            type=float,
            default=0.05,
            help="The threshold for filter low counts, default 0.05")

    parser.add_argument("--figure", "-f",
            help="Output matrix figure if specified.")

    return parser


def main(input, output_prefix, percentage, figure):
    hicmat = load_hicmat(args.input)
    ice(hicmat, percentage=percentage)

    hicmat.save(output_prefix)

    if figure:
        hicmat.plot()
        plt.savefig(figure, dpi=600)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    input = args.input
    output_prefix = args.output_prefix
    percentage = args.percentage
    figure = args.figure

    main(input, output_prefix, percentage, figure)
