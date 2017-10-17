# -*- coding: utf-8 -*-

"""
ice_adjust.py
~~~~~~~~~~~~~

Perform ICE(Iterative Correction and Eigenvector decomposition) on hic result matrix.

"""

import sys
import argparse

import matplotlib
matplotlib.use("AGG") # use AGG backend, avoid system don't have GUI.
import matplotlib.pyplot as plt
import numpy as np
from iced import filter, normalization

from hicmatrix import HicMatrix, HicChrMatrix


def argument_parser():
    parser = argparse.ArgumentParser(" Perform ICE(Iterative Correction and Eigenvector decomposition) on hic result matrix.")

    parser.add_argument("--input", "-i",
            required=True,
            help="prefix of input matrix files: .npy and [.pyobj] file")

    parser.add_argument("--output", "-o",
            help="prefix of output files.")

    parser.add_argument("--percentage", "-p",
            type=float,
            default=0.05,
            help="The threshold for filter low counts, default 0.05")

    parser.add_argument("--figure", "-f",
            help="Output matrix figure if specified.")

    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()

    hicmat = HicMatrix.load(args.input)
    hicmat.matrix = hicmat.matrix.astype(np.float)
    if isinstance(hicmat, HicChrMatrix):
        filtered_mat = filter.filter_low_counts(hicmat.matrix,\
                lengths=hicmat.lengths, percentage=args.percentage)
    else:
        filtered_mat = filter.filter_low_counts(hicmat.matrix,\
                percentage=args.percentage)
    normed_mat = normalization.ICE_normalization(filtered_mat)
    hicmat.matrix = normed_mat

    HicMatrix.save(hicmat, args.output)

    if args.figure:
        hicmat.plot()
        plt.savefig(args.figure, dpi=500)


if __name__ == "__main__":
    main()
