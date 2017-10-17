# -*- coding: utf-8 -*-

"""
creat_matrix.py
~~~~~~~~~~~~~~~

Creat interaction matrix from pairs file in bedpe format.

"""

import sys
import argparse

import matplotlib
matplotlib.use("AGG") # use AGG backend, avoid system don't have GUI.
import matplotlib.pyplot as plt

from hicmatrix import HicChrMatrix


def argument_parser():
    parser = argparse.ArgumentParser(
            description = "Creat interaction matrix from pairs file in bedpe format.")

    parser.add_argument("--input", "-i",
            type=argparse.FileType(mode="r"),
            default=sys.stdin,
            help="Input pairs file in bedpe format, default stdin.")

    parser.add_argument("--output", "-o",
            required=True,
            help="Output matrix file.")

    parser.add_argument("--binsize", "-b",
            type=float,
            default=1.0,
            help="Length of bins, unit in 'M'(100000), default 1")

    parser.add_argument("--chr_len", "-g",
            required=True,
            type=argparse.FileType(mode="r"),
            help="Tab split file which storege the length of each chromosome.")

    parser.add_argument("--figure", "-f",
            help="Output matrix figure if specified.")

    return parser



def main():
    parser = argument_parser()
    args = parser.parse_args()
    binsize = int(args.binsize * 1000000)
    hicmat = HicChrMatrix.load_from_file(args.chr_len, binsize)
    print("chromosome axis: " + str(hicmat.axis))
    print("total bins:" + str(hicmat.num_bins))
    with args.input as f:
        for line in f:
            items = line.strip().split("\t")
            chr_a, x_s, x_e, chr_b, y_s, y_e = items[:6]
            x_s, x_e, y_s, y_e = [int(i) for i in (x_s, x_e, y_s, y_e)]
            x, y = (x_s + x_e)//2, (y_s + y_e)//2
            try:
                hicmat.locate(chr_a, x, chr_b, y)
            except KeyError as ke:
                sys.stderr.write("unknow chr:" + str(ke) + "\n")
                continue
    if args.figure:
        hicmat.plot()
        plt.savefig(args.figure, dpi=500)
    HicChrMatrix.save(hicmat, args.output)
    

if __name__ == "__main__":
    main()
