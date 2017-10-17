# -*- coding: utf-8 -*-

"""
creat_matrix.py
~~~~~~~~~~~~~~~

Creat interaction matrix from pairs file in bedpe format.

"""

from __future__ import print_function
import sys
import argparse

import matplotlib.pyplot as plt

import dlo_hic
from dlo_hic.genomes import supported_genomes
from dlo_hic.utils.parse_text import is_comment, parse_line_bedpe, parse_line_short


def argument_parser():
    parser = argparse.ArgumentParser(
            description = "Creat interaction matrix from interaction file format.")

    parser.add_argument("input",
            help="Input pairs file in bedpe format, default stdin.")

    parser.add_argument("output_prefix",
            help="Output matrix file prefix.")

    parser.add_argument("--binsize", "-b",
            type=float,
            default=1.0,
            help="Length of bins, unit in 'M'(100000), default 1")

    parser.add_argument("--genome", "-g",
            required=True,
            help="The reference genome you used, like 'hg19', 'hg38'... \n"
                 "or use the "
                 "tab split file which storege the length of each chromosome.")

    parser.add_argument("--figure", "-f",
            help="Output matrix figure if specified.")

    parser.add_argument("--format",
            default='bedpe',
            choices=['bedpe', 'short'],
            help="bedpe(default) or short\n"
                 "bedpe:\n"
                 "\"<chr_a>\t<start_a>\t<end_a>\t<chr_b>\t<start_b>\t<end_b>\t<val>\"\n"
                 "short:\n"
                 "\"<chr_a>\t<pos_a>\t<chr_b>\t<pos_b>\t<val>\"\n")

    return parser


def main(input, output_prefix, binsize, genome, figure, format_):
    binsize = int(binsize * 1000000)

    if genome in supported_genomes:
        chr_len_file = supported_genomes
    else:
        chr_len_file = genome

    chr_len = dlo_hic.IO.load_chr_len(chr_len_file)

    if format_ == 'bedpe':
        parse_func = parse_line_bedpe
    elif format_ == 'short':
        parse_func = parse_line_short

    print("chromosome axis: " + str(hicmat.axis))
    print("total bins:" + str(hicmat.num_bins))
    with open(input) as f:
        for line in f:
            if is_comment(line):
                continue
            chr_a, pos_a, chr_b, pos_b, val = parse_func(line)
            try:
                hicmat.locate(chr_a, pos_a, chr_b, pos_b, val)
            except KeyError as ke:
                sys.stderr.write("unknow chr:" + str(ke) + "\n")
                continue
    if figure:
        hicmat.plot()
        plt.savefig(figure, dpi=600)
    hicmat.save(output_prefix)
    

if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    input = args.input
    output_prefix = args.output_prefix
    binsize = args.binsize
    genome = args.genome
    figure = args.figure
    format_ = args.format

    main(input, output_prefix, binsize, genome, figure, format_)
