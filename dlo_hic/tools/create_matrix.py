# -*- coding: utf-8 -*-

"""
creat_matrix.py
~~~~~~~~~~~~~~~

Creat interaction matrix from pairs file in bedpe/short format.

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
            help="Input pairs file in bedpe format.")

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

    parser.add_argument("--chromosome",
            help="Only generate one chromosome's interaction matrix, if specified.")

    return parser


def main(input, output_prefix, binsize, genome, figure, format_, chromosome):
    binsize = int(binsize * 1000000)

    if genome in supported_genomes:
        chr_len = supported_genomes[genome]
    else:
        chr_len_file = genome
        chr_len = dlo_hic.IO.load_chr_len(chr_len_file)

    # select line parse function
    if format_ == 'bedpe':
        parse_func = parse_line_bedpe
    elif format_ == 'short':
        parse_func = parse_line_short

    if chromosome: # select target chrosomosome
        chr_len = [(chr_, len_) for (chr_, len_) in chr_len if chr_ == chromosome]
    hicmat = dlo_hic.HicChrMatrix(chr_len, binsize)

    print("chromosome axis: " + str(hicmat.axis))
    print("total bins:" + str(hicmat.num_bins))

    with open(input) as f:
        for line in f:
            if is_comment(line):
                continue
            chr_a, pos_a, chr_b, pos_b, val = parse_func(line)

            if chromosome: # skip other chromosome
                if chr_a != chromosome or chr_b != chromosome:
                    continue

            try:
                hicmat.locate(chr_a, pos_a, chr_b, pos_b, val)
            except KeyError as ke:
                print("unknow chr:" + str(ke), file=sys.stderr)
                continue
            except IndexError as ie:
                print(str(ie), file=sys.stderr)
                #print(line.strip(), file=sys.stderr)
                continue

    hicmat.save(output_prefix)

    if figure:
        hicmat.plot()
        plt.savefig(figure, dpi=600)
    

if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    input = args.input
    output_prefix = args.output_prefix
    binsize = args.binsize
    genome = args.genome
    figure = args.figure
    format_ = args.format
    chromosome = args.chromosome

    main(input, output_prefix, binsize, genome, figure, format_, chromosome)
