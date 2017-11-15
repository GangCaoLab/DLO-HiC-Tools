import sys
import argparse
import subprocess

from dlo_hic.utils import read_args
from dlo_hic.utils.tabix_wrap import sort_pairs, index_pairs
from dlo_hic.utils.parse_text import Bedpe


def argument_parser():
    parser = argparse.ArgumentParser(
            description="transform bedpe format file to pairs format."
                        "about pairs format: \n"
                        "https://github.com/4dn-dcic/pairix/blob/master/pairs_format_specification.md \n"
                        "and index it use pairix")

    parser.add_argument("input",
            help="Input bedpe file.")

    parser.add_argument("output",
            help="Output pairs file.")

    return parser


def bedpe2pairs(input, output):
    with open(input) as fi, open(output, 'w') as fo:
        for line in fi:
            bpe = Bedpe(line)
            pairs_line = bpe.to_pairs_line()
            fo.write(pairs_line + "\n")


def main(input, output):
    # sort input file firstly
    tmp0 = input + '.tmp.0'
    bedpe2pairs(input, tmp0)
    sort_pairs(tmp0, output)
    index_pairs(output)

    subprocess.check_call(['rm', tmp0]) # remove tmp files


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    read_args(args, globals())

    main(input, output)