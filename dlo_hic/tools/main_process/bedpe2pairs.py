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

    parser.add_argument("--keep",
            action="store_true",
            help="keep non compressed pairs file," + \
                 " if you need create .hic file use this option.")

    return parser


def bedpe2pairs(input, output):
    with open(input) as fi, open(output, 'w') as fo:
        for line in fi:
            bpe = Bedpe(line)
            pairs_line = bpe.to_pairs_line()
            fo.write(pairs_line + "\n")


def add_pairs_header(input):
    """ add header to pairs file. """
    header = "## pairs format v1.0\\n" +\
             "#columns: readID chr1 position1 chr2 position2 strand1 strand2"
    tmp0 = input + '.tmp'
    tmp1 = ".header"
    cmd = "echo \"{}\" > {}".format(header, tmp1)
    subprocess.check_call(cmd, shell=True)
    cmd = "cat {} {} > {}".format(tmp1, input, tmp0)
    subprocess.check_call(cmd, shell=True)
    cmd = "mv {} {}".format(tmp0, input)
    subprocess.check_call(cmd, shell=True)
    cmd = "rm {}".format(tmp1)
    subprocess.check_call(cmd, shell=True)


def main(input, output, keep):
    # sort input file firstly
    tmp0 = input + '.tmp.0'
    bedpe2pairs(input, tmp0)
    sort_pairs(tmp0, output)
    subprocess.check_call(['rm', tmp0]) # remove tmp files
    index_pairs(output)
    if keep:
        add_pairs_header(output) # if keep uncompressed file, add header to it
    else:
        subprocess.check_call(['rm', output]) # remove uncompressed file


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    read_args(args, globals())

    main(input, output, keep)