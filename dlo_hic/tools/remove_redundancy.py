# -*- coding: utf-8 -*-

"""
remove_redundancy.py
~~~~~~~~~~~~~~~~~~~~

Remove the redundancy within pairs.

If pairs both ends's distance,
small than the threshold distance at same time,
consider them as the redundancy.

for example:
    reads1    <---  ...  --->
           |-----|      |-----|
    reads2   <--- ... --->

    reads2 can be consider as the replection of reads1.
    reads2 will be remove.


"""

import sys
import argparse
import subprocess

from dlo_hic.utils import read_args
from dlo_hic.utils.tabix_wrap import sort_bedpe_reads1
from dlo_hic.utils.parse_text import parse_line_bedpe

class Bedpe:
    """ The abstract of bedpe record. """
    def __init__(self, line):
        self.line = line.strip()
        items = parse_line_bedpe(line)
        self.items = items
        self.chr1 = items[0]
        self.chr2 = items[3]
        self.start1, self.start2 = items[1], items[4]
        self.end1, self.end2 = items[2], items[5]
        self.center1 = (self.start1 + self.end1) // 2
        self.center2 = (self.start2 + self.end2) // 2

    def is_rep_with(self, another, dis=50):
        """ Judge another bedpe record is replection of self or not. """
        if (self.chr1 != another.chr1) or (self.chr2 != another.chr2):
            return False
        else:
            center1_close = abs(self.center1 - another.center1) <= dis
            center2_close = abs(self.center2 - another.center2) <= dis
            return center1_close and center2_close

    def __str__(self):
        return self.line



def argument_parser():
    parser = argparse.ArgumentParser(
            description="Remove the redundancy within pairs.")

    parser.add_argument("input",
            help="Input bedpe file.")

    parser.add_argument("output",
            help="Output bedpe file.")

    parser.add_argument("--distance", "-d",
            type=int,
            default=50,
            help="The threshold of distance, if pairs both ends's distance,"
                 "small than this at same time, consider them as the redundancy.")

    return parser


def main(input, output, distance):
    # sort input file firstly
    tmp = input + '.tmp'
    sort_bedpe_reads1(input, tmp)

    with open(tmp, 'r') as f, open(output, 'w') as fo:
       base = Bedpe(f.readline())
       while True:
            for line in f:
                another = Bedpe(line)
                if base.is_rep_with(another, distance): # is replication, check next line.
                    continue
                else: # not replication, output base line and change base line.
                    out_line = str(base) + "\n"
                    fo.write(out_line)
                    base = another 
                    break
            else: # arrive at end of file.
                out_line = str(base) + "\n"
                fo.write(out_line)
                break

    subprocess.check_call(['rm', tmp]) # remove tmp file


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    read_args(args, globals())

    main(input, output, distance)