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

import pandas as pd


class Bedpe:
    """ The abstract of bedpe record. """
    def __init__(self, line):
        self.line = line.strip()
        items = self.line.split("\t")
        self.items = items 
        self.chr1 = items[0]
        self.chr2 = items[3]
        self.start1, self.start2 = int(items[1]), int(items[4])
        self.end1, self.end2 = int(items[2]), int(items[5])
        self.center1 = (self.start1 + self.end1) // 2
        self.center2 = (self.start2 + self.end2) // 2

    def is_rep_with(self, another, dis=50):
        """ Judge another bedpe record is replection of self or not. """
        if (self.chr1 != another.chr1) or (self.chr2 != another.chr2):
            return False
        else:
            return abs(self.center1 - another.center1) <= dis and abs(self.center2 - another.center2) <= dis

    def __str__(self):
        return self.line



def argument_parser():
    parser = argparse.ArgumentParser(
            description="Remove the redundancy within pairs.")

    parser.add_argument("--input", "-i",
            type=argparse.FileType(mode="r"),
            default=sys.stdin,
            help="Input bedpe file, column 1-3 must be sorted.")

    parser.add_argument("--distance", "-d",
            type=int,
            default=50,
            help="The threshold of distance, if pairs both ends's distance,"
                 "small than this at same time, consider them as the redundancy.")

    return parser


def main():
    parser = argument_parser()
    args = parser.parse_args()
    dis = args.distance

    with args.input as f:
        base = Bedpe(f.readline())
        while True:
            for line in f:
                another = Bedpe(line)
                if base.is_rep_with(another, dis): # is replication, check next line.
                    continue
                else: # not replication, output base line and change base line.
                    print(str(base))
                    base = another 
                    break
            else: # arrive at end of file.
                print(str(base))
                break


if __name__ == "__main__":
    main()
