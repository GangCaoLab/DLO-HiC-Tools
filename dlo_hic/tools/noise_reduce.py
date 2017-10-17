# -*- coding: utf-8 -*-

"""
noise_reduce.py
~~~~~~~~~~~~~~~

Remove DLO-HiC noise (self-ligation).

There are two kinds of self-ligation in intra chromosomal interacting pair:
    1. circle ligation.
    2. non-interacting restriction enzyme cutting site.

In first situation, 
there shold no restriction site within pair, in the genome. 
But sometime resriction enzyme can not cutting complemently, 
so maybe there are few restriction sites within pair:

normal:
                         many restriction site
             
               left hooker     |  ...  |      right hooker
                    V          V       V            V     
    genome    <-+++..+++--..---*---..--*----..--+++..+++->

abnormal(self-ligation-1):
                        no or few restriction site
    
               left hooker     |              right hooker
                    V          V                    V     
    genome    <-+++..+++--..---*---..-------..--+++..+++->


In second situation,
left and right hookers will overlap like:

abnormal(self-ligation-2):
            
                              +++++++..++ right hooker
             left hooker ++..+++++++
    genome    <---..--------------------------------..--->

Therefore, with the propose of clean the noise, here shold remove the pair,
which there are no restriction site or just very few restriction(e.g. one)
within it, or the left and right hookers overlapped.

"""


import sys
import argparse


def argument_parser():
    parser = argparse.ArgumentParser(
            description="Remove DLO-HiC noise (self-ligation).")

    parser.add_argument("--input", "-i",
            type=argparse.FileType(mode="r"),
            default=sys.stdin,
            help="Input file, in bedpe file format.")

    parser.add_argument("--output", "-o",
            type=argparse.FileType(mode="w"),
            default=sys.stdout,
            help="Output file, default stdout.")

    parser.add_argument("--restriction", "-r",
            required=True,
            type=argparse.FileType(mode="r"),
            help="bed file which recorded the position of all restriction sites"
                 "in the genome.")

    parser.add_argument("--threshold_num_rest", "-n",
            type=int,
            default=1,
            help="Threshold of number of restriction sites with pair, default 1.")

    parser.add_argument("--threshold_span", "-s",
            type=int,
            default=2000,
            help="Threshold of span, check the pair is noise or not,"
                 "only the pair span less than this parameter, default 2000.")

    parser.add_argument("--debug",
            action="store_true",
            help="If specified, program will write abnormal pairs to stderr.")

    return parser


def load_rest_sites(sites_file):
    """ load all restriction site """
    rest_sites = {}
    with sites_file as f:
        for line in f:
            line = line.strip()
            chr_, start, end, _, _, strand = line.split("\t")
            start, end = int(start), int(end)
            assert start <= end
            rest_sites.setdefault(chr_, list())
            rest_sites[chr_].append((start, end))
    return rest_sites


def find_interval_rests(span, chr_, rest_sites):
    """ find all rests within span. """
    intervals = []
    try:
        rest_sites = rest_sites[chr_]    
    except KeyError as e:
        sys.stderr.write(e)
        raise IOError("There are no chr: {} in resriction bed file.".format(chr_))
    for start, end in rest_sites:
        if start > span[0] and end < span[1]:
            intervals.append((start, end))
    return intervals


def stream_processing(file_in, file_out,
        rest_sites,
        threshold_num_rest,
        threshold_span,
        debug=False):
    """ processing input file and output result """
    with file_in as fi, file_out as fo:
        for line in fi:
            chr1, start1, end1, chr2, start2, end2 = line.strip().split("\t")[:6]
            start1, end1, start2, end2 = [int(i) for i in 
                    (start1, end1, start2, end2)]
            if chr1 != chr2:
                # inter chromalsomal interaction
                fo.write(line)
            else:
                # inter chromalsomal interaction
                if end1 < start2:
                    """ 
                              hooker1         hooker2
                    <--...--s1--..--e1--..--s2--..--e2--..-->
                    """
                    if start2 - end1 > threshold_span:
                        # safe span: interaction distance enough long
                        fo.write(line)
                    else:
                        # unsafe span, need to be check
                        interval_rests = find_interval_rests((end1, start2), chr1, rest_sites)
                        if len(interval_rests) <= threshold_num_rest:
                            fo.write(line)
                        elif debug:
                            sys.stderr.write("[abnormal-1]\t" + line)
                elif end2 < start1:
                    """ 
                              hooker2        hooker1
                    <--...--s2--..--e2--..--s1--..--e1--..-->
                    """
                    if start1 - end2 > threshold_span:
                        # safe span: interaction distance enough long
                        fo.write(line)
                    else:
                        # unsafe span, need to be check
                        interval_rests = find_interval_rests((end2, start1), chr1, rest_sites)
                        if len(interval_rests) <= threshold_num_rest:
                            fo.write(line)
                        elif debug:
                            sys.stderr.write("[abnormal-1]\t" + line)
                else:
                    """ 
                    overlap:

                              hooker1         hooker2
                    <--...--s1--..--s2--..--e1--..--e2--..-->

                    et.al

                    """
                    if debug:
                        sys.stderr.write("[abnormal-2]\t" + line)
                    else:
                        continue


def main():
    parser = argument_parser()
    args = parser.parse_args()
    # load restriction sites
    rest_sites = load_rest_sites(args.restriction)
    # processing data
    stream_processing(args.input, args.output,
            rest_sites,
            args.threshold_num_rest,
            args.threshold_span,
            args.debug)
    

if __name__ == "__main__":
    main()
