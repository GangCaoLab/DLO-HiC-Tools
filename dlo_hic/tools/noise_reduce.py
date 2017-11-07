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
But sometime restriction enzyme can not cutting complemently, 
so maybe there are few restriction sites within pair:

normal:
                         many restriction site
             
                 left PET      |  ...  |        right PET
                    V          V       V            V     
    genome    <-+++..+++--..---*---..--*----..--+++..+++->

abnormal(self-ligation-1):
                     no or few restriction site

                 left PET      |                right PET
                    V          V                    V     
    genome    <-+++..+++--..---*---..-------..--+++..+++->


In second situation,
left and right PETs will overlap like:

abnormal(self-ligation-2):
            
                              +++++++..++ right PET
             left PET ++..+++++++
    genome    <---..--------------------------------..--->

Therefore, with the propose of clean the noise, here shold remove the pair,
which there are no restriction site or just very few restriction(e.g. one)
within it, or the left and right PETs overlapped.

"""

from __future__ import print_function
import os
import sys
import time
import cPickle
import argparse
import signal
from Queue import Empty
import multiprocessing
from multiprocessing import Queue, Process

from intervaltree import Interval, IntervalTree

from dlo_hic.utils import read_args
from dlo_hic.utils.parse_text import parse_line_bed6, parse_line_bedpe
from dlo_hic.utils.tabix_wrap import query_bed6


TIME_OUT = 1
CHUNK_SIZE = 1000


def argument_parser():
    parser = argparse.ArgumentParser(
            description="Remove DLO-HiC noise (self-ligation).")

    parser.add_argument("input",
            help="Input file, in bedpe file format.")

    parser.add_argument("output",
            help="Output file, default stdout.")

    parser.add_argument("--restriction", "-r",
            required=True,
            help="bed file which recorded the position of all restriction sites"
                 "in the genome.")

    parser.add_argument("--processes", "-p",
            type=int,
            default=1,
            help="Use how many processes to run.")

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


def find_interval_rests(sites_file, chr_, span):
    """ find all rests within span. """
    start, end = span
    intervals = list(query_bed6(sites_file, chr_, start, end))
    return intervals


def read_chunk(file, chunk_size=CHUNK_SIZE):
    """ read and parse bedpe file,
    return a list of bedpe items
    """
    chunk = []
    for i in range(chunk_size):
        line = file.readline()
        if not line:
            if i == 0:
                raise StopIteration
            else:
                return chunk
        items = parse_line_bedpe(line)
        chunk.append(items)
    return chunk


def bedpe_type(restriction, bedpe_items):
    """ judge interaction type,
    types:
        "normal"
        "abnormal-1"
        "abnormal-2"
    """
    chr1, start1, end1, chr2, start2, end2 = bedpe_items[:6]
    # inter chromalsomal interaction
    if chr1 != chr2:
        return "normal"

    # intra chromalsomal interaction
    if end1 < start2:
        """
                  PET1            PET2
        <--...--s1--..--e1--..--s2--..--e2--..-->
        """
        if start2 - end1 > threshold_span:
            # safe span: interaction distance enough long
            return "normal"
        else:
            # unsafe span, need to be check
            interval_rests = find_interval_rests(restriction, chr1, (end1, start2))
            if len(interval_rests) <= threshold_num_rest:
                return "normal"
            else:
                return "abnormal-1"
    elif end2 < start1:
        """
                  PET2            PET1
        <--...--s2--..--e2--..--s1--..--e1--..-->
        """
        if start1 - end2 > threshold_span:
            # safe span: interaction distance enough long
            return "normal"
        else:
            # unsafe span, need to be check
            interval_rests = find_interval_rests(restriction, chr1, (end1, start2))
            if len(interval_rests) <= threshold_num_rest:
                return "normal"
            else:
                return "abnormal-1"
    else:
        """
        overlap:

                  PET1            PET2
        <--...--s1--..--s2--..--e1--..--e2--..-->

        et.al

        """
        return "abnormal-2"


def worker(task_queue, output_queue, err_queue, restriction):
    while 1:
        try:
            chunk = task_queue.get(timeout=TIME_OUT)
        except Empty:
            break
        out_chunk = []
        err_chunk = []
        for items in chunk:
            type_ = bedpe_type(restriction, items)
            if type_ == 'normal':
                out_chunk.append(items)
            else:
                err_chunk.append((type_, items))
        output_queue.put(out_chunk)
        err_queue.put(err_chunk)


def outputer(out_file, output_p):
    f = open(out_file, 'w')
    def signal_handeler(signal, frame):
        f.close()
        sys.exit(0)
    signal.signal(signal.SIGTERM, signal_handeler)
    while 1:
        chunk = output_p.get()
        for items in chunk:
            items = [str(i) for i in items]
            line = "\t".join(items) + "\n"
            f.write(line)


def main(input, output,
         restriction, processes,
         threshold_num_rest, threshold_span,
         debug):
    task_queue = Queue()
    output_queue = Queue()
    err_queue = Queue()

    workers = [Process(target=worker,
                       args=(task_queue, output_queue, err_queue, restriction))
               for i in range(processes)]

    output_p = Process(target=outputer, args=(output, output_queue))
    
    for w in workers:
        w.start()

    output_p.start()

    with open (input) as f:
        while 1:
            try:
                task_queue.put(read_chunk(f))
            except StopIteration:
                break

    for w in workers:
        w.join()
    while not output_queue.empty():
        time.sleep(TIME_OUT)
    output_p.terminate()

    if debug:
        while not err_queue.empty():
            chunk = err_queue.get()
            for type_, items in chunk:
                items = [str(i) for i in items]
                line = "\t".join(items)
                print("[{}]\t{}".format(type_, line), file=sys.stderr)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()
    read_args(args, globals())
    main(input, output,
         restriction, processes,
         threshold_num_rest, threshold_span,
         debug)