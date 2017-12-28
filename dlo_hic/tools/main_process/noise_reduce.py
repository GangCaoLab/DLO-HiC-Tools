import os
import sys
import time
import logging
import subprocess
import signal
from queue import Empty
import multiprocessing as mp

import click

from dlo_hic.utils import read_args
from dlo_hic.utils.parse_text import parse_line_bed6, parse_line_bedpe
from dlo_hic.utils.wrap.tabix import query_bed6


TIME_OUT = 1
CHUNK_SIZE = 10000

log = logging.getLogger(__name__)


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


def bedpe_type(restriction, bedpe_items, threshold_span, threshold_num_rest):
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
            interval_rests = find_interval_rests(restriction, chr1, (end2, start1))
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


def worker(task_queue,
           restriction,
           output_file, err_file,
           threshold_span, threshold_num_rest):
    current = mp.current_process().pid
    output_f = open(output_file, 'w')
    err_f = open(err_file, 'w')
    while 1:
        try:
            chunk = task_queue.get()
            if chunk is None:
                raise Empty
        except Empty:
            output_f.close()
            err_f.close()
            log.debug("Process-%d done."%current)
            break

        for items in chunk:
            type_ = bedpe_type(restriction, items, threshold_span, threshold_num_rest)
            if type_ == 'normal':
                # normal, will output this line
                out_line = "\t".join([str(i) for i in items]) + "\n"
                output_f.write(out_line)
            else:
                out_line = "\t".join([type_] + [str(i) for i in items]) + "\n"
                err_f.write(out_line)


@click.command(name="noise_reduce")
@click.argument("bedpe")
@click.argument("output")
@click.option("--restriction", "-r",
    required=True,
    help="bed file which recorded the position of all restriction sites"
        "in the genome.")
@click.option("--processes", "-p", 
    default=1,
    help="Use how many processes to run.")
@click.option("--threshold-num-rest", "-n",
    default=1,
    help="Threshold of number of restriction sites with pair, default 1.")
@click.option("--threshold-span", "-s",
    default=50,
    help="Threshold of pair span, default 50.")
def _main(bedpe, output,
         restriction, processes,
         threshold_num_rest, threshold_span):
    """
    Remove DLO-HiC noise (self-ligation).

    \b
    There are two kinds of self-ligation in intra chromosomal interacting pair:
        1. circle ligation.
        2. non-interacting restriction enzyme cutting site.

    \b
    In first situation, 
    there shold no restriction site within pair, in the genome. 
    But sometime restriction enzyme can not cutting complemently, 
    so maybe there are few restriction sites within pair:

    \b
    normal:
                             many restriction site

                     left PET      |  ...  |        right PET
                        V          V       V            V     
        genome    <-+++..+++--..---*---..--*----..--+++..+++->

    \b
    abnormal(self-ligation-1):
                         no or few restriction site

                     left PET      |                right PET
                        V          V                    V     
        genome    <-+++..+++--..---*---..-------..--+++..+++->


    \b
    In second situation,
    left and right PETs will overlap like:

    \b
    abnormal(self-ligation-2):
                                  +++++++..++ right PET
                 left PET ++..+++++++
        genome    <---..--------------------------------..--->

    Therefore, with the propose of clean the noise, here shold remove the pair,
    which there are no restriction site or just very few restriction(e.g. one)
    within it, or the left and right PETs overlapped.

    """

    log.info("noise reduce on tile %s"%bedpe)

    if restriction.endswith(".gz"): # remove .gz suffix
        restriction = restriction.replace(".gz", "")

    task_queue = mp.Queue() 

    workers = [mp.Process(target=worker,
                          args=(task_queue,
                                restriction,
                                output+".tmp.%d"%i, output+".tmp.e.%d"%i,
                                threshold_num_rest, threshold_span))
               for i in range(processes)]

    log.info("%d workers spawned for noise refuce"%len(workers))

    for w in workers:
        w.start()

    with open(bedpe) as f:
        while 1:
            try:
                task_queue.put(read_chunk(f))
            except StopIteration:
                break
    
    for w in workers:
        task_queue.put(None)

    for w in workers:
        w.join()

    log.info("merging tmporary files.")
    # merge tmp files
    tmp_files = [output+".tmp.%d"%i for i in range(processes)]
    cmd = "cat " + " ".join(tmp_files) + " > " + output
    subprocess.check_call(cmd, shell=True)
    cmd = "rm " + " ".join(tmp_files)
    subprocess.check_call(cmd, shell=True)
    # merge error tmp files
    err_files = [output+".tmp.e.%d"%i for i in range(processes)]
    cmd = "cat " + " ".join(err_files) + " > " + output + ".err"
    subprocess.check_call(cmd, shell=True)
    cmd = "rm " + " ".join(err_files)
    subprocess.check_call(cmd, shell=True)


main = _main.callback


if __name__ == "__main__":
    _main()