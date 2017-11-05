"""
Extract all restriction sites from fasta file, save to BED6 file format.
"""

from __future__ import print_function
import sys
import re
import argparse
from multiprocessing import Process, Queue
from Queue import Empty
import subprocess
import tempfile

import pyfaidx

from dlo_hic.utils import read_args
from dlo_hic.utils import reverse_complement as rc


TIME_OUT = 1


def argument_parser():
    parser = argparse.ArgumentParser(
        description='Extract all restriction sites from fasta file, save to BED6 file format.')

    parser.add_argument('fasta', help='The input fasta file.')

    parser.add_argument('--rest-seq', '-r',
        dest="rest", required=True,
        help='The sequence of restriction site')

    parser.add_argument('output', help='Name for the resulting bed file.')

    parser.add_argument("--processes", "-p",
        type=int,
        default=1,
        help='Use how many processes to run.')

    return parser


def worker(task_queue, output_queue, rest, fasta):
    while 1:
        try:
            chr_ = task_queue.get(timeout=TIME_OUT)
        except Empty:
            break
        print(chr_, file=sys.stderr)
        seq = fasta[chr_][:].seq
        for match in re.finditer(rest, seq, re.IGNORECASE):
            output_queue.put(
                (chr_, str(match.start()), str(match.end()), '.', '0', '+')
            )
        rc_rest = rc(rest)
        if rc_rest != rest: # find reverse complement restriction site
            for match in re.finditer(rest_rc, seq, re.IGNORECASE):
                output_queue.put(
                    (chr_, str(match.start()), str(match.end()), '.', '0', '-')
                )


def outputer(output_file, output_queue):
    """ output extracted results """
    while 1:
        try:
            record = output_queue.get(timeout=TIME_OUT)
        except Empty:
            break
        line = "\t".join(record) + "\n"
        output_file.write(line)
    output_file.flush()


def main(fasta, rest, output, processes):
    fasta = pyfaidx.Fasta(fasta)
    chrs = fasta.keys()
    task_queue   = Queue()
    output_queue = Queue()
    processes = min(processes, len(chrs))
    workers = [Process(target=worker,
                       args=(task_queue, output_queue, rest, fasta))
               for i in range(processes)]

    for chr_ in chrs:
        task_queue.put(chr_)

    with tempfile.NamedTemporaryFile() as tmp:
        output_p = Process(target=outputer, args=(tmp, output_queue))

        for w in workers:
            w.start()
        output_p.start()

        for w in workers:
            w.join()
        output_p.join()

        # sort output bed file
        print("sorting bed file ...", file=sys.stderr)
        cmd = "sort -k1,1 -k2,2n -u {} > {}".format(tmp.name, output)
        subprocess.check_call(cmd, shell=True)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()
    read_args(args, globals())
    main(fasta, rest, output, processes)