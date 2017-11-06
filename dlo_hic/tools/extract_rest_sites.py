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
import signal

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


def worker(task_queue, output_queue, count_queue, rest, fasta):
    rc_rest = rc(rest)
    faidx = pyfaidx.Fasta(fasta)
    while 1:
        c = 0
        try:
            chr_ = task_queue.get(timeout=TIME_OUT)
        except Empty:
            break
        print(chr_, file=sys.stderr)
        seq = faidx[chr_][:].seq
        for match in re.finditer(rest, seq, re.IGNORECASE):
            output_queue.put(
                (chr_, str(match.start()), str(match.end()), '.', '0', '+')
            )
            c += 1
        if rc_rest != rest: # find reverse complement restriction site
            for match in re.finditer(rest_rc, seq, re.IGNORECASE):
                output_queue.put(
                    (chr_, str(match.start()), str(match.end()), '.', '0', '-')
                )
                c += 1
        count_queue.put(c)


def outputer(output_file, output_queue):
    """ output extracted results """
    def signal_handeler(signal, frame):
        output_file.flush()
        sys.exit(0) 
    signal.signal(signal.SIGTERM, signal_handeler)
    while 1:
        record = output_queue.get()
        line = "\t".join(record) + "\n"
        output_file.write(line)


def main(fasta, rest, output, processes):
    faidx = pyfaidx.Fasta(fasta)
    chrs = faidx.keys()
    task_queue   = Queue()
    output_queue = Queue()
    count_queue = Queue()
    processes = min(processes, len(chrs))
    workers = [Process(target=worker,
                       args=(task_queue, output_queue, count_queue, rest, fasta))
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

        c = 0
        while not count_queue.empty():
            c += count_queue.get()
        print("%d restriction sites found"%c, file=sys.stderr)

        while not output_queue.empty():
            time.sleep(TIME_OUT)
        output_p.terminate()

        # sort output bed file
        print("sorting bed file ...", file=sys.stderr)
        cmd = "sort -k1,1 -k2,2n -u {} > {}".format(tmp.name, output)
        subprocess.check_call(cmd, shell=True)


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()
    read_args(args, globals())
    main(fasta, rest, output, processes)