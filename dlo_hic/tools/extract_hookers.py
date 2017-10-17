# -*- coding: utf-8 -*-

"""
extract_hookers.py
~~~~~~~~~~~~~~~~~~

Extract the hookers sequences on both sides of linker sequence.
In DLO HiC, 'hooker' is a piece of short sequence,
which is used for mapping to the genome.


This script accept a fastq or fasta file,
and output two hookers files in fasta or fastq file format.

"""

from __future__ import print_function
import sys
import json
import gzip
import argparse
from itertools import islice
from multiprocessing import Process, Queue

from Bio import SeqIO
import regex

from dlo_hic.seq_utils import reverse_complement as rc


QUEUE_TIME_OUT = 3


def argument_parser():
    parser = argparse.ArgumentParser(
            description="Extract the hookers sequences in both left and right side.")

    subparsers = parser.add_subparsers(
            title="commands",
            dest="command",
            metavar="")

    # "pair-end" sub command
    PE_parser = subparsers.add_parser("PE",
            help="option for pair-end sequencing")
    PE_parser.add_argument("input1",
            help="reads1 fastq file")
    PE_parser.add_argument("input2",
            help="reads2 fastq file")

    # "single-end" sub command
    SE_parser = subparsers.add_parser("SE",
            help="option for single-end sequencing")
    SE_parser.add_argument("input",
            help="The input fastq or fastq.gz file.")

    # other arguments
    parser.add_argument("--out1",
            required=True,
            help="output1: left side hookers file")

    parser.add_argument("--out2",
            required=True,
            help="output2: right side hookers file")

    parser.add_argument("--linker-A",
            dest="linker_a",
            required=True,
            help="linker-A")

    parser.add_argument("--linker-B",
            dest="linker_b",
            required=True,
            help="linker-B")

    parser.add_argument("--mismatch",
            type=int,
            default=3,
            help="threshold of linkers base mismatch number, default 3")

    parser.add_argument("--rest",
            type=str,
            default="A*AGCTT",
            help="The sequence of restriction enzyme recognition site, "
                 "default HindIII: 'A*AGCT*T' ")

    parser.add_argument("--phred",
            type=int,
            choices=[33, 64],
            default=33,
            help="The Phred score encode offset type, 33 or 64. default 33")

    parser.add_argument("--processes", "-p",
            type=int,
            default=1,
            help="Use how many processes do calculation. default 1")

    return parser


def parse_rest(rest_str):
    """ parse restriction enzyme site sequence,
    return left, nick and right sequence. """
    left, nick, right = regex.split("[*^]", rest_str)
    return left, nick, right


def load_linkers(linker_a, linker_b):
    """
    compose linker A-A, A-B, B-A, B-B
    """
    linkers = {}
    linkers['A-A'] = linker_a + linker_a
    linkers['A-B'] = linker_a + linker_b
    linkers['B-A'] = linker_b + linker_a
    linkers['B-B'] = linker_b + linker_b
    return linkers


def log_linkers(linkers):
    print("linkers:", file=sys.stderr)
    for key, linker in linkers.items:
        print("{}\t{}".format(key, linker), file=sys.stderr)


def open_file(fname, mode='r'):
    if fname.endswith('.gz'):
        fh = gzip.open(fname, mode=mode+'b')
    else:
        fh = open(fname, mode=mode)
    return fh


def read_fastq(fastq_iter, n=1000):
    """
    read fastq file, return records chunk
    default chunk size 10000
    """
    chunk_iter = itertools.islice(fastq_iter, n)
    chunk = list(chunk_iter)
    if len(chunk) == 0:
        raise StopIteration
    return chunk


def output(fastq_writer, output_queue):
    """output extracted results"""
    while 1:
        try:
            records = output_queue.get(timeout=QUEUE_TIME_OUT)
            fastq_writer.write_records(records)
        except Empty:
            break


def extract_hooker(record, matched_start, matched_end):
    """
    extract hooker in matched record.
    """
    hooker = record[matched_end:]
    return hooker


def worker(records, linkers, rest):
    pass


def stream_processing(fastq_iter, fastq_writer, processes):
    task_queue = Queue()
    output_queue = Queue()
    workers = [Process(target=worker, args=())
            for i in range(processes)]
    output_p = Process(target=output, 
            args=(fastq_writer, output_queue))

    for w in workers:
        w.start()
    output_p.start()

    try:
        while 1:
            task_queue.put(read_fastq(file_in))
    except StopIteration:
        pass

    for w in worker:
        w.join()
    output_p.join()


def mainPE(input1, input2, out1, out2,
        linekr_a, linker_b,
        mismatch, rest, phred, processes):
    # parse restriction enzyme site
    rest_site = parse_rest(rest)

    # load linkers
    linkers = load_linkers(linekr_a, linker_b)
    log_linkers(linkers)

    # open input file
    file_in1 = open_file(input1)
    file_in2 = open_file(input2)
    # open output files
    file_out1 = open_file(out1, "w")
    file_out2 = open_file(out2, "w")

    if phred == 33:
        fastq_iter_1 = SeqIO.parse(file_in_1, 'fastq')
        fastq_iter_2 = SeqIO.parse(file_in_2, 'fastq')
        fastq_writer_1 = SeqIO.QualityIO.FastqPhredWriter(file_out)
        fastq_writer_2 = SeqIO.QualityIO.FastqPhredWriter(file_out)
    else:
        fastq_iter_1 = SeqIO.parse(file_in_1, 'fastq-illumina')
        fastq_iter_2 = SeqIO.parse(file_in_2, 'fastq-illumina')
        fastq_writer_1 = SeqIO.QualityIO.FastqIlluminaWriter(file_out)
        fastq_writer_2 = SeqIO.QualityIO.FastqIlluminaWriter(file_out)

    fastq_writer_1.write_header()
    stream_processing(fastq_iter_1, fastq_writer_1, processes)
    fastq_writer_1.write_footer()

    fastq_writer_2.write_header()
    stream_processing(fastq_iter_2, fastq_writer_2, processes)
    fastq_writer_2.write_footer()

    # close files
    file_in_1.close()
    file_in_2.close()
    file_out1.close()
    file_out2.close()


def mainSE(input, out1, out2,
        linekr_a, linker_b,
        mismatch, rest, phred, processes):
    # parse restriction enzyme site
    rest_site = parse_rest(rest)

    # load linkers
    linkers = load_linkers(linekr_a, linker_b)
    log_linkers(linkers)

    # open input file
    file_in = open_file(input1)
    # open output files
    file_out1 = open_file(out1, "w")
    file_out2 = open_file(out2, "w")

    if phred == 33:
        fastq_iter_1 = SeqIO.parse(file_in, 'fastq')
        fastq_writer_1 = SeqIO.QualityIO.FastqPhredWriter(file_out)
        fastq_writer_2 = SeqIO.QualityIO.FastqPhredWriter(file_out)
    else:
        fastq_iter_1 = SeqIO.parse(file_in, 'fastq-illumina')
        fastq_writer_1 = SeqIO.QualityIO.FastqIlluminaWriter(file_out)
        fastq_writer_2 = SeqIO.QualityIO.FastqIlluminaWriter(file_out)

    fastq_iter_2 = reverse_complement_iter(fastq_iter_1)

    fastq_writer_1.write_header()
    stream_processing(fastq_iter_1, fastq_writer_1, processes)
    fastq_writer_1.write_footer()

    fastq_writer_2.write_header()
    stream_processing(fastq_iter_2, fastq_writer_2, processes)
    fastq_writer_2.write_footer()

    # close files
    file_in_1.close()
    file_out1.close()
    file_out2.close()


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    out1, out2 = args.out1, args.out2
    linker_a = args.linker_a
    linker_b = args.linker_b
    rest = args.rest
    mismatch = args.mismatch
    outtype = args.outtype
    gzip = args.gzip
    phred = args.phred
    processes = args.processes

    command = args.command
    if command == 'SE':
        input = args.input
        mainSE(input, out1, out2,
                linker_a, linker_b,
                rest, mismatch, gzip, phred, processes)

    elif command == 'PE':
        input1 = args.input1
        input2 = args.input2
        mainPE(input1, input2, out1, out2,
                linker_a, linker_b,
                rest, mismatch, gzip, phred, processes)
