# -*- coding: utf-8 -*-

"""
extract_PETs.py
~~~~~~~~~~~~~~~~~~

Extract the PETs sequences on both sides of linker sequence.
In DLO HiC, 'PET'(Pair End Tag) is a piece of short sequence,
which is used for mapping to the genome.


This script accept a fastq or fasta file,
and output two PETs files in fasta or fastq file format.

"""

from __future__ import print_function
import re
import sys
import json
import gzip
import argparse
import itertools
from math import floor
from copy import copy
from multiprocessing import Process, Queue, Manager

import regex
from Bio import SeqIO
from cutadapt._align import Aligner

from dlo_hic.utils import reverse_complement as rc
from dlo_hic.utils import read_args


TIME_OUT = 1
CHUNK_SIZE = 1000


def argument_parser():
    parser = argparse.ArgumentParser(
            description="Extract the PET(Pair End Tag) sequences in both left and right side.")

    subparsers = parser.add_subparsers(
            title="commands",
            dest="command",
            metavar="")

    def add_arguments(parser):
        parser.add_argument("--out1", '-o1',
                required=True,
                help="output1: left side PET fastq file")

        parser.add_argument("--out2", '-o2',
                required=True,
                help="output2: right side PET fastq file")

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
                default=4,
                help="threshold of linkers base mismatch(and gap open extends) number, default 3")

        parser.add_argument("--allow-gap",
                action="store_true",
                dest="allow_gap",
                help="allow gap open when align linkers to sequence.")

        parser.add_argument("--rest",
                type=str,
                default="A*AGCT*T",
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

    # "pair-end" sub command
    PE_parser = subparsers.add_parser("PE",
            help="option for pair-end sequencing")
    PE_parser.add_argument("input1",
            help="reads1 fastq file")
    PE_parser.add_argument("input2",
            help="reads2 fastq file")
    add_arguments(PE_parser)

    # "single-end" sub command
    SE_parser = subparsers.add_parser("SE",
            help="option for single-end sequencing")
    SE_parser.add_argument("input",
            help="The input fastq or fastq.gz file.")
    SE_parser.add_argument("--PET-len",
            dest="PET_len",
            type=int,
            default=0,
            help="The expected length of PET sequence,"
                 "if PET_len==0 (default) will not limit length.")
    add_arguments(SE_parser)

    return parser


def parse_rest(rest_str):
    """ parse restriction enzyme site sequence,
    return left, nick and right sequence. """
    left, nick, right = re.split("[*^]", rest_str)
    return left, nick, right


def reverse_complement_record(record):
    """ reverse complement a fastq record """
    rc_record = copy(record)
    rc_record.seq = rc(str(rc_record.seq))
    qua = rc_record.letter_annotations['phred_quality']
    qua.reverse()
    return rc_record


def load_linkers(linker_a, linker_b):
    """
    compose linker A-A, A-B, B-A, B-B
    """
    linkers = {}
    linkers['A-A'] = linker_a + rc(linker_a)
    linkers['A-B'] = linker_a + rc(linker_b)
    linkers['B-A'] = linker_b + rc(linker_a)
    linkers['B-B'] = linker_b + rc(linker_b)
    return linkers


def get_linker_regex(linker, mismatch_threshold, exact_flank_len=4):
    """
    get compiled linker regular expression object.
    """
    # construct global linker2regex for storage regex objects
    global_vals = globals()
    if 'linker2regex' not in global_vals:
        global_vals['linker2regex'] = {}

    linker_regex = linker2regex.get((linker, mismatch_threshold, exact_flank_len), False)
    if not linker_regex: # linker regex object haven't built
        ex = exact_flank_len
        if mismatch_threshold == 0:
            linker_regex =  regex.compile("%s"%linker)
        else:
            linker_regex = regex.compile("%s(%s){e<=%s}%s"%(
                linker[:ex],
                linker[ex:-ex],
                mismatch_threshold,
                linker[-ex:]
            ))
        linker2regex[(linker, mismatch_threshold, exact_flank_len)] = linker_regex
    return linker_regex


def log_linkers(linkers, file=sys.stderr):
    print("linkers:", file=file)
    for key, linker in linkers.items():
        print("{}\t{}".format(key, linker), file=file)


def log_counts(counts, file=sys.stderr):
    print("Quality Control:", file=file)
    for k, v in counts.items():
        print("\t{}\t{}".format(k, v), file=file)
    ratio = counts['intra-molecular'] / float(counts['all'])
    print("Valid reads ratio: {}".format(ratio))


def open_file(fname, mode='r'):
    if fname.endswith('.gz'):
        fh = gzip.open(fname, mode=mode+'b')
    else:
        fh = open(fname, mode=mode)
    return fh


def read_fastq(fastq_iter, n=CHUNK_SIZE):
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
    """ output extracted results """
    from Queue import Empty
    while 1:
        try:
            records = output_queue.get(timeout=TIME_OUT)
        except Empty:
            break
        c = fastq_writer.write_records(records)

    fastq_writer.write_footer()
    fastq_writer.handle.close()


def add_base_to_PET(PET, base):
    quality = PET.letter_annotations['phred_quality']
    quality.append(38)
    PET.letter_annotations = {}
    PET.seq = PET.seq + base
    PET.letter_annotations['phred_quality'] = quality


def extract_PET(record, span, rest, PET_len=None, SE=False):
    """
    extract PET in matched record.
    """
    start, end = span
    if PET_len:
        # PET_len decrease 1 because will add one base at end
        PET_len = PET_len - 1
        s = start - PET_len
        e = end + PET_len
        if s < 0:
            s = 0
        if e > len(record):
            e = len(record)
        record = record[s:e+1]

    if not SE:
        PET = record[:start]
        # add the end base to PET sequence
        add_base_to_PET(PET, rest[2])
        return PET
    else:
        PET1 = record[:start]
        PET2 = record[end+1:]
        PET2 = reverse_complement_record(PET2)
        add_base_to_PET(PET1, rest[2])
        add_base_to_PET(PET2, rest[2])
        return PET1, PET2


def align_linker(seq, linker, mismatch_threshold):
    """
    align linker within seq, use local alignment.
    """
    err_rate = float(mismatch_threshold)/len(linker)
    aligner = Aligner(seq, err_rate)
    alignment = aligner.locate(linker)
    if not alignment:
        # linker can't alignment to seq
        return False
    start, end, s_, e_, m, err = alignment
    return (start, end)


def regex_search_linker(seq, linker, mismatch_threshold):
    """
    match linker within seq, use regular expression(regex) do fuzzy match
    """
    reg = get_linker_regex(linker, mismatch_threshold)
    m = reg.search(seq)
    if not m:
        return False
    return m.span()


def match_linker(seq, linker, mismatch_threshold, allow_gap, seed_ratio=0.25):
    """ linker match algorithm """
    seed_len = int(floor(seed_ratio * len(linker)))
    seed = linker[:seed_len]
    if seed not in seq:
        return False

    if mismatch_threshold == 0:
        start = seq.find(linker)
        end = start + len(linker)
        span = (start, end)
        return span

    if allow_gap:
        span = align_linker(seq, linker, mismatch_threshold)
    else:
        span = regex_search_linker(seq, linker, mismatch_threshold)
    return span


def worker(task_queue, output_queue, counter_queue, linkers, mismatch, allow_gap, rest):
    """ stream processing(PET extract) task """
    from Queue import Empty
    from time import time 
    while 1:
        all, inter, intra, unmatch = 0,0,0,0 # variables for count reads
        try:
            records = task_queue.get(timeout=TIME_OUT)
        except Empty:
            break
        PETs = []
        for r in records:
            all += 1
            seq = str(r.seq)
            for ltype, linker in linkers.items():
                span = match_linker(seq, linker, mismatch, allow_gap)
                if span:
                    # linker matched
                    if   (ltype == 'A-A') or (ltype == 'B-B'):
                        # intra-molcular interaction
                        intra += 1
                        PET = extract_PET(r, span, rest)
                        PETs.append(PET)
                    elif (ltype == 'A-B') or (ltype == 'B-A'):
                        # inter-molcular interaction
                        inter += 1
                    break
            else:
                # all linkers can't match
                unmatch += 1
        output_queue.put(PETs)
        counter_queue.put((all, inter, intra, unmatch))


def worker_SE(task_queue, output_queue_1, output_queue_2, counter_queue, linkers, mismatch, allow_gap, rest, PET_len):
    """ stream processing(PET extract) task for SE mode """
    from Queue import Empty
    from time import time 
    while 1:
        all, inter, intra, unmatch = 0,0,0,0 # variables for count reads
        try:
            records = task_queue.get(timeout=TIME_OUT)
        except Empty:
            break
        PETs_1 = []
        PETs_2 = []
        for r in records:
            all += 1
            seq = str(r.seq)
            for ltype, linker in linkers.items():
                span = match_linker(seq, linker, mismatch, allow_gap)
                if span:
                    # linker matched
                    if   (ltype == 'A-A') or (ltype == 'B-B'):
                        # intra-molcular interaction
                        intra += 1
                        PET_1, PET_2 = extract_PET(r, span, rest, PET_len, SE=True)
                        PETs_1.append(PET_1)
                        PETs_2.append(PET_2)
                    elif (ltype == 'A-B') or (ltype == 'B-A'):
                        # inter-molcular interaction
                        inter += 1
                    break
            else:
                # all linkers can't match
                unmatch += 1
        output_queue_1.put(PETs_1)
        output_queue_2.put(PETs_2)
        counter_queue.put((all, inter, intra, unmatch))


def stream_processing(fastq_iter, fastq_writer, processes, linkers, mismatch, allow_gap, rest_site):
    """ assign stream processing tasks using multiprocessing. """
    task_queue = Queue()
    output_queue = Queue()

    # a global queue for count record how many:
    #   intra-molecular(A-A B-B)
    #   inter-molecular(A-B B-A)
    #   unmatchable
    # reads were captured
    # element in queue: (all, intra-molecular, inter-molecular, unmatchable)
    counter_queue = Queue()
    counter = {
        'all' : 0,
        'inter-molecular' : 0,
        'intra-molecular' : 0,
        'unmatchable': 0,
    }

    workers = [Process(target=worker, 
                         args=(task_queue, output_queue,
                               counter_queue, linkers, mismatch, allow_gap, rest_site))
               for i in range(processes)]
    output_p = Process(target=output, 
            args=(fastq_writer, output_queue))

    for w in workers:
        w.start()
    
    fastq_writer.write_header()
    output_p.start()

    try:
        while 1:
            task_queue.put(read_fastq(fastq_iter))
    except StopIteration:
        pass

    for w in workers:
        w.join()

    output_p.join()

    while not counter_queue.empty():
        all, inter, intra, un = counter_queue.get()
        counter['all'] += all
        counter['inter-molecular'] += inter
        counter['intra-molecular'] += intra
        counter['unmatchable'] += un

    return counter


def stream_processing_SE(fastq_iter, fastq_writer_1, fastq_writer_2,
        processes, linkers, mismatch, allow_gap, rest_site, PET_len):
    """ assign stream processing tasks using multiprocessing. for SE mode """
    task_queue = Queue()
    output_queue_1 = Queue()
    output_queue_2 = Queue()

    # a global queue for count record how many:
    #   intra-molecular(A-A B-B)
    #   inter-molecular(A-B B-A)
    #   unmatchable
    # reads were captured
    # element in queue: (all, intra-molecular, inter-molecular, unmatchable)
    counter_queue = Queue()
    counter = {
        'all' : 0,
        'inter-molecular' : 0,
        'intra-molecular' : 0,
        'unmatchable': 0,
    }

    workers = [Process(target=worker_SE, 
                         args=(task_queue, output_queue_1, output_queue_2,
                               counter_queue, linkers, mismatch, allow_gap, rest_site, PET_len))
               for i in range(processes)]
    output_p_1 = Process(target=output, 
            args=(fastq_writer_1, output_queue_1))
    output_p_2 = Process(target=output, 
            args=(fastq_writer_2, output_queue_2))

    for w in workers:
        w.start()
    
    fastq_writer_1.write_header()
    output_p_1.start()
    fastq_writer_2.write_header()
    output_p_2.start()

    try:
        while 1:
            task_queue.put(read_fastq(fastq_iter))
    except StopIteration:
        pass

    for w in workers:
        w.join()

    output_p_1.join()
    output_p_2.join()

    while not counter_queue.empty():
        all, inter, intra, un = counter_queue.get()
        counter['all'] += all
        counter['inter-molecular'] += inter
        counter['intra-molecular'] += intra
        counter['unmatchable'] += un

    return counter


def fastq_iter(file_in, phred):
    """ return a fastq iterator """
    if phred == 33:
        fastq_iter = SeqIO.parse(file_in, 'fastq')
    else:
        fastq_iter = SeqIO.parse(file_in, 'fastq-illumina')
    return fastq_iter


def fastq_writer(file_out, phred):
    """ return a fastq writer """
    if phred == 33:
        fastq_writer = SeqIO.QualityIO.FastqPhredWriter(file_out)
    else:
        fastq_writer = SeqIO.QualityIO.FastqIlluminaWriter(file_out)
    return fastq_writer


def mainPE(input1, input2, out1, out2,
        linekr_a, linker_b,
        mismatch, allow_gap, rest, phred, processes):
    # parse restriction enzyme site
    rest_site = parse_rest(rest)

    # load linkers
    linkers = load_linkers(linekr_a, linker_b)
    log_linkers(linkers)

    with open_file(input1) as file_in1, open_file(input2) as file_in2,\
         open_file(out1, 'w') as file_out1, open_file(out2, 'w') as file_out2:

        fastq_iter_1 = fastq_iter(file_in1, phred)
        fastq_iter_2 = fastq_iter(file_in2, phred)
        fastq_writer_1 = fastq_writer(file_out1, phred)
        fastq_writer_2 = fastq_writer(file_out2, phred)

        counts_1 = stream_processing(fastq_iter_1, fastq_writer_1,
            processes, linkers, mismatch, allow_gap, rest_site)

        counts_2 = stream_processing(fastq_iter_2, fastq_writer_2,
            processes, linkers, mismatch, allow_gap, rest_site)
    
    return (counts_1, counts_2)


def mainSE(input, out1, out2,
        linekr_a, linker_b,
        mismatch, allow_gap, rest, phred, processes, PET_len):
    # parse restriction enzyme site
    rest_site = parse_rest(rest)

    # load linkers
    linkers = load_linkers(linekr_a, linker_b)
    log_linkers(linkers)

    with open_file(input) as file_in1, open_file(input) as file_in2,\
         open_file(out1, "w") as file_out1, open_file(out2, 'w') as file_out2:

        fq_iter = fastq_iter(file_in1, phred)
        fastq_writer_1 = fastq_writer(file_out1, phred)
        fastq_writer_2 = fastq_writer(file_out2, phred)

        counts = stream_processing_SE(fq_iter, fastq_writer_1, fastq_writer_2,
            processes, linkers, mismatch, allow_gap, rest_site, PET_len)
    
    return counts


if __name__ == "__main__":
    parser = argument_parser()
    args = parser.parse_args()

    read_args(args, globals())
    
    if command == 'SE':
        counts = mainSE(input, out1, out2,
                linker_a, linker_b,
                mismatch, allow_gap, rest, phred, processes, PET_len)
        log_counts(counts)

    elif command == 'PE':
        counts_1, counts_2 = mainPE(input1, input2, out1, out2,
                linker_a, linker_b,
                mismatch, allow_gap, rest, phred, processes)
        log_counts(counts_1)
        log_counts(counts_2)