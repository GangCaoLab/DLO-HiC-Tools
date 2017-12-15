from __future__ import print_function
import re
import sys
import json
import gzip
import itertools
import subprocess
from math import floor
from copy import copy
import time
import multiprocessing
from multiprocessing import Process, Manager, Queue

import click
from Bio import SeqIO

from dlo_hic.utils.align import Aligner
from dlo_hic.utils import reverse_complement as rc
from dlo_hic.utils import read_args


TIME_OUT = 1
CHUNK_SIZE = 1000


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
    if linker_b:
        linkers['A-A'] = linker_a + rc(linker_a)
        linkers['A-B'] = linker_a + rc(linker_b)
        linkers['B-A'] = linker_b + rc(linker_a)
        linkers['B-B'] = linker_b + rc(linker_b)
    else:
        linkers['A-A'] = linker_a + rc(linker_a)
    return linkers


def log_linkers(linkers, file=sys.stderr):
    print("linkers:", file=file)
    for key, linker in linkers.items():
        print("{}\t{}".format(key, linker), file=file)


def log_counts(counts, file=sys.stderr):
    print("Quality Control:", file=file)
    for k, v in counts.items():
        print("\t{}\t{}".format(k, v), file=file)
    try:
        ratio = counts['intra-molecular'] / float(counts['all'])
    except ZeroDivisionError:
        ratio = 0
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
    default chunk size 1000
    """
    chunk_iter = itertools.islice(fastq_iter, n)
    chunk = list(chunk_iter)
    if len(chunk) == 0:
        raise StopIteration
    return chunk


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
    return (start, end-1)


def match_linker(seq, linker, mismatch_threshold, seed_ratio=0.25):
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

    span = align_linker(seq, linker, mismatch_threshold)
    return span


def worker_SE(task_queue, out1, out2, phred, counter, linkers, mismatch, rest, PET_len):
    """ stream processing(PET extract) task for SE mode """
    from Queue import Empty
    current = multiprocessing.current_process().pid

    file_out1 = open_file(out1, 'w')
    file_out2 = open_file(out2, 'w')
    fastq_writer_1 = fastq_writer(file_out1, phred)
    fastq_writer_2 = fastq_writer(file_out2, phred)

    fastq_writer_1.write_header()
    fastq_writer_2.write_header()

    all, inter, intra, unmatch = 0,0,0,0 # variables for count reads
    while 1:
        try:
            records = task_queue.get(timeout=TIME_OUT)
            if records is None:
                raise Empty
        except Empty:
            fastq_writer_1.write_footer()
            fastq_writer_1.handle.close()
            fastq_writer_2.write_footer()
            fastq_writer_2.handle.close()
            # update counter dict
            counter['all'] += all
            counter['inter-molecular'] += inter
            counter['intra-molecular'] += intra
            counter['unmatchable'] += unmatch
            # log
            print("Process-%d"%current , "exit.", file=sys.stderr)
            break

        for r in records:
            all += 1
            seq = str(r.seq)
            for ltype, linker in linkers.items():
                span = match_linker(seq, linker, mismatch)
                if span:
                    # linker matched
                    if   (ltype == 'A-A') or (ltype == 'B-B'):
                        # intra-molcular interaction
                        intra += 1
                        PET_1, PET_2 = extract_PET(r, span, rest, PET_len, SE=True)
                        fastq_writer_1.write_record(PET_1)
                        fastq_writer_2.write_record(PET_2)
                    elif (ltype == 'A-B') or (ltype == 'B-A'):
                        # inter-molcular interaction
                        inter += 1
                    break
            else:
                # all linkers can't match
                unmatch += 1


def fastq_iter(file_in, phred):
    """ return a fastq iterator """
    if str(phred) == '33':
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


@click.command(name="extract_PET")
@click.argument("fastq", nargs=1)
@click.option("--out1", '-o1', required=True,
    help="output1: right side PET fastq file")
@click.option("--out2", '-o2', required=True,
    help="output2: right side PET fastq file")
@click.option("--linker-A", required=True,
    help="The sequence of linkerA")
@click.option("--linker-B",
    help="The sequence of linkerB")
@click.option("--mismatch", default=4,
    help="threshold of linkers base mismatch(and gap open extends) number, default 3")
@click.option("--rest", default="A*AGCT*T",
    help="The sequence of restriction enzyme recognition site, " +\
         "default HindIII: 'A*AGCT*T' ")
@click.option("--phred", default='33', type=click.Choice(['33', '64']),
    help="The Phred score encode offset type, 33 or 64. default 33")
@click.option("--processes", "-p", default=1,
    help="Use how many processes do calculation. default 1")
@click.option("--PET-len", default=0, 
    help="The expected length of PET sequence," +\
         "if PET_len==0 (default) will not limit length.")
def _main(fastq, out1, out2,
        linekr_a, linker_b,
        mismatch, rest, phred, processes, PET_len):
    """
    Extract the PETs sequences on both sides of linker sequence.

    This script accept a fastq or fasta file,
    and output two PETs files in fastq file format.

    """
    # parse restriction enzyme site
    rest_site = parse_rest(rest)

    # load linkers
    linkers = load_linkers(linekr_a, linker_b)
    log_linkers(linkers)

    manager = Manager()
    task_queue = Queue()
    counter = manager.dict() # a global queue for count record how many:
    counter['all'] = 0
    counter['inter-molecular'] = 0
    counter['intra-molecular'] = 0
    counter['unmatchable'] = 0

    workers = [Process(target=worker_SE, 
                         args=(task_queue, out1+".tmp.%d"%i, out2+".tmp.%d"%i, phred,
                               counter, linkers, mismatch, rest_site, PET_len))
               for i in range(processes)]

    for w in workers:
        w.start()

    # put reads in queue
    with open_file(fastq) as file_in:
        fq_iter = fastq_iter(file_in, phred)
        try:
            while 1:
                task_queue.put(read_fastq(fq_iter))
        except StopIteration:
            pass

    # wait subprocesses end
    for w in workers:
        w.join()

    # merge all tmp files
    tmpfiles_1 = [out1+".tmp.%d"%i for i in range(processes)]
    tmpfiles_2 = [out2+".tmp.%d"%i for i in range(processes)]
    cmd = "cat " + " ".join(tmpfiles_1) + " > " + out1
    subprocess.check_call(cmd, shell=True)
    cmd = "cat " + " ".join(tmpfiles_2) + " > " + out2
    subprocess.check_call(cmd, shell=True)
    cmd = "rm " + " ".join(tmpfiles_1 + tmpfiles_2)
    subprocess.check_call(cmd, shell=True)

    counts = dict(counter)

    return counts


main = _main.callback


if __name__ == "__main__":
    counts = _main()
    log_counts(counts)
