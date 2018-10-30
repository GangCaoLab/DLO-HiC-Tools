import re
import io
import sys
import json
import gzip
import logging
import itertools
import subprocess
from math import floor
from copy import copy
import time
import multiprocessing as mp

import click
from Bio import SeqIO

from dlo_hic.utils.align import Aligner
from dlo_hic.utils import reverse_complement as rc
from dlo_hic.utils import read_args

CHUNK_SIZE = 1000

log = logging.getLogger(__name__)


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
    # convert linker sequence from unicode to str
    linker_a = str(linker_a)
    linker_b = str(linker_b) if linker_b else None
    if linker_b:
        linkers['A-A'] = linker_a + rc(linker_a)
        linkers['A-B'] = linker_a + rc(linker_b)
        linkers['B-A'] = linker_b + rc(linker_a)
        linkers['B-B'] = linker_b + rc(linker_b)
    else:
        log.warning("using single linker.")
        linkers['A-A'] = linker_a + rc(linker_a)
    return linkers


def log_linkers(linkers):
    log.info("linkers:")
    for key, linker in linkers.items():
        log.info("{}\t{}".format(key, linker))


def log_counts(counts, log_file=None):
    if not log_file:
        log.info("Quality Control:")
        for k, v in counts.items():
            log.info("\t{}\t{}".format(k, v))
        try:
            ratio = counts['intra-molecular'] / float(counts['all'])
        except ZeroDivisionError:
            ratio = 0
        log.info("Valid reads ratio: {}".format(ratio))
    else:
        with open(log_file, 'w') as f:
            for k, v in counts.items():
                outline = "\t".join([str(k), str(v)]) + "\n"
                f.write(outline)


def open_file(fname, mode='r'):
    if fname.endswith('.gz'):
        fh = gzip.open(fname, mode=mode+'b')
        fh = io.TextIOWrapper(fh)
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


def add_base_to_PET(PET, base, base_qual=38):
    """ add one base to fastq record object """
    quality = PET.letter_annotations['phred_quality']
    quality.append(base_qual)
    PET.letter_annotations = {}
    PET.seq = PET.seq + base
    PET.letter_annotations['phred_quality'] = quality


def extract_PET(record, span, rest,
                PET_len=None, adapter=("", 0)):
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

    PET1 = record[:start]
    PET2 = record[end+1:]

    adapter, mismatch = adapter
    if adapter:
        PET2 = cut_adapter(PET2, adapter, mismatch_threshold=mismatch)

    PET2 = reverse_complement_record(PET2)
    add_base_to_PET(PET1, rest[2])
    add_base_to_PET(PET2, rest[2])
    return PET1, PET2


def align_(seq, pattern, mismatch_threshold):
    """
    align pattern within seq, use local alignment.
    if not matched return False.
    """
    err_rate = float(mismatch_threshold)/len(pattern)
    aligner = Aligner(seq, err_rate)
    aligner.min_overlap = len(pattern) - mismatch_threshold
    alignment = aligner.locate(pattern)
    if not alignment:
        # linker can't alignment to seq
        return False
    start, end, s_, e_, m, err = alignment
    return (start, end-1)


def match_(seq, pattern, mismatch_threshold):
    """
    match algorithm.
    if not matched return False.
    """
    if mismatch_threshold == 0:
        start = seq.find(pattern)
        if start == -1:
            return False
        end = start + len(pattern)
        span = (start, end)
        return span

    span = align_(seq, pattern, mismatch_threshold)
    return span


def cut_adapter(seq, adapter_pattern, mismatch_threshold):
    """ cut the adapter sequence """
    matched = match_(str(seq.seq), adapter_pattern, mismatch_threshold)
    if not matched:
        clean_seq = seq
    else:
        start, end = matched
        clean_seq = seq[:start]
    return clean_seq


def worker(task_queue, counter, lock, # objects for multi process work
           out1, out2, # output file name
           phred, linkers, mismatch, rest, PET_len, adapter):  # parameters
    """ stream processing(PET extract) task """
    from queue import Empty
    current = mp.current_process().pid

    file_out1 = open_file(out1, 'w')
    file_out2 = open_file(out2, 'w')
    fastq_writer_1 = fastq_writer(file_out1, phred)
    fastq_writer_2 = fastq_writer(file_out2, phred)

    fastq_writer_1.write_header()
    fastq_writer_2.write_header()

    all, inter, intra, unmatch = 0,0,0,0  # variables for count reads
    while 1:
        records = task_queue.get()
        if records is None:
            fastq_writer_1.write_footer()
            fastq_writer_1.handle.close()
            fastq_writer_2.write_footer()
            fastq_writer_2.handle.close()

            # update counter dict
            lock.acquire()
            counter['all'] += all
            counter['inter-molecular'] += inter
            counter['intra-molecular'] += intra
            counter['unmatchable'] += unmatch
            log.debug("conter increse, %d, %d, %d, %d"%(all, inter, intra, unmatch))
            log.debug("%d, %d, %d, %d"%(
                counter['all'],
                counter['inter-molecular'],
                counter['intra-molecular'],
                counter['unmatchable'],
            ))
            lock.release()

            log.debug("Process-%d done"%current)
            break

        for r in records:
            all += 1
            seq = str(r.seq) # extract reocrd's sequence
            for ltype, linker in linkers.items():
                span = match_(seq, linker, mismatch)
                if span: # linker matched
                    if   (ltype == 'A-A') or (ltype == 'B-B'):
                        # intra-molcular interaction
                        intra += 1
                        PET_1, PET_2 = extract_PET(r, span, rest, PET_len, adapter)
                        fastq_writer_1.write_record(PET_1)
                        fastq_writer_2.write_record(PET_2)
                    elif (ltype == 'A-B') or (ltype == 'B-A'):
                        # inter-molcular interaction
                        inter += 1
                    break
            else: # all linkers can't match
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
    if str(phred) == '33':
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
    help="threshold of linkers base mismatch(and gap open extends) number, default 4")
@click.option("--rest", default="A*AGCT*T",
    help="The sequence of restriction enzyme recognition site, " +\
         "default HindIII: 'A*AGCT*T' ")
@click.option("--phred", default='33', type=click.Choice(['33', '64']),
    help="The Phred score encode offset type, 33 or 64. default 33")
@click.option("--processes", "-p", default=1,
    help="Use how many processes do calculation. default 1")
@click.option("--PET-len", 'PET_len', default=0, 
    help="The expected length of PET sequence," +\
         "if 0 (default) will not limit length," +\
         "NOTE: Adapter for sequencing must be trimmed in this situation.")
@click.option("--cut-adapter", "adapter",
    help="If specified, Cut the adapter sequence in the PET2.")
@click.option("--mismatch-adapter", "mismatch_adapter", default=3,
    help="mismatch threshold in alignment in cut adapter step.")
@click.option("--log-file",
    default="PET_count.txt",
    help="Sperate log file record reads count information. default PET_count.txt")
def _main(fastq, out1, out2,
        linker_a, linker_b,
        mismatch, rest, phred, processes, PET_len,
        adapter, mismatch_adapter, log_file):
    """
    Extract the PETs sequences on both sides of linker sequence.

    Input:
        fastq file, support gziped file.

    Output:
        two fastq files contain PETs sequences.

    """
    log.info("Extract PETs from file %s"%fastq)

    # parse restriction enzyme site
    rest_site = parse_rest(rest)

    log.info("enzyme cutting site: %s"%rest)
    # load linkers
    linkers = load_linkers(linker_a, linker_b)
    log_linkers(linkers)

    manager = mp.Manager()
    task_queue = mp.Queue()
    counter = manager.dict() # a global queue for count record how many:
    lock = mp.Lock()
    # init counter
    counter['all'] = 0
    counter['inter-molecular'] = 0
    counter['intra-molecular'] = 0
    counter['unmatchable'] = 0

    workers = [mp.Process(target=worker, 
                          args=(task_queue, counter, lock, out1+".tmp.%d"%i, out2+".tmp.%d"%i, phred,
                                linkers, mismatch, rest_site, PET_len, (adapter, mismatch_adapter)))
               for i in range(processes)]

    for w in workers:
        w.start()

    log.info("%d worker process spawned for extract PETs."%len(workers))

    # put reads in queue
    with open_file(fastq) as file_in:
        fq_iter = fastq_iter(file_in, phred)
        try:
            while 1:
                task_queue.put(read_fastq(fq_iter))
        except StopIteration:
            pass
    
    for w in workers:
        task_queue.put(None)

    # wait subprocesses end
    for w in workers:
        w.join()

    # merge all tmp files
    tmpfiles_1 = [out1+".tmp.%d"%i for i in range(processes)]
    tmpfiles_2 = [out2+".tmp.%d"%i for i in range(processes)]
    # merge tmp files
    log.info("merging temporary files ...")
    cmd = "cat " + " ".join(tmpfiles_1) + " > " + out1
    subprocess.check_call(cmd, shell=True)
    cmd = "cat " + " ".join(tmpfiles_2) + " > " + out2
    subprocess.check_call(cmd, shell=True)
    # delete tmp files
    log.info("delete temporary files.")
    cmd = "rm " + " ".join(tmpfiles_1 + tmpfiles_2)
    subprocess.check_call(cmd, shell=True)

    counts = dict(counter)
    log_counts(counts)

    # write counts info to sperate file
    log_counts(counts, log_file)

    return counts


main = _main.callback


if __name__ == "__main__":
    main()
