import re
import logging
import itertools
import subprocess
from copy import copy
import multiprocessing as mp

import click

from dlo_hic.utils.align import Aligner
from dlo_hic.utils import reverse_complement as rc
from dlo_hic.utils import guess_fq_phred
from dlo_hic.utils.filetools import open_file
from dlo_hic.utils.fastqio import read_fastq, write_fastq, Fastq

CHUNK_SIZE = 1000

log = logging.getLogger(__name__)


def parse_rest(rest_str):
    """
    parse restriction enzyme site sequence
    """
    left, nick, right = re.split("[*^]", rest_str)
    return left, nick, right


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


def read_fastq_chunk(fastq_iter, n=CHUNK_SIZE):
    """
    read fastq file, return records chunk
    default chunk size 1000
    """
    chunk_iter = itertools.islice(fastq_iter, n)
    chunk = list(chunk_iter)
    if len(chunk) == 0:
        raise StopIteration
    return chunk


def add_base_to_PET(PET, base, pos='end'):
    """ add one base to fastq record object """
    if pos == 'end':
        qual = PET.quality[-1]
        new = Fastq(PET.seqid, PET.seq + base, PET.quality + qual)
    else:
        qual = PET.quality[0]
        new = Fastq(PET.seqid, base + PET.seq, qual + PET.quality)
    return new


def extract_PET(record, span, rest, adapter=("", 0)):
    """
    extract PET in matched record.
    """
    start, end = span
    PET1 = record[:start]
    PET2 = record[end+1:]

    adapter, mismatch = adapter
    if adapter:
        PET2 = cut_adapter(PET2, adapter, mismatch_threshold=mismatch)

    pet1_end = rest[0] + rest[1]
    pet2_start = rest[1] + rest[2]
    if PET1.seq[-len(pet1_end):] == pet1_end:
        PET1 = add_base_to_PET(PET1, rest[2], pos='end')
    if PET2.seq[:len(pet2_start)] == pet2_start:
        PET2 = add_base_to_PET(PET2, rest[0], pos='start')
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


def cut_adapter(rec, adapter_pattern, mismatch_threshold):
    """ cut the adapter sequence """
    matched = match_(rec.seq, adapter_pattern, mismatch_threshold)
    if not matched:
        clean_seq = rec
    else:
        start, end = matched
        clean_seq = rec[:start]
    return clean_seq


def cut_PET(PET1, PET2, length_range, PET_cut_len):
    lower, upper = length_range
    if len(PET1) < lower:
        PET1 = False
    elif len(PET1) > upper:
        PET1 = PET1[-PET_cut_len:]

    if len(PET2) < lower:
        PET2 = False
    elif len(PET2) > upper:
        PET2 = PET2[:PET_cut_len]

    return PET1, PET2


def worker(task_queue, counter, lock, # objects for multi process work
           out1, out2, # output file name
           linkers, mismatch, rest, PET_len_range, PET_cut_len, adapter):  # parameters
    """ stream processing(PET extract) task """
    from queue import Empty
    current = mp.current_process().pid

    file_out1 = open_file(out1, 'w')
    file_out2 = open_file(out2, 'w')

    all, inter, intra, unmatch = 0,0,0,0  # variables for count reads
    while 1:
        records = task_queue.get()
        if records is None:

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
                if span:  # linker matched
                    if (ltype == 'A-A') or (ltype == 'B-B'):
                        # intra-molcular interaction
                        PET_1, PET_2 = extract_PET(r, span, rest, adapter)
                        PET_1, PET_2 = cut_PET(PET_1, PET_2, PET_len_range, PET_cut_len)
                        if (not PET_1) or (not PET_2):
                            # PET too short
                            unmatch += 1
                            continue
                        intra += 1
                        write_fastq(PET_1, file_out1)
                        write_fastq(PET_2, file_out2)
                    elif (ltype == 'A-B') or (ltype == 'B-A'):
                        # inter-molcular interaction
                        inter += 1
                    break
            else: # all linkers can't match
                unmatch += 1


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
@click.option("--processes", "-p", default=1,
    help="Use how many processes do calculation. default 1")
@click.option("--PET-len-range", 'PET_len_range',
    default=(10, 22), show_default=True,
    help="The expected length range of PET sequence," +\
         "if the PET_length exceed the upper bound will cut the exceeded sequence, " +\
         "if it lower than the lower bound will treat the sequence as the unmatched sequence. ")
@click.option("--PET-cut-len", 'PET_cut_len',
    default=20, show_default=True,
    help="If PET length large than the upper len range, will cut to this length.")
@click.option("--cut-adapter", "adapter",
    help="If specified, Cut the adapter sequence in the PET2.")
@click.option("--mismatch-adapter", "mismatch_adapter", default=3,
    help="mismatch threshold in alignment in cut adapter step.")
@click.option("--log-file",
    default="PET_count.txt",
    help="Sperate log file record reads count information. default PET_count.txt")
def _main(fastq, out1, out2,
        linker_a, linker_b,
        mismatch, rest, processes, PET_len_range, PET_cut_len,
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
                          args=(task_queue, counter, lock, out1+".tmp.%d"%i, out2+".tmp.%d"%i,
                                linkers, mismatch, rest_site, PET_len_range, PET_cut_len, (adapter, mismatch_adapter)))
               for i in range(processes)]

    for w in workers:
        w.start()

    log.info("%d worker process spawned for extract PETs."%len(workers))

    # put reads in queue
    with open_file(fastq) as file_in:
        fq_iter = read_fastq(fastq)
        try:
            while 1:
                task_queue.put(read_fastq_chunk(fq_iter))
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
    eval("_main()")
