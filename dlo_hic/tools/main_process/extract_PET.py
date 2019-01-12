import re
import logging
from itertools import repeat
from concurrent.futures import ProcessPoolExecutor
from collections import OrderedDict

import click
import numpy as np

from dlo_hic.utils import reverse_complement as rc
from dlo_hic.utils.fastqio import read_fastq, write_fastq
from dlo_hic.utils.filetools import open_file
from dlo_hic.utils.linker_trim import process_chunk, COUNT_ITEM_NAMES


log = logging.getLogger(__name__)


def parse_rest(rest_str):
    """
    parse restriction enzyme site sequence
    """
    left, nick, right = re.split("[*^]", rest_str)
    return left, nick, right


def load_linkers(linker_a, linker_b):
    """
    compose linkers

    Return
    ------
    linkers : tuple
        (AA, BB, AB, BA)
    """
    linkers = {}
    # convert linker sequence from unicode to str
    linker_a = str(linker_a)
    linker_b = str(linker_b) if linker_b else None
    if linker_b:
        AA = linker_a + rc(linker_a)
        BB = linker_b + rc(linker_b)
        AB = linker_a + rc(linker_b)
        BA = linker_b + rc(linker_a)
        linkers = (AA, BB, AB, BA)
    else:
        AA = linker_a + rc(linker_a)
        linkers = (AA, None, None, None)
    return linkers


def log_linkers(linkers):
    log.info("linkers:")
    linkers_name = ("AA", "BB", "AB", "BA")
    for name, linker in zip(linkers_name, linkers):
        if linker is None:
            break
        log.info("\t{}:\t{}".format(name, linker))


def log_counts(counts, log_file=None):
    if not log_file:
        log.info("Quality Control:")
        for k, v in zip(COUNT_ITEM_NAMES, counts):
            log.info("\t{}:\t{}".format(k, v))
        try:
            ratio = counts[-2] / float(counts[-1])
        except ZeroDivisionError:
            ratio = 0
        log.info("Valid reads ratio: {}".format(ratio))
    else:
        with open(log_file, 'w') as f:
            for k, v in zip(COUNT_ITEM_NAMES, counts):
                outline = "\t".join([str(k), str(v)]) + "\n"
                f.write(outline)


def chunking(fq_iter, chunk_size=10000):
    chunk = []
    for idx, fq_rec in enumerate(fq_iter):
        if idx != 0 and idx % chunk_size == 0:
            yield chunk
            chunk = []
        chunk.append(fq_rec)
    yield chunk


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
    default="",
    help="If specified, Cut the adapter sequence in the PET2.")
@click.option("--mismatch-adapter", "mismatch_adapter", default=3,
    help="mismatch threshold in alignment in cut adapter step.")
@click.option("--log-file",
    default="PET_count.txt",
    help="Sperate log file record reads count information. default PET_count.txt")
@click.option("--flag-file",
    help="If specified, write the addition information(like flag) in linker trim to a file.")
@click.option("--chunk-size",
    default=10000,
    show_default=True,
    help="How many records in one chunk, when do parallel computing.")
def _main(fastq, out1, out2,
          linker_a, linker_b,
          mismatch, rest, processes, PET_len_range, PET_cut_len,
          adapter, mismatch_adapter, log_file, flag_file, chunk_size):
    """
    Extract the PETs sequences on both sides of linker sequence.

    \b
    Input:
        fastq file, support gziped file.

    \b
    Output:
        Two fastq files contain PETs sequences.

    \b
    About flag:
        Bit    Description
        ------------------
        1      linker unmatchable
        2      intra-molecular linker (AA, BB)
        4      adapter unmatchable
        8      added base to left PET
        16     added base to right PET
        32     left PET len less than threshold
        64     left PET len large than threshold
        128    right PET len less than threshold
        256    right PET len large than threshold

    """
    log.info("Extract PETs from file %s"%fastq)

    # parse restriction enzyme site
    rest_site = parse_rest(rest)
    log.info("enzyme cutting site: %s"%rest)

    # load linkers
    if linker_b == linker_a:
        linker_b = None
    linkers = load_linkers(linker_a, linker_b)
    log_linkers(linkers)
    log.info("linker alignment mismatch threshold: {}".format(mismatch))

    if adapter:
        log.info("adapter: {}".format(adapter))
        log.info("adapter alignment mismatch threshold: {}".format(mismatch_adapter))

    log.info("Expect PET length range: [{}, {}]".format(*PET_len_range))


    # Perform linker trimer on fastq file

    fq_iter = read_fastq(fastq)
    fq_pet1 = open_file(out1, mode='w')
    fq_pet2 = open_file(out2, mode='w')
    if flag_file:
        flag_fh = open_file(flag_file, "w")
        header = "Read_ID\tFlag\n"
        flag_fh.write(header)

    counts = np.zeros(len(COUNT_ITEM_NAMES), dtype=np.int)
    with ProcessPoolExecutor(max_workers=processes) as exc:
        map_ = exc.map if processes > 1 else map
        args = linkers, adapter, rest_site, mismatch, mismatch_adapter, PET_len_range, PET_cut_len
        for out_chunk, counts_ in map_(process_chunk, chunking(fq_iter, chunk_size), repeat(args)):
            counts += counts_
            for fq_rec, flag, PET1, PET2 in out_chunk:
                write_fastq(PET1, fq_pet1)
                write_fastq(PET2, fq_pet2)
                if flag_file:
                    out_line = "{}\t{}\n".format(fq_rec.seqid, flag)
                    flag_fh.write(out_line)

    fq_pet1.close()
    fq_pet2.close()
    if flag_file:
        flag_fh.close()

    # log counts
    log_counts(counts)
    log_counts(counts, log_file)

    return counts


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
