import sys
import re
import logging
import multiprocessing as mp
from queue import Empty

import click

from dlo_hic.utils import reverse_complement as rc
from dlo_hic.utils.parse_text import parse_rest
from dlo_hic.utils.fastqio import read_fastq, write_fastq
from dlo_hic.utils.filetools import open_file, merge_tmp_files
from dlo_hic.utils.linker_trim import process_chunk, COUNT_ITEM_NAMES, init_counts
from dlo_hic.utils.infer_adapter import infer_adapter_seq


log = logging.getLogger(__name__)




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
        for k, v in zip(COUNT_ITEM_NAMES, counts['flag']):
            log.info("\t{}:\t{}".format(k, v))
        try:
            ratio = counts['flag'][-2] / float(counts['flag'][-1])
        except ZeroDivisionError:
            ratio = 0
        log.info("Valid reads ratio: {}".format(ratio))
    else:
        with open(log_file, 'w') as f:
            f.write("# Flag counts(quality control), columns: [item, number]\n")
            for k, v in zip(COUNT_ITEM_NAMES, counts['flag']):
                outline = "\t".join([str(k), str(v)]) + "\n"
                f.write(outline)

            f.write("\n")
            f.write("# PET1 length distribution, columns: [length, number]\n")
            for k, v in counts['pet1_len'].most_common():
                f.write("{}\t{}\n".format(k, v))
            f.write("\n")
            f.write("# PET2 length distribution, columns: [length, number]\n")
            for k, v in counts['pet2_len'].most_common():
                f.write("{}\t{}\n".format(k, v))

            f.write("\n")
            f.write("# Linker match score distribution, columns: [score, number]\n")
            for k, v in counts['linker_match_score'].most_common():
                f.write("{}\t{}\n".format(k, v))

            f.write("\n")
            f.write("# Adapter match score distribution, columns: [score, number]\n")
            for k, v in counts['adapter_match_score'].most_common():
                f.write("{}\t{}\n".format(k, v))


def chunking(fq_iter, chunk_size=10000):
    chunk = []
    for idx, fq_rec in enumerate(fq_iter):
        if idx != 0 and idx % chunk_size == 0:
            yield chunk
            chunk = []
        chunk.append(fq_rec)
    yield chunk


def worker(task_queue, out1, out2, flag_file, lock, counter, args):
    current = mp.current_process().pid
    fq_pet1 = open_file(out1, mode='w')
    fq_pet2 = open_file(out2, mode='w')
    if flag_file:
        flag_fh = open_file(flag_file, "w")

    counts = init_counts()

    while 1:
        try:
            chunk = task_queue.get()
            if chunk is None:
                raise Empty
        except Empty:
            fq_pet1.close()
            fq_pet2.close()
            if flag_file:
                flag_fh.close()
            log.debug("Process-%d done."%current)
            lock.acquire()
            counter[current] = counts
            lock.release()
            break

        out_chunk, counts_ = process_chunk(chunk, args)
        for cnt_itms in counts.keys():
            counts[cnt_itms] += counts_[cnt_itms]

        for fq_rec, flag, PET1, PET2, align, align_ada in out_chunk:
            if (flag & 1 == 0) and (flag & 32 == 0) and (flag & 128 == 0):
                write_fastq(PET1, fq_pet1)
                write_fastq(PET2, fq_pet2)
            if flag_file:  # write flag file
                out_itms = [fq_rec.seqid, flag]
                out_itms.append(fq_rec.seq)

                if align:
                    lin_s, lin_e, lin_q_s, lin_q_e, lin_match, lin_cost = align
                    out_itms.extend([lin_s, lin_e, lin_q_s, lin_q_e, lin_match, lin_cost])
                    p2_ = fq_rec.seq[lin_e:]
                else:
                    out_itms.extend(['.']*6)
                    p2_ = '.'

                if PET1:
                    out_itms.append(PET1.seq)
                    out_itms.append(p2_)
                else:
                    out_itms.append('.')
                    out_itms.append('.')

                if align_ada:
                    ada_s, ada_e, ada_q_s, ada_q_e, ada_match, ada_cost = align_ada
                    out_itms.extend([ada_s, ada_e, ada_q_s, ada_q_e, ada_match, ada_cost])
                else:
                    out_itms.extend(['.'] * 6)

                if PET2:
                    out_itms.append(PET2.seq)
                else:
                    out_itms.append('.')

                out_line = "\t".join([str(i) for i in out_itms]) + "\n"
                flag_fh.write(out_line)


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
@click.option("--mismatch", default=0.12,
    show_default=True,
    help="Threshold of linkers base mismatch(and gap open extends) number. "
         "If specify int type will fix the mismatch size, "
         "otherwise, specify float type, will assign the ratio mismatch/len(linker)")
@click.option("--rest", default="A^AGCTT",
    show_default=True,
    help="The sequence of restriction enzyme recognition site. Use '^' to indicate the cutting site.")
@click.option("--processes", "-p", default=1,
    show_default=True,
    help="Use how many processes do calculation. default 1")
@click.option("--PET-len-range", 'PET_len_range',
    default=(10, 22), show_default=True,
    help="The expected length range of PET sequence, "
         "if the PET_length exceed the upper bound will cut the exceeded sequence, "
         "if it lower than the lower bound will treat the sequence as the unmatched sequence. ")
@click.option("--PET-cut-len", 'PET_cut_len',
    default=20, show_default=True,
    help="If PET length large than the upper len range, will cut to this length.")
@click.option("--cut-adapter", "adapter",
    default="auto",
    help="Specify the sequence of adapter. If specified, Cut the adapter sequence in the PET2."
         "If specify 'auto', will inference adapter sequence using mafft.")
@click.option("--mismatch-adapter", "mismatch_adapter", default=0.3,
    show_default=True,
    help="Mismatch threshold in alignment in cut adapter step. If specify int type will fix the mismatch size, "
         "otherwise, specify float type, will assign the ratio mismatch/len(adapter)")
@click.option("--log-file",
    default="PET_count.txt",
    show_default=True,
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
    if mismatch.is_integer():
        mismatch = int(mismatch)
    else:
        mismatch = int(len(linkers[0]) * mismatch)
    log.info("linker alignment mismatch threshold: {}, ratio: {:.2f}".format(mismatch, mismatch/len(linkers[0])))

    if adapter:
        if adapter == 'auto':
            try:
                adapter = infer_adapter_seq(fastq)
                if adapter is None:
                    raise ValueError()
            except:
                log.error("Fail to inference adapter, please specify the adapter sequence directly.")
                sys.exit()
            adapter = adapter.upper()
            log.info("adapter sequence found")
        if not mismatch_adapter.is_integer():
            mismatch_adapter = int(mismatch_adapter * len(adapter))
        else:
            mismatch_adapter = int(mismatch_adapter)

        msg = "adapter: {}, length: {}, ".format(adapter, len(adapter))
        msg += "mismatch-N: {}, mismatch-ratio: {:.2f}".format(mismatch_adapter, mismatch_adapter/len(adapter))
        log.info(msg)

    log.info("Expect PET length range: [{}, {}]".format(*PET_len_range))

    # Perform linker trimer on fastq file

    args = linkers, adapter, rest_site, mismatch, mismatch_adapter, PET_len_range, PET_cut_len
    tmp_files_o1 = [out1+".tmp."+str(i) for i in range(processes)]
    tmp_files_o2 = [out2+".tmp."+str(i) for i in range(processes)]
    if flag_file:
        tmp_files_flag = [flag_file+".tmp."+str(i) for i in range(processes)]
    else:
        tmp_files_flag = ["" for _ in range(processes)]

    lock = mp.Lock()
    counter = mp.Manager().dict()
    task_queue = mp.Queue(maxsize=processes)
    workers = [mp.Process(target=worker,
        args=(task_queue, o1, o2, fg, lock, counter, args))
        for o1, o2, fg in zip(tmp_files_o1, tmp_files_o2, tmp_files_flag)
    ]

    log.info("%d workers spawned for extract PETs"%len(workers))
    for w in workers:
        w.start()

    fq_iter = read_fastq(fastq)
    for chunk in chunking(fq_iter, chunk_size):
        task_queue.put(chunk, block=True)

    for w in workers:
        task_queue.put(None)

    for w in workers:
        w.join()

    merge_tmp_files(tmp_files_o1, out1)
    merge_tmp_files(tmp_files_o2, out2)
    if flag_file:
        header_itms = ['read_ID', 'flag', 'read_seq', 'linker_start']
        header_itms.extend(['linker_end', 'linker_query_start', 'linker_query_end', 'linker_match', 'linker_cost'])
        header_itms.extend(['PET1_seq', 'PET2_seq_with_adapter'])
        header_itms.extend(['adapter_start', 'adapter_end', 'adapter_query_start', 'adapter_query_end', 'adapter_match', 'adapter_cost', 'pet2_seq'])
        header = "\t".join(header_itms) + "\n"
        merge_tmp_files(tmp_files_flag, flag_file, header)

    counts = init_counts()
    for pid in counter.keys():
        for cnt_itm in counts.keys():
            counts[cnt_itm] += counter[pid][cnt_itm]

    # log counts
    log_counts(counts)
    log_counts(counts, log_file)

    return counts


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
