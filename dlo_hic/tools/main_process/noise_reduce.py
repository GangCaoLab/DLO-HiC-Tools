import io
import logging
import subprocess
from queue import Empty
import multiprocessing as mp
from collections import defaultdict

import click
import numpy as np
import h5py

from dlo_hic.utils.parse_text import parse_line_bed6, parse_line_bedpe, parse_rest
from dlo_hic.utils.filetools import merge_tmp_files, open_file


TIME_OUT = 1
CHUNK_SIZE = 10000

log = logging.getLogger(__name__)


def load_rest_sites(frag_file):
    """
    load restriction cutting sites(positions) into memory.
    """
    rest_sites = {}
    if frag_file.endswith(".gz") or frag_file.endswith(".bed"):
        with io.TextIOWrapper(open_file(frag_file)) as f:
            for idx, line in enumerate(f):
                line = line.strip()
                if line.startswith("# rest_seq"):
                    rest_seq = line.split()[-1]
                    continue
                chr_, start, end, _, _, _ = parse_line_bed6(line)
                if start == 0:  # remove first
                    continue
                rest_sites.setdefault(chr_, [])
                rest_sites[chr_].append(start)
        for k in rest_sites.keys():
            rest_list = rest_sites[k]
            rest_list.pop()  # remove last
            rest_sites[k] = np.asarray(rest_list)
    else:  # hdf5 file
        with h5py.File(frag_file, 'r') as f:
            rest_seq = f.attrs['rest_seq']
            chromosomes = list(f['chromosomes'].keys())
            for chr_ in chromosomes:
                fragments = f['chromosomes'][chr_][()]
                rest_sites[chr_] = fragments[1:-1]
    return rest_sites, rest_seq


def find_frag(rest_sites, start, end, rest_site_len):
    """
    Assign genome interval to fragment.

    Parameters
    ----------
    rest_sites : `numpy.ndarray`
        restriction sites(positions) of a chromosome.
    start : int
        PET start position
    end : int
        end position
    rest_site_len : int
        Length of restriction site.

    Return
    ------
    n : int
        Index number of fragment.
    pos : {'s', 't'}
        fragment start or end.
    dis : int
        distance to closest cutting site.
    """
    search_s = start + 1  # start search point
    search_e = end - rest_site_len - 1  # end search point
    sites = rest_sites
    frag_idx_s = sites.searchsorted(search_s, side='right')
    frag_idx_e = sites.searchsorted(search_e, side='right')

    if sites.shape[0] == 0:  # no sites in this chromosome
        frag_idx = 0
        pos = 's'
        dis = 0
    elif frag_idx_s <= 0:  # in first fragment
        frag_end = sites[frag_idx_e]
        frag_idx = frag_idx_s
        pos = 't'
        dis = end - frag_end
    elif frag_idx_e >= sites.shape[0]:  # in last fragment
        frag_start = sites[frag_idx_s - 1]
        frag_idx = frag_idx_e
        pos = 's'
        dis = start - frag_start
    else:
        frag_start = sites[frag_idx_s - 1]
        frag_end = sites[frag_idx_e]
        left_span  = start - frag_start
        right_span = end - frag_end
        if abs(left_span) < abs(right_span):
            pos = 's'
            frag_idx = frag_idx_s
            dis = left_span
        else:
            pos = 't'
            frag_idx = frag_idx_e
            dis = right_span

    return (frag_idx, pos, dis)


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


def bedpe_type(rest_sites, bedpe_items, threshold_span, rest_site_len):
    """
    judge interaction type,

    Parameters
    ----------
    rest_sites : dict
        Each chromosome's restriction start positions.
    bedpe_items : list
        bedpe fields.
    threshold_span : int
        If span large than this, will not check.
        Use -1 to force assign PET to fragment.
    rest_site_len : int
        Length of restriction site.

    Return
    ------
    types: {"normal", "self-ligation", "re-ligation"}
    frag1: {(int, {'s', 't'}, int), None}
    frag2: {(int, {'s', 't'}, int), None}
    """
    chr1, start1, end1, chr2, start2, end2 = bedpe_items[:6]
    # inter chromosome interaction
    if chr1 != chr2:
        if threshold_span != -1:
            return "normal", None, None
        else:
            frag1 = find_frag(rest_sites[chr1], start1, end1, rest_site_len)
            frag2 = find_frag(rest_sites[chr2], start2, end2, rest_site_len)
            return "normal", frag1, frag2

    # intra chromosome interaction
    span =  abs(start1 - start2)
    if span > threshold_span and threshold_span != -1:
        # safe span: interaction distance enough long
        if threshold_span != -1:
            return "normal", None, None
    else:
        # unsafe span, need to check
        sites = rest_sites[chr1]
        frag1 = find_frag(sites, start1, end1, rest_site_len)
        frag2 = find_frag(sites, start2, end2, rest_site_len)

        if frag1[0] == frag2[0]:
            return "self-ligation", frag1, frag2
        elif (frag2[0] - frag1[0] == 1) and frag2[1] == 's' and frag1[1] == 't':
            return "re-ligation", frag1, frag2
        else:
            return "normal", frag1, frag2


def worker(task_queue,
           rest_sites, rest_site_len,
           output_file, err_sel_file, err_re_file,
           threshold_span, counts, lock):
    current = mp.current_process().pid
    output_f = open(output_file, 'w')
    err_sel_f = open(err_sel_file, 'w')
    err_re_f = open(err_re_file, 'w')
    l_t_counts = { k: 0 for k in ("normal", "self-ligation", "re-ligation") }
    l_p_counts = { k: 0 for k in [ oc + pc for oc in [i+j for i in "+-" for j in "+-"]
                                           for pc in [i+j for i in "st" for j in "st"] ] }
    l_d_counts = { "PET1": defaultdict(lambda:0), "PET2": defaultdict(lambda:0) }

    while 1:
        try:
            chunk = task_queue.get()
            if chunk is None:
                raise Empty
        except Empty:
            output_f.close()
            err_sel_f.close()
            err_re_f.close()
            lock.acquire()
            for key, l_counts in zip(['type', 'position'], [l_t_counts, l_p_counts]):
                counts_ = counts[key]
                for k, v in l_counts.items():  # update 
                    counts_[k] += v
                counts[key] = counts_  # save to manager dict
            counts_ = counts['distance']
            for key in "PET1", "PET2":
                for k, v in l_d_counts[key].items():
                    counts_[key].setdefault(k, v)
            counts['distance'] = counts_
            lock.release()
            log.debug("Process-%d done."%current)
            break

        for items in chunk:
            try:
                type_, frag1, frag2 = bedpe_type(rest_sites, items, threshold_span, rest_site_len)
            except KeyError as e:
                log.warning(str(e))
                continue

            out_line = "\t".join([str(i) for i in items])
            if frag1:
                out_line = out_line + "\t{}-{}\t{}\t{}-{}\t{}".format(frag1[0], frag1[1], frag1[2], frag2[0], frag2[1], frag2[2])
            out_line += "\n"

            if type_ == 'normal':
                output_f.write(out_line)
                l_t_counts['normal'] += 1
            elif type_ == 'self-ligation':
                err_sel_f.write(out_line)
                l_t_counts['self-ligation'] += 1
            elif type_ == 're-ligation':
                err_re_f.write(out_line)
                l_t_counts['re-ligation'] += 1

            if frag1:
                position = items[8] + items[9] + frag1[1] + frag2[1]
                l_p_counts[position] += 1
                l_d_counts['PET1'][frag1[2]] += 1
                l_d_counts['PET2'][frag2[2]] += 1


def log_counts(counts, log_file):
    log.info("Noise reduce result count:")
    t_counts = counts['type']
    total = sum(list(t_counts.values()))
    n = t_counts['normal']
    s = t_counts['self-ligation']
    r = t_counts['re-ligation']
    n_r = 0 if total == 0 else n / total
    s_r = 0 if total == 0 else s / total
    r_r = 0 if total == 0 else r / total
    msg1 = "normal\t{}\tpercent\t{:.2%}".format(n, n_r)
    msg2 = "self-ligation\t{}\tpercent\t{:.2%}".format(s, s_r)
    msg3 = "re-ligation\t{}\tpercent\t{:.2%}".format(r, r_r)
    msg4 = "total\t{}".format(total)
    log.info("\t" + msg1)
    log.info("\t" + msg2)
    log.info("\t" + msg3)
    log.info("\t" + msg4)
    with open(log_file, 'w') as f:
        f.write("# type count\n")
        f.write("normal\t{}\n".format(n))
        f.write("self-ligation\t{}\n".format(s))
        f.write("re-ligation\t{}\n".format(r))
        f.write("total\t{}\n".format(total))
        f.write("\n")
        f.write("# position count\n")
        for k, v in counts['position'].items():
            f.write("{}\t{}\n".format(k, v))
        f.write("\n")
        f.write("# PET1 distance distribution\n")
        for k, v in sorted(counts['distance']['PET1'].items(), key=lambda t:t[0]):
            f.write("{}\t{}\n".format(k, v))
        f.write("\n")
        f.write("# PET2 distance distribution\n")
        for k, v in sorted(counts['distance']['PET2'].items(), key=lambda t:t[0]):
            f.write("{}\t{}\n".format(k, v))


@click.command(name="noise_reduce")
@click.argument("bedpe")
@click.argument("output")
@click.option("--restriction", "-r",
    required=True,
    help="File which recorded the position of all restriction sites(DNA fragments)"
        "in the genome. BED6 format or hdf5 file(default) is supported."
        "Which can be produced by 'dlohic extract_fragments' command")
@click.option("--processes", "-p", 
    default=1,
    help="Use how many processes to run.")
@click.option("--threshold-span", "-s",
    default=-1,
    show_default=True,
    help="Threshold of pair span. Use -1 to force check.")
@click.option("--log-file",
    default="noise_reduce.log",
    help="Sperate log file for storage some count information.")
def _main(bedpe, output,
         restriction, processes,
         threshold_span,
         log_file):
    """
    Remove noise in DLO-HiC data.

    \b
    There are two kinds of removable abnormal intra chromosomal interacting pair:
        1. self-ligation.
        2. re-ligation

    \b
    In first situation, 
    PET1 and PET2 to will be assigned to same fragment.

    \b
    normal:
                           (frag n)         (frag m, m > n+1)
                           PET1             PET2
                           --------       --------
        genome    <---------------*---..--*--------------->
                                  ^       ^
                                  restriction sites

    \b
    abnormal-1 (self-ligation):
                            (frag n)   (frag n)
                            PET1       PET2
                            --------   --------
        genome    <---------*-------...-------*---------->
                            ^                 ^
                            restriction sites


    \b
    In second situation,
    PET1 and PET2 will be assigned to adjcent fragment.

    \b
    abnormal-2 (re-ligation):
                                      (frag n+1)
                                      PET2
                                      --------
                               --------
                               PET1 (frag n)
        genome    <-------------------*--------------------->


    """

    log.info("noise reduce on file %s"%bedpe)

    log.info("loading restriction sites from file: {}".format(restriction))
    rest_sites, rest = load_rest_sites(restriction)
    _, rest_seq = parse_rest(rest)
    rest_site_len = len(rest_seq)

    # init counts
    lock = mp.Lock()
    counts = mp.Manager().dict()
    counts['type'] = {
        k: 0 for k in ("normal", "self-ligation", "re-ligation")
    }
    counts['position'] = {  # permutation of "+-" and "st"
        k: 0 for k in [
            oc + pc for oc in [i+j for i in "+-" for j in "+-"]
                    for pc in [i+j for i in "st" for j in "st"]
        ]
    }
    counts['distance'] = {
        "PET1": {},
        "PET2": {},
    }

    task_queue = mp.Queue()
    workers = [mp.Process(target=worker,
                          args=(task_queue,
                                rest_sites, rest_site_len,
                                output+".tmp.%d"%i, output+".tmp.sel.%d"%i, output+".tmp.re.%d"%i,
                                threshold_span,
                                counts, lock))
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
    merge_tmp_files(tmp_files, output)
    tmp_sel_files = [output+".tmp.sel.%d"%i for i in range(processes)]
    merge_tmp_files(tmp_sel_files, output+".sel")
    tmp_re_files = [output+".tmp.re.%d"%i for i in range(processes)]
    merge_tmp_files(tmp_re_files, output+".re")

    log_counts(counts, log_file)
    log.info("Noise reduce done.")


main = _main.callback


if __name__ == "__main__":
    eval("_main()")
