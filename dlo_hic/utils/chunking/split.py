from subprocess import Popen, PIPE
from os.path import splitext, dirname, basename
import os

from dlo_hic.utils.filetools import is_gziped


def split_by_n_chunk(fq_path, n):
    if is_gziped(fq_path):
        raise IOError("Input file is gziped, please gunzip it firstly.")

    cmd = ['wc', '-l', fq_path]
    p = Popen(cmd, stdout=PIPE)
    out, _ = p.communicate()
    p.kill()
    nlines = int(out.decode('utf-8').split()[0])
    if nlines % 4 != 0:
        raise IOError("Truncated file. Numer of lines: {}".format(nlines))

    n_recs = nlines // 4
    recs_per_file = (n_recs + n - 1) // n
    lines_per_file = recs_per_file * 4
    last_file_recs = n_recs - (recs_per_file * (n-1))

    prefix, suffix = splitext(fq_path)
    cmd = ['split', fq_path, '-l', str(lines_per_file), prefix+".chunk.", "--additional-suffix="+suffix]
    p = Popen(cmd)
    p.wait()
    p.kill()

    dir_ = dirname(fq_path)
    splited_files = sorted([f for f in os.listdir(dir_) if basename(prefix)+".chunk" in f])
    assert len(splited_files) == n

    return splited_files

