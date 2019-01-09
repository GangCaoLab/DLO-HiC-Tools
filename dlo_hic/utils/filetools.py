import gzip
import os
import io
import subprocess


def sample_file(path, output, n):
    """
    Sample random lines from a file.

    Parameters
    ----------
    path : str
        Path to text or gziped text file.
    output : str
        Path to output file.
    n : int
        How many lines in sampled file.
    """
    if path.endswith(".gz"):
        cmd = "zcat {} | shuf -n {} > {}".format(path, n, output)
    else:
        cmd = "shuf -n {} {} > {}".format(n, path, output)
    subprocess.check_call(cmd, shell=True)


def open_file(fname, mode='r'):
    if fname.endswith('.gz'):
        fh = gzip.open(fname, mode=mode+'b')
        fh = io.TextIOWrapper(fh)
    else:
        fh = open(fname, mode=mode)
    return fh

