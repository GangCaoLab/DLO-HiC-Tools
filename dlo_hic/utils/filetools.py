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


def merge_tmp_files(tmp_files, merged_file, header=""):
    if os.path.exists(merged_file):
        os.remove(merged_file)
    if header:
        with open(merged_file, 'w') as f:
            f.write(header)
    # merge tmp files
    cmd = "cat " + " ".join(tmp_files) + " >> " + merged_file
    subprocess.check_call(cmd, shell=True)
    # remove tmp files
    cmd = "rm " + " ".join(tmp_files)
    subprocess.check_call(cmd, shell=True)


def is_gziped(path):
    GZ_MAGIC_HEX = "1f8b"
    with open(path, 'rb') as f:
        if f.read(2).hex() == GZ_MAGIC_HEX:
            return True
        else:
            return False
