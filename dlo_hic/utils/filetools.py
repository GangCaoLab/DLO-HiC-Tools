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


def infer_interaction_file_type(path):
    """
    Inference a interaction file in Pairs or Bedpe format.
    """
    from .parse_text import Bedpe, Pairs, is_comment
    f_st = os.stat(path)
    if f_st.st_size == 0:
        raise IOError("Empty file")
    with open_file(path) as f:
        while True:  # skip header lines
            line = f.readline()
            if not is_comment(line):
                break
        else:
            raise IOError("{} not have any contents.")
    try:  # try Bedpe
        Bedpe(line)
        return Bedpe
    except:
        pass
    try:  # try Pairs
        Pairs(line)
        return Pairs
    except:
        pass
    raise NotImplementedError("Only support pairs and bedpe file format.")

            
