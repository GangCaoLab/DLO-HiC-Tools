import pandas as pd

from .parse_text import Bedpe, Pairs, infer_interaction_file_type, is_comment
from .filetools import open_file


def skip_header(fh):
    while True:  # skip headers
        line = fh.readline()
        if not is_comment(line):
            break
    return fh


def pairs2df(path):
    """
    Construct datafrom from pairs file.
    """
    fields = Pairs.fields
    with open_file(path) as fh:
        fh = skip_header(fh)
        df = pd.read_table(fh, header=None) 
        df.columns = fields
    return df


def bedpe2df(path, pos1="start1", pos2="start2"):
    """
    Construct dataframe from bedpe file.
    """
    fields = Bedpe.fields
    with open_file(path) as fh:
        fh = skip_header(fh)
        df = pd.read_table(fh, header=None) 
        columns = list(fields)
        ex_idx = 1
        while len(df.columns) > columns:  # add extend fields
            columns.append("ex-{}".format(ex_idx))
            ex_idx += 1
        df.columns = columns
        df['pos1'] = df[pos1]  # mimic pairs
        df['pos2'] = df[pos2]
    return df


def interactions_dataframe(path):
    """
    Interaction dataframe.
    """
    fmt = infer_interaction_file_type(path)
    if fmt == Pairs:
        df = pairs2df(path)
    else:
        df = bedpe2df(path)
    return df
