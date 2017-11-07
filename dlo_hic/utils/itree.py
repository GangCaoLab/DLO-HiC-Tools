import os
import cPickle

from intervaltree import IntervalTree

from .parse_text import parse_line_bed6


def build_bed6_itree(bedfile, idxfile):
    """
    build interval tree for bedfile, and save to disk.
    """
    chr2itree = {}
    with open(bedfile) as f:
        for line in f:
            chr_, start, end, name, score, strand = parse_line_bed6(line)
            chr2itree.setdefault(chr_, IntervalTree())
            itree = chr2itree[chr_]
            itree[start:end] = (name, score, strand)
    with open(idxfile, 'w') as f:
        cPickle.dump(chr2itree, f)
    return chr2itree


def load_bed6_itree(idxfile):
    """
    load bedfile's interval tree object from disk.
    """
    with open(idxfile) as f:
        return cPickle.load(f)