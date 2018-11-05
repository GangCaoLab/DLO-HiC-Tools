"""

stream processing functions,
using the `yield` syntax.

"""

import logging

log = logging.getLogger(__name__)


def beds2bedpe(bed1_path, bed2_path):
    """
    merge two ends bed file to bedpe.
    """
    import subprocess as subp
    cmd = "join -j 4 {} {}".format(bed1_path, bed2_path)
    p = subp.Popen(cmd, shell=True, stdout=subp.PIPE)
    #
    # format joined bed bedpe
    for line in p.stdout:
        line = line.decode("utf-8")
        items = line.strip().split()
        name, chr_a, s_a, e_a, score_a, strand_a, \
        chr_b, s_b, e_b, score_b, strand_b = items
        outitems = [chr_a, s_a, e_a,
                    chr_b, s_b, e_b,
                    name, '0', strand_a, strand_b]
        outline = "\t".join(outitems)
        yield outline


def upper_triangle(line_iterator, fmt='bedpe'):
    """
    transform bedpe file's all line to upper trangle form.

    Arguments
    ---------
    fmt : str
        The input line format, 'bedpe' or 'pairs'
    """
    from dlo_hic.utils.parse_text import Bedpe, Pairs
    itr = line_iterator
    for line in itr:
        if fmt == 'bedpe':
            o = Bedpe(line)
            o.to_upper_trangle()
        elif fmt == 'pairs':
            o = Pairs(line)
            o.to_upper_trangle()
        else:
            raise ValueError("fmt only 'bedpe' or 'pairs'.")
        outline = str(o)
        yield outline


def write_to_file(line_iterator, output_path):
    itr = line_iterator
    with open(output_path, 'w') as of:
        for line in itr:
            of.write(line + "\n")
